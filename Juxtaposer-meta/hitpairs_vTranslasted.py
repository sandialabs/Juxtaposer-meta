#!/usr/bin/env python3

import sys
import os
import re
from collections import defaultdict, Counter
import subprocess

# Usage check
if len(sys.argv) < 5:
	print(f"Usage: python {sys.argv[0]} outputDirectory referencePathPrefix maxReadLength palflag(0 or 1)")
	sys.exit(1)

outDir, refdir, readmax_str, palflag_str = sys.argv[1:5]
readmax = int(readmax_str)
palflag = int(palflag_str)

delta = '^'
xslen = 20000000
endwindow = 5
uniqlen = 8
spanlen = readmax * 2

# define data structures
dnas,seqs,mobs,tns,juxtas = {},{},{},{},{}
mobloc = defaultdict(dict)
uids = defaultdict(int)
cts = defaultdict(int)
cts['overlap']=defaultdict(int)

# defline and check required files
lens_file,mobile_bed,fasta_file,blast_file = [f"{refdir}.lens", f"{refdir}.mobile.bed", os.path.join(outDir, "candidates.fa"),os.path.join(outDir, "candidates.blast")]
for file in [lens_file,mobile_bed,fasta_file,blast_file]:
		if not os.path.isfile(file):
				sys.exit(f"Could not open {file}")
tn_bed = f"{refdir}.tn.bed" #this file is not required

def readfile(file,obj_type=open):
	with obj_type(file,'r') as f:
		for line in f:
			yield line.strip()			
def revcomp(seq):
	complement = str.maketrans('ACGTacgt', 'TGCAtgca')
	return seq.translate(complement)[::-1]

def register(pos, length):
	pos = pos % length
	if pos < 1:
		pos += length
	return pos

def load_ref_files():
	for line in readfile(lens_file):
		line = line.rstrip('\n')
		m = re.match(r'^(\S+).*\t(\d+)\s*', line)
		if m:
			key = m.group(1)
			dnas[key] = {'linear': line, 'len': int(m.group(2))} # Store the entire line as 'linear' and length as int
	# Load read sequences from fasta
	fa_obj=readfile(fasta_file)
	seq=[]
	defline=next(fa_obj).lstrip('>').split()[0]
	while True:
		try:
			line=next(fa_obj)
			if not line.startswith('>'):
				if seq:
					seq.append(line)
				else:
					seq=[line]
			else:
				seqs[defline]=''.join(seq)
				seq,defline=[],line.lstrip('>').split()[0]
		except StopIteration: # handle last seq
			seqs[defline]=''.join(seq)
			break
	mobSeg=defaultdict(list)
	for line in readfile(mobile_bed):
		flds = line.split('\t')
		if len(flds) < 4:
			continue
		dna, L, R, mob_id = flds[0], int(flds[1]), int(flds[2]), flds[3]
		mobs[mob_id] = {'dna': dna, 'L': L+1, 'R': R}
		if not re.search(r'Integrase',mob_id,re.I):
			continue
		mobloc[dna][L+1]=mob_id
		mobloc[dna][R]=mob_id
		if L+1 <= 5746:
			contigLeft,remainsLeft=1,abs(L+1-5746)
		else:
			contigLeft,remainsLeft=L+1-5746,None
		if R+5746 > dnas.get(dna,{}).get('len'):
			contigRight,remainsRight=dnas.get(dna,{}).get('len'),5746-(dnas.get(dna,{}).get('len')-R)
		else:
			contigRight,remainsRight=R+5746,None
		mobSeg[dna].append({'contigLeft':contigLeft,'remainsLeft':remainsLeft,'contigRight':contigRight,
		'remainsRight':remainsRight,'id':mob_id,'L': L+1, 'R': R, 'midpoint':int((L+1+R)/2)})
	if os.path.isfile(tn_bed):
		for line in readfile(tn_bed):
			line = line.rstrip('\n')
			flds = line.split('\t')
			if len(flds) < 4:
				continue
			dna, L, R, tn_id = flds[0], int(flds[1]), int(flds[2]), flds[3]
			uids[tn_id] += 1
			count = uids[tn_id]
			tns[f"{tn_id}.{count}L"] = {'dna': dna, 'L': L+1, 'R': L+5, 'term': L+1, 'end': -1}
			tns[f"{tn_id}.{count}R"] = {'dna': dna, 'L': R-5+1, 'R': R, 'term': R, 'end': 1}
	return


# Load hits into reads
def load_read(read, temphits):
	recordreads = False
	for hit in temphits:
		dna = hit['dna']
		if dna not in mobloc:
			continue
		for intcoord in mobloc[dna].keys():
			if abs(hit['sL'] - intcoord) <= 5746 or abs(hit['sR'] - intcoord) <= 5746:
				recordreads = True
				break
		if recordreads:
			break
	if recordreads:
		process_reads(read,temphits)

def per_read():
	read = None
	ReadHits = []
	for line in readfile(blast_file):
		line = line.rstrip('\n')
		flds = line.split('\t')
		if len(flds) < 13:
			continue
		current_read = flds[0]
		if read is not None and current_read != read:
			# Check seq defined in hits
			load_read(read, ReadHits)
			ReadHits = []
		read = current_read
		dna, id_, btop = [flds[i] for i in [1,2,12]] #string fields
		qL, qR, sS, sE = [int(flds[i]) for i in [6,7,8,9]] #all numeric fields
		sco = float(flds[11]) #sometimes a float
		dir_ = 1
		clos = 0
		sL, sR = sS, sE
		if sL > sR:
			sL, sR = sE, sS
			dir_ = -1

		# Process btop string
#		if not re.match(r'^\d{0,}[\-ACTG]{0,2}\d{0,}$', btop): # one or fewer snps
		if not re.match(r'^\d+$', btop): # must be perfect
			matches = re.findall(r'(\d+|-[ACGT]|[ACGT]-|[ACGT]{2})', btop)
			shift_start = 0
			max_length = 0
			new_btop = 0
			current_length = 0
			current_position = 0
			dels = []
			ins = []
			for i, match in enumerate(matches):
				if re.match(r'^\d+$', match):
					current_length = int(match)
					if current_length > max_length:
						max_length = current_length
						shift_start = current_position
						new_btop = i
					current_position += current_length
				else:
					current_length = 0
					current_position += 1
					if re.match(r'-[ACGT]', match):
						dels.append(i)
					if re.match(r'[ACGT]-', match):
						ins.append(i)
			if current_length > max_length:
				max_length = current_length
				shift_start = current_position
				new_btop = len(matches) - 1
			delcount = sum(1 for d in dels if d < new_btop)
			incount = sum(1 for i_ in ins if i_ < new_btop)
			bqL, bqR = qL, qR
			if shift_start:
				qL = qL + (shift_start - delcount)
				sL = sL + (shift_start - incount)
			qR = qL + max_length - 1
			sR = sL + max_length - 1
			if dir_ == -1:
				sS, sE = sR, sL
			else:
				sS, sE = sL, sR
			if qR > 150:
				print(f"before:{bqL} {bqR} {btop} shift: {shift_start} extends: {max_length} insertions: {incount} deletions: {delcount} after:{qL} {qR} {sL} {sR} {matches[new_btop]}")
			btop = matches[new_btop]

		if qL == 1 and qR >= readmax:
			continue

		# Handle closure hits

		fullseq = seqs.get(read, '')
		length = qR - qL + 1
		offset = qL - 1
		subseqs = []
		if offset != 0:
			subseqs = [fullseq[0:offset], fullseq[offset:offset+length], fullseq[qR:]]
		else:
			subseqs = ['', fullseq[offset:offset+length], fullseq[qR:]]

		seq = ''
		for i in range(3):
			if i == 1 and len(subseqs[i]) > 0:
				seq += subseqs[i].upper()
			elif len(subseqs[i]) > 0:
				seq += subseqs[i].lower()
		ReadHits.append({'dna': dna, 'qL': qL,'qR': qR,'sL': sL,'sR': sR, 'id': id_,'len': length,'seq': seq,
			'sco': sco,'sS': sS,'sE': sE,'dir': dir_,'offset': offset,'btop': btop,'clos': clos})
	# After loop ends, load last read hits
	if read is not None and ReadHits:
		for hit in ReadHits:
			if 'seq' not in hit or not hit['seq']:
				print(f"Undefined seq in hit for read {read}\n{hit}")
				sys.exit(1)
		load_read(read, ReadHits)
	cts['nonstandard']=len(seqs)
	cts['mobs']=len(mobs)
	cts['tn']=len(tns)

def process_reads(read,hits):
	#for read, hits in reads.items():
	print(len(seqs))
	cts['genomic'] += 1
	hitct = len(hits)
	if hitct <= 1:
		return #continue
	else:
		cts['multi'] += 1
		if hitct == 2:
			cts['double'] += 1
		checklog=[]
		for i in range(len(hits)):
			for j in range(len(hits)):
				if i==j:
					continue
				elif sorted([i,j]) in checklog:
					continue
				checklog.append(sorted([i,j]))
		for n in checklog:
			print(read,n)
		# Find top scoring hits
		tophit = [0]  # [0] top score, [1..] indices of hits with that score
		hitdict=defaultdict(list)
		for i, hit in enumerate(hits):
			hitdict[hit['sco']].append(i)
		tophit=sorted(list(hitdict.items()),key= lambda s:float(s[0]),reverse=True)[0] #tophit==(top score(float),[first index of hit in hits with top score,...,])
		hitpairs = []
		# Remove leading spaces from tophit[0] if string (should be float, but keep for compatibility)
		if isinstance(tophit[0], str):
			tophit[0] = tophit[0].lstrip()
		# Sort tophit indices by dna(lexical) and sL if dna is the same
		sorted_tophit_indices = sorted(tophit[1], key=lambda i: (hits[i]['dna'], hits[i]['sL']))
		for idx_i in sorted_tophit_indices: # for all tophits
			for j in range(len(hits)): #for all remaining hits
				if j == idx_i:
					continue
				# Sort pair by qL
				gnmOrder = sorted([idx_i, j], key=lambda x: hits[x]['qL'])
				pair = [hits[gnmOrder[0]], hits[gnmOrder[1]]]
				p0, p1 = pair[0], pair[1]
				if p0['qL'] > (p1['qL'] - uniqlen) or p0['qR'] > (p1['qR'] - uniqlen): # Reject non-shifted configuration
					#print(f"REJECTED:{i},{j} of {read}, no shift. LEFT0,LEFT1:", p0['qL'],p1['qL'],"RIGHT0,RIGHT1", p0['qR'],p1['qR'])
					#print(p0['dna'],p0['dir'],p0['seq'])
					#print(p1['dna'],p1['dir'],p1['seq'])
					continue
				sco = p0['sco'] + p1['sco']
				if hitpairs and sco < hitpairs[0]: # Reject lower score than previous
					continue
				overlap = p0['qR'] - p1['qL'] + 1 # length of sequence shared between pair
				goff0 = roff0 = goff1 = roff1 = 0
				if overlap > 0:
					#print(f"{i},{j} of {read}, positive {overlap}.")
					# Adjust for overlap imperfections using btop strings
					btop0 = p0['btop']
					last = 0
					if not re.search(r'^[0-9]+$', btop0):
						#print(f"START btop0 parsing.dir:{p0['dir']}, sE: {p0['sE']},qR: {p0['qR']}")
						while btop0:
							m = re.search(r'([0-9]+)$', btop0)
							if m:
								#print(m)
								last = int(m.group(1))
								#print(f"btop0,-hit:{btop0}, {-len(m.group(1))}")
								btop0 = btop0[: -len(m.group(1))]
								#print(f"btop0[: -len(m.group(1)):{btop0}")
								if last + roff0 > overlap:
									break
							if len(btop0) < 2:
								break
							pair_chars = btop0[-2:]
							#print(f"pair_chars:{pair_chars}")
							btop0 = btop0[:-2]
							#print(f"btop0[:-2]:{btop0}")
							roff0 += last + 1
							goff0 += last + 1
							last = 0
							#print(f"roff,goff:{roff0},{goff0}")
							if pair_chars[0] == '-':
								roff0 -= 1
							elif pair_chars[1] == '-':
								goff0 -= 1
							#print(f"roff,goff:{roff0},{goff0}")
							if roff0 > overlap:
								break
						p0['sE'] -= goff0 * p0['dir']
						p0['qR'] -= roff0
						#print(f"END btop0 parsing.dir:{p0['dir']}, sE: {p0['sE']},qR: {p0['qR']}")

					btop1 = p1['btop']
					last = 0

					if not re.search(r'^[0-9]+$', btop1):
						#print(f"START btop1 parsing.dir:{p1['dir']}, sS: {p1['sS']},qL: {p1['qL']}")
						while btop1:
							m = re.match(r'^([0-9]+)', btop1)
							if m:
								#print(m)
								last = int(m.group(1))
								#print(f"btop1,-hit:{btop1}, {-len(m.group(1))}")
								btop1 = btop1[len(m.group(1)):]
							#	print(f"btop1[: -len(m.group(1)):{btop1}")
								if last + roff1 > overlap:
									break
							if len(btop1) < 2:
								break
							pair_chars = btop1[:2]
						#	print(f"pair_chars:{pair_chars}")
							btop1 = btop1[2:]
						#	print(f"btop1[:-2]:{btop1}")
							roff1 += last + 1
							goff1 += last + 1
							last = 0
						#	print(f"roff,goff:{roff1},{goff1}")
							if pair_chars[0] == '-':
								roff1 -= 1
							elif pair_chars[1] == '-':
								goff1 -= 1
						#	print(f"roff,goff:{roff0},{goff0}")
							if roff1 > overlap:
								break
						p1['sS'] -= goff1 * p1['dir']
						p1['qL'] += roff1
						#print(f"END btop1 parsing.dir:{p1['dir']}, sS: {p1['sS']},qL: {p1['qL']}")
					overlap -= (goff0 + goff1)
				#print("before hit_idx loop")
				for hit_idx in [0, 1]:
					dna = pair[hit_idx]['dna']
					for key in ['sE', 'sS', 'sL', 'sR']:
						if dna not in dnas or 'len' not in dnas[dna]:
							raise ValueError(f"DNA length missing for {dna}")
						pair[hit_idx][key] = register(pair[hit_idx][key], dnas[dna]['len'])
				shift = [1, pair[0]['sE'], pair[1]['sS'], overlap, 0, -1]
				comp = ( (pair[0]['dna'] > pair[1]['dna']) - (pair[0]['dna'] < pair[1]['dna']),
						 shift[1] - shift[2] )
				if comp[0] == 1 or (comp[0] == 0 and comp[1] == 1):
					# Reverse order
					gnmOrder.reverse()
					pair.reverse()
					p0, p1 = pair
					shift = [-1, shift[2], shift[1], overlap, 0, -1]

				if comp[0] == 0:
					shift[4] = 1
					shift[5] = shift[2] - shift[1]
					otherway = dnas[p0['dna']]['len'] - shift[5]
					if otherway < shift[5]:
						shift[4] = -1
						shift[5] = otherway

				config = ['+', '+', '', '']
				if p0['dir'] * shift[0] == -1:
					config[0] = '-'
				if -1 * p1['dir'] * shift[0] == -1:
					config[1] = '-'
				config[2] = config[0] + config[1]

				# Convert '+'/'-' to 1/-1 for config[0] and config[1]
				config[0] = 1 if config[0] == '+' else -1
				config[1] = 1 if config[1] == '+' else -1

				tnend = [[], []]
				for hit_idx in [0, 1]:
					for tn_id, tn_info in tns.items():
						hit = pair[hit_idx]
						if (hit['sR'] < tn_info['L'] or hit['sL'] > tn_info['R'] or
							tn_info['dna'] != hit['dna'] or config[hit_idx] * tn_info['end'] == -1):
							continue
						tnend[hit_idx].append([tn_id, config[hit_idx] * (shift[hit_idx + 1] - tn_info['term'])])
					if tnend[hit_idx]:
						tnend[hit_idx].sort(key=lambda x: abs(x[1]))
						tnend[hit_idx].extend(tnend[hit_idx][0])
					else:
						tnend[hit_idx].extend(['.', '.'])

				if tnend[1][-1] != '.' and tnend[1][-1] > 1000:
					raise ValueError(f"{tnend[1][-2]} {tnend[1][-1]} {shift[2]} {config[1]} {tns[tnend[1][-2]]['term']}")

				seq = p0['seq']
				length = xslen
				align = []
				for hit in pair:
					hitseq = list(hit['seq'][hit['offset']:hit['offset'] + hit['len']])
					btop,offset = hit['btop'],0
					while btop:
						m = re.match(r'^(\d+)', btop)
						if m:
							offset += int(m.group(1))
							btop = btop[m.end():]
							continue
						if len(btop) < 2:
							break
						pair_chars = btop[:2]
						btop = btop[2:]
						if pair_chars[0] == '-':
							hitseq.insert(offset, '-')
						elif pair_chars[1] == '-':
							hitseq[offset] = delta
						else:
							hitseq[offset] = hitseq[offset].lower()
						offset += 1
					align.append(''.join(hitseq))
				if p0['dna'] == p1['dna']:
					if (config[2] == '-+' and shift[4] == 1) or (config[2] == '+-' and shift[4] == -1):
						config[3] = 'circleJxn'
					elif (config[2] == '-+' and shift[4] == -1) or (config[2] == '+-' and shift[4] == 1):
						config[3] = 'scar'
					elif (config[2] in ['--', '++']) and shift[5] < 10000:
						config[3] = 'palindrome'
					length = shift[5]

				if palflag and config[3] == 'palindrome':
					continue

				mobct = sum(1 for h in [0, 1] if tnend[h][-2] == '.')
				if not config[3] and mobct == 1:
					config[3] = 'transpose'
				if not config[3]:
					config[3] = '.'

				cts[config[3]] += 1

				second = config[2][1]
				second = '+' if second == '-' else '-'
				config[2] = config[2][0] + second

				order = sorted([p0['offset'], p0['offset'] + p0['len'], p1['offset'], p1['offset'] + p1['len']])
				jL, jR = order[1], order[2]

				jxn = seq[jL - 5:jL] + '/' + seq[jL:jR] + '/' + seq[jR:jR + 5]

				seq = seq[:p1['qL'] - 1] + '-' + seq[p1['qL'] - 1:p1['qR']] + '-' + seq[p1['qR']:]

				if shift[0] == -1:
					seq = revcomp(seq)
					jxn = revcomp(jxn)
					align = [revcomp(a) for a in align]

				ct = 0
				if hitpairs and sco == hitpairs[0]:
					ct = hitpairs[-1] + 1

				if not hitpairs or sco > hitpairs[0] or length < hitpairs[1]:
					hitpairs = [sco, length, shift[1], shift[2], config[2], seq,
								tnend[0][-2], tnend[0][-1], tnend[1][-2], tnend[1][-1],
								shift[3], align[0], align[1], jxn, config[3], shift[4],
								gnmOrder[0], gnmOrder[1], 0]
				else:
					hitpairs[-1] += 1


		if not hitpairs:
			return #continue
		cts['shift'] += 1
		cts['overlap'][hitpairs[10]]+= 1

		if hitpairs[1] == xslen:
			hitpairs[1] = -1

		hit1 = hits[hitpairs[-3]]
		hit2 = hits[hitpairs[-2]]
		label = f"{hit1['dna']}/{hitpairs[2]}/{hit2['dna']}/{hitpairs[3]}/{hitpairs[4]}"
		if label not in juxtas or juxtas[label]['sco'] < hitpairs[0]:
			ct = 1
			if label in juxtas:
				ct += juxtas[label]['ct']
			juxtas[label] = {
				'hit1dna':hit1['dna'],'hit2dna':hit2['dna'],'hit1id':hit1['id'],'hit2id':hit2['id'],
				'read': read, 'ct': ct,
				'hit1': hitpairs[-3],
				'hit2': hitpairs[-2],
				'readseq': hitpairs[5],
				'sites': hitpairs[-1],
				'coord': hitpairs[2],
				'coordp': hitpairs[3],
				'config': hitpairs[4],
				'q1R': hitpairs[7],
				'q2L': hitpairs[8],
				'overlap': hitpairs[10],
				'sco': hitpairs[0],
				'mob1': hitpairs[6],
				'mob1off': hitpairs[7],
				'mob2': hitpairs[8],
				'mob2off': hitpairs[9],
				'len': hitpairs[1],
				'hitseq1': hitpairs[11],
				'hitseq2': hitpairs[12],
				'jxn': hitpairs[13],
				'call': hitpairs[14],
				'origin': hitpairs[15],
			}
		else:
			juxtas[label]['ct'] += 1
def write_juxtas_bed():
	bed_file = os.path.join(outDir, "juxtas.bed")
	with open(bed_file, 'w') as out:
		overlaps = {}
		for label, line in juxtas.items():
			#hit1 = reads[line['read']][line['hit1']]
			dna = line['hit1dna']
			L = line['coord'] - 1
			R = line['coordp']
			if line['origin'] == -1:
				out.write(f"{dna}\t{L}\t{dnas[dna]['len']}\t{label}\n")
				out.write(f"{dna}\t0\t{R}\t{label}\n")
			else:
				if R < L:
					R = L
				out.write(f"{dna}\t{L}\t{R}\t{label}\n")
	return bed_file

def run_bedtools_intersect(bed_file):
	cmd = ['bedtools', 'intersect', '-wo', '-a', bed_file, '-b', f"{refdir}.mobile.bed"]
	result = subprocess.run(cmd, capture_output=True, text=True)
	if result.returncode != 0:
		raise RuntimeError(f"bedtools intersect failed: {result.stderr}")
	return result.stdout.splitlines()

def process_intersections(lines):
	overlaps = {}
	for line in lines:
		flds = line.split('\t')
		label = flds[3]
		mobile = flds[7]
		Lj, Rj = int(flds[1]), int(flds[2])
		Lm, Rm = int(flds[5]), int(flds[6])
		dL = Lm - Lj
		dR = Rj - Rm
		pct = 100
		lenM = Rm - Lm
		low = Lm - Lj
		if low > dR:
			low = dR
		if dL < 0:
			pct += 100 * dL / lenM
		if dR < 0:
			pct += 100 * dR / lenM
		if low < 0:
			if label not in overlaps or overlaps[label][1] < pct:
				overlaps[label] = [0, pct, mobile]
		else:
			if label not in overlaps or low < overlaps[label][0]:
				overlaps[label] = [low, pct, mobile]
	return overlaps

def write_juxtas_txt(overlaps):
	txt_file = os.path.join(outDir, "juxtasint_py.txt")
	with open(txt_file, 'w') as out:
		header = ['dnaLo', 'CoordLo', 'dnaHi', 'CoordHi', 'config', 'call', 'overlap', 'origin',
				  'idLo', 'idHi', 'tnLo', 'tnOffsetLo', 'tnHi', 'tnOffsetHi', 'sampleRead',
				  'overlapSeq', 'readseq', 'seqLo', 'seqHi', 'count', 'sites', 'Mob', 'pctMob', 'shortDistMob']
		out.write('#' + '\t'.join(header) + '\n')
		for label in sorted(juxtas.keys(), key=lambda x: juxtas[x]['ct'], reverse=True):
			line = juxtas[label]
			#hit1 = reads[line['read']][line['hit1']]
			#hit2 = reads[line['read']][line['hit2']]
			mobile = ''
			if label in overlaps:
				mobile = '\t'.join(map(str, reversed(overlaps[label])))
			fields = [
				line['hit1dna'], line['coord'], line['hit2dna'], line['coordp'], line['config'], line['call'],
				line['overlap'], line['origin'], line['hit1id'], line['hit2id'], line['mob1'], line['mob1off'],
				line['mob2'], line['mob2off'], line['read'], line['jxn'], line['readseq'], line['hitseq1'],
				line['hitseq2'], line['ct'], line['sites'], mobile
			]
			out.write('\t'.join(map(str, fields)) + '\n')
		out.write("# ok\n")

def print_counts():
	for key in ['nonstandard', 'genomic', 'multi', 'double', 'shift', 'circleJxn', 'scar', 'palindrome', 'transpose', 'tn']:
		#if key not in cts:
		#	cts[key] = 0
		print(f"{key}={cts[key]}", end=' ')
	print()
	if 'overlap' in cts:
		for k in sorted(list(cts['overlap'].keys())):
			print(f"overlap {k}={cts['overlap'][k]}")

# Main execution
if __name__ == "__main__":
	load_ref_files()
	per_read()
	#process_reads()
	#print("JUXTAS",'\n'.join([f"{key}\t{juxtas[key]}" for key in juxtas]))
	#sys.exit()
	bed_file = write_juxtas_bed()
	intersect_lines = run_bedtools_intersect(bed_file)
	os.unlink(bed_file)
	overlaps = process_intersections(intersect_lines)
	write_juxtas_txt(overlaps)
	print_counts()





