#!/usr/bin/env python3

import sys, os, re, pickle
from collections import defaultdict

# Usage check
if len(sys.argv) < 5:
	print(f"Usage: python {sys.argv[0]} outputDirectory referencePathPrefix maxReadLength palflag(0 or 1)")
	sys.exit(1)

outDir, refdir, readmax_str, palflag_str = sys.argv[1:5]
intfilter=False
if len(sys.argv)==6:
	intfilter=True
runcyc=2
readmax = int(readmax_str)
palflag = int(palflag_str)

xslen = 20000000
uniqlen = 8 # length of distally unshared sequence required for a read to pass
min_cont_len=500
mo_patterns=[re.compile(pattern) for pattern in [r'([0-9]+)$',r'^([0-9]+)']]
oriswitch={"+":"-","-":"+",1:-1,-1:1,"L":"R","R":"L","B":"E","E":"B"}

# define data structures
dnas,seqs,tns,juxtas,identical,genomebins = {},{},{},{},{},{}
mobSeg=defaultdict(list)
cts,uids = defaultdict(int),defaultdict(int)
cts['overlap']=defaultdict(int)
cycle=0
# defline and check required files
lens_file,mobile_bed,fasta_file,blast_file = [f"{refdir}.lens", f"{refdir}.mobile.bed", os.path.join(outDir, "candidates.fa"),os.path.join(outDir, "candidates.blast")] 
for file in [lens_file,mobile_bed,fasta_file,blast_file]:
	if not os.path.isfile(file):
		sys.exit(f"Could not open {file}")
juxtasout=os.path.join(outDir,"juxtas.pkl")
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

def juxtalabel2dict(label):
	dnaLo,coordLo,dnaHi,coordHi,config=label.split('/')
	return {dnaLo:int(coordLo),dnaHi:int(coordHi)}

def load_ref_files():
	longestcontig=0
	for line in readfile(lens_file):
		line = line.rstrip('\n')
		m = re.match(r'^(\S+).*\t(\d+)\s*', line)
		if m:
			key = m.group(1)
			conlen=int(m.group(2))
			dnas[key] = {'linear': line, 'len': conlen} # Store the entire line as 'linear' and length as int
			if conlen>longestcontig:
				longestcontig=conlen
	bin=1
	for i in range(1,longestcontig+1):
		if i % readmax ==0:
			bin+=1
		genomebins[str(i)]=bin
	# Load read sequences from fasta
	fa_obj=readfile(fasta_file)
	seq=[]
	try:
		defline=next(fa_obj).lstrip('>').split()
	except StopIteration:
		print(f"{fasta_file} is empty")
		sys.exit()
	while True:
		try:
			line=next(fa_obj)
			if not line.startswith('>'):
				if seq:
					seq.append(line)
				else:
					seq=[line]
			else:
				if len(defline)>1:
					identical[defline[0]]=len(defline[1::])
				seqs[defline[0]]=''.join(seq)
				seq,defline=[],line.lstrip('>').split()
		except StopIteration: # handle last seq
			seqs[defline[0]]=''.join(seq)
			break
	for line in readfile(mobile_bed):
		flds = line.split('\t')
		if len(flds) < 4:
			continue
		dna, L, R, mob_id, dir= flds[0], int(flds[1]), int(flds[2]), flds[3], flds[5]
		if not re.search(r'PhageIntegrase|Recombinase',mob_id,re.I): 
			continue
		cts['mobs']+=1
		mobSeg[dna].append({'dna':dna,'id':mob_id,'L': L+1, 'R': R, 'dir':dir, 'midpoint':int((L+1+R)/2),'history':[]})
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
	cts['nonstandard']=len(seqs)
	cts['tn']=len(tns)
	return

def findnearbyInt(dna,left,right):
	'''requires the dna and two coordinates'''
	hitmp=int((left+right)/2)
	mobgenes=[]
	thresh=3801 if intfilter else 3801 if cycle==0 else readmax
	if dna in mobSeg:
		try:
			for pos,intgene in sorted(enumerate(mobSeg[dna]), key=lambda x: abs(x[1]['midpoint']-hitmp)):
				intToTheRight = intgene["midpoint"] > hitmp # False==0, True==1
				gapsize=[thresh-(right-intgene['R']), thresh-(intgene['L']-left)][(intToTheRight)] #the more positive the gap, the closer the int gene is
				if gapsize>0:
					if intgene['history']:
						for links in intgene['history']: # this means the putative site is secondary
							for link in links:
								 #dnaLo/coordLo/dnaHi/coordHi -> {dnaLo:coordLo,dnaHi:coordHi}
								juxpos=juxtalabel2dict(link).get(dna)
								distance=min([abs(juxpos-coord) for coord in [left,right]])
								if distance<=readmax:
									mobgenes.append((dna,pos,intToTheRight,gapsize))
								elif dnas.get(dna) and min([left,right,juxpos])<=150:
									proxtermd=min([abs(readmax-i) for i in [left,right,juxpos]])
									dlen=dnas.get(dna)["len"]
									disttermd=min([abs(dlen-i) for i in [left,right,juxpos]])
									if (proxtermd+disttermd)<=readmax:
										mobgenes.append((dna,pos,intToTheRight,gapsize))

					else:
						mobgenes.append((dna,pos,intToTheRight,gapsize))



		except Exception as e:
			print(e,mobSeg[dna])
			sys.exit()
			return

	return mobgenes

# Load hits into reads
def preprocess(hits):
	''' if there are many perfect or identical hits from the same read, then we pass'''
	left,right=True,True
	for side in "L R".split():
		errors={}
		freq=0
		for i,poshit in hits[side].items():
			for hit in poshit:
				error=re.findall("[A-Z]+",hit['btop'])
				if error:
					errors[hit['dna']]=error
				else:
					freq+=1
		if freq and errors:
			pass
		elif freq>3:
			if side=="L":
				left=False
			else:
				right=False
		else:
			uniqer=set([tuple(val) for key,val in errors.items()])
			if len(errors)>3 and len(uniqer)==1:
				if side=="L":
					left=False
				else:
					right=False
	return left and right
def per_read():
	read = None
	ReadHits = []
	orderedhits={'L':defaultdict(list),'R':defaultdict(list)}
	perreadct=0
	print("Running per_read")
	for line in readfile(blast_file):
		line = line.rstrip('\n')
		flds = line.split('\t')
		if len(flds) < 13:
			continue
		current_read = flds[0]
		if read is not None and current_read != read:
			if perreadct>1:
				prepass=preprocess(orderedhits)
				if prepass:
					process_read(read, orderedhits)
			else:
				cts['genomic']+=1
			orderedhits={'L':defaultdict(list),'R':defaultdict(list)}
			perreadct=0
		read = current_read
		perreadct+=1
		dna, id_, btop = [flds[i] for i in [1,2,12]] #string fields
		if dnas[dna]['len'] < min_cont_len:
			continue
		qL, qR, sS, sE = [int(flds[i]) for i in [6,7,8,9]] #all numeric fields
		sco = float(flds[11]) #sometimes a float
		dir_ = 1
		clos = 0
		sL, sR = sS, sE
		if sL > sR:
			sL, sR = sE, sS
			dir_ = -1
		fullseq = seqs.get(read, '')
		if not fullseq:
			print(f"Could not retrieve sequence for {read}, not in fasta file.")
			continue
		readlen=len(fullseq)
		left = qL - 1 < readlen - qR

		if qL == 1 and qR >= readlen:
			continue


		length = qR - qL + 1
		offset = qL - 1
		subseqs = [fullseq[0:offset].lower() if offset!= 0 else '', fullseq[offset:offset+length].upper(), fullseq[qR:].lower() if qR<readmax else '']
		seq= ''.join(subseqs)
		orderedhits['L' if left else 'R'][qL-1 if left else readlen-qR].append({'dna': dna, 'qL': qL,'qR': qR,'sL': sL,'sR': sR, 'id': float(id_),'len': length,'seq': seq,
			'sco': sco,'sS': sS,'sE': sE,'dir': dir_,'offset': offset,'btop': btop,'clos': clos})

	if read is not None and (orderedhits['L'] or orderedhits['R']):
		if  perreadct>1:
			prepass=preprocess(orderedhits)
			if prepass:
				process_read(read, orderedhits)
		else:
			cts['genomic']+=1


def process_read(read,hits):

	weight= 1 + identical.get(read,0) #number of reads with identical sequence as current read
	hitpairs = {}
	hitct = len(hits['L'])+ len(hits['R'])
	cts['genomic'] += weight
	if hitct > 1:
		cts['multi'] += weight
	if hitct == 2:
		cts['double'] += weight
	if hitct <= 1:
		return
	bestleft=sorted(list(hits['L'].items()), key= lambda x: x[0])
	bestright=sorted(list(hits['R'].items()),key= lambda x: x[0])
	for ld in bestleft:
		for l_id,lhit in enumerate(ld[1]):
			for rd in bestright:
				for r_id,rhit in enumerate(rd[1]):
					gnmOrder=[('L',ld[0],l_id,lhit['dir']),('R',rd[0],r_id,rhit['dir'])]
					pair = [lhit, rhit]
					p0, p1 = pair
					if p0['qL'] > (p1['qL'] - uniqlen) or p0['qR'] > (p1['qR'] - uniqlen): # Reject non-shifted configuration
						continue
					if hitpairs and (hitpairs['hit1'][1]==0 and hitpairs['hit2'][1]==0 and (ld[0]!=0 or rd[0]!=0)): # if hitpair for read already goes to end, don't search for 2nd site
						continue
					overlap = p0['qR'] - p1['qL'] + 1 # length of sequence shared between pair
					sco = p0['sco'] + p1['sco']

					if overlap > 0:
						adjustments=[] #minimize overlap
						for pos in range(2):
							btop=pair[pos]['btop']
							pi=pair[pos]['id']
							if int(pi)<97:
								if not re.search(r'^[0-9]+[ACTGN\-]{0,2}[0-9]+$',btop):
									mismatches=re.findall(r'[ACTGN\-]{2}',btop)
									if mismatches:
										ambiguous=[i for i in range(len(mismatches)) if 'N' in mismatches[i]]
										if len(mismatches)-len(ambiguous)<=1:
											adjustments.append((0,0,btop))
											continue
								adjustments.append(minimize_overlap(pos,overlap,btop))
							else:
								adjustments.append((0,0,btop))
						(goff0,roff0,nbtop0),(goff1,roff1,nbtop1)=adjustments
						b4=[', '.join([f"{key}:{p[key]}" for key in ['qL','qR','sS','sE']]) for p in [p0,p1]]
						b4b=[p0['btop'],p1['btop']]
						p0['sE'] -= goff0 * p0['dir']
						p0['qR'] -= roff0
						p1['sS'] += goff1 * p1['dir'] #this was -=, but should shift forward if dir='+'
						p1['qL'] += roff1
						overlap -= (goff0 + goff1)
						p0['offset']+=roff0
						p1['offset']+=roff1
						p0['btop']=nbtop0
						p1['btop']=nbtop1
					if overlap < 0 or (hitpairs and overlap < hitpairs['overlap'] and int(hitpairs['sumPID'])>int(p0['id']+p1['id'])):
						continue
					for hit_idx in [0, 1]:
						dna = pair[hit_idx]['dna']
						for key in ['sE', 'sS', 'sL', 'sR']:
							before=f"{read}.{key}:{pair[hit_idx][key]}"
							if dna not in dnas or 'len' not in dnas[dna]:
								raise ValueError(f"DNA length missing for {dna}")
							pair[hit_idx][key] = register(pair[hit_idx][key], dnas[dna]['len'])
							after=f"{read}.{key}:{pair[hit_idx][key]}"
							if before!=after:
								print(f"before= {before}\nafter= {after}")
# shift:[0] 1 if pair[0] on left end of read, else =-1; [1] pair[0] jxn coord; [2] pair[1] jxn coord; [3] overlap; [4] shortest span over origin? [5] shortest span around circular replicon  
# config [0] hit jxn orient 0; [1] ... 1; [2] concat [0].[1]; [3] circle/scar/pal call;     # Config describes first the lower genomic-position hit, then the higher; + means L$

					#shift=[hit direction on genome,leftHalfGenomeCoord,RightHalfGenomeCoord,overlap,samecontig?,]
					shift = [1, pair[0]['sE'], pair[1]['sS'], overlap, 0, -1] 
				# boolean values convert to integers: False->0, True->1
					comp = ( int(pair[0]['dna'] > pair[1]['dna']) - int(pair[0]['dna'] < pair[1]['dna']),
							 1 if shift[1] > shift[2] else 0 if shift[1]==shift[2] else -1 )
					if comp[0] == 1 or (comp[0] == 0 and comp[1] == 1):
						# Reverse order
						gnmOrder.reverse()
						pair.reverse()
						p0, p1 = p1,p0 #p0, p1 = pair
						shift = [-1, shift[2], shift[1], overlap, 0, -1]
					if comp[0] == 0: # this means the two hitpairs are on the same contig
						shift[4] = 1
						shift[5] = shift[2] - shift[1]
						otherway = dnas[p0['dna']]['len'] - shift[5]
						if otherway < shift[5]:
							shift[4] = -1
							shift[5] = otherway
					config = ['+', '+', '', '']
					#orientation of hit * order in genome
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
								hitseq[offset] = "^"
							else:
								try:
									hitseq[offset] = hitseq[offset].lower()
								except Exception as e:
									print(f"Stopped due to {e}\noffset:{offset}\nhitseq length:{len(hitseq)}\nbtop: {hit['btop']},sequence so far: {''.join(hitseq)}")
									sys.exit()
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
						prox,dist=[func([p[key] for p in [p0,p1] for key in ['sS','sE']]) for func in [min,max]]
						if abs(dist-prox)<2000:
							continue
					if palflag and config[3] == 'palindrome':
						continue

					mobct = sum(1 for h in [0, 1] if tnend[h][-2] == '.')
					if not config[3] and mobct == 1:
						config[3] = 'transpose'
					if not config[3]:
						config[3] = '.'

					cts[config[3]] += 1
					second = config[2][1] #this part of the code is switching orient of second?
					second = '+' if second == '-' else '-'
					config[2] = config[2][0] + second
					#offset is 0-based index adjusted qL: offset=(qL-1)
					order = sorted([p0['offset'], p0['offset'] + p0['len'], p1['offset'], p1['offset'] + p1['len']])
					#workinghere
					jL, jR = order[1], order[2] #segment of shared sequence

					unsharedL=order[1]-order[0]
					unsharedR=order[3]-order[2]

					jxn = seq[jL - 5:jL] + '/' + seq[jL:jR] + '/' + seq[jR:jR + 5]

					seq = seq[:p1['qL'] - 1] + '-' + seq[p1['qL'] - 1:p1['qR']] + '-' + seq[p1['qR']:]

					if shift[0] == -1:
						seq = revcomp(seq)
						jxn = revcomp(jxn)
						align = [revcomp(a) for a in align]

					ct = 0
					if hitpairs and sco == hitpairs['sco']:
						ct = hitpairs['sites'] + 1

					read_termini = ["1" if i else "0" for i in [lhit['qL']==1,
							rhit['qR']==len(rhit['seq'])]]
					if p0!=lhit:
						read_termini.reverse()
					dna_termini=["0","0"]
					#for i,hit in [(0,lhit),(1,rhit)]:
					for i,hit in enumerate(pair):
						minc,maxc=[func([hit["sL"],hit["sR"]]) for func in [min,max]]
						if minc==1:
							dna_termini[i]="1"
						if abs(dnas[hit['dna']]['len']-maxc)<=1:
							dna_termini[i]="1"
					if not hitpairs or sco > hitpairs['sco']: # or length < hitpairs['len'] or int(hitpairs['sumPID'])<int(p0['id']+p1['id']):
						if hitpairs:
							losingsite.add('/'.join([str(hitpairs[lab]) for lab in ['dna1','coord','dna2','coordp','config']]))
						else:
							losingsite=set()
						hitpairs = {'sco': sco,'len':length,'dna1':p0['dna'],'coord':shift[1],'coordp':shift[2],
								'dna2':p1['dna'],'sumPID':p0['id']+p1['id'],'config':config[2],
								'readseq':seq,'mob1':tnend[0][-2],'mob1off':tnend[0][-1],'mob2':tnend[1][-2],
								'mob2off':tnend[1][-1],'overlap':shift[3],'hitseq1':align[0],'hitseq2':align[1],
								'jxn':jxn,'call':config[3],'origin':shift[4],'hit1':gnmOrder[0],
								'hit2':gnmOrder[1],'sites':ct,'losingsites':losingsite,'hitpair':pair,'cycle':cycle,'IntID':{},
								'read_termini':','.join(read_termini), 'dna_termini':','.join(dna_termini)}
					else:
						hitpairs['sites'] += 1
						hitpairs['losingsites'].add('/'.join([p0['dna'],str(shift[1]),p1['dna'],str(shift[2]),config[2]]))
	if hitpairs:
		cts['shift'] += weight
		cts['overlap'][hitpairs['overlap']]+= weight

		if hitpairs['len'] == xslen:
			hitpairs['len'] = -1
		(dna1,id1,sL1,sR1,dir1),(dna2,id2,sL2,sR2,dir2)=[[hitpairs['hitpair'][i][key] for key in ['dna','id','sL','sR','dir']] for i in range(2)]
		label = f"{dna1}/{hitpairs['coord']}/{dna2}/{hitpairs['coordp']}/{hitpairs['config']}"
		if label not in juxtas or juxtas[label]['sco'] < hitpairs['sco']:
			ct = 0 if intfilter else 1
			if label in juxtas:
				if not intfilter or (intfilter and cycle==3):
					ct += juxtas[label]['ct']
			else:
				juxtas[label]={}
			savekeys=['sco','coord', 'coordp', 'config', 'call','overlap', 'origin', 'mob1', 'mob1off',
				'mob2', 'mob2off', 'jxn', 'readseq', 'hitseq1','hitseq2', 'sites','read_termini', 
				'dna_termini','hit1','hit2','losingsites']
			for key in savekeys:
				juxtas[label][key]=hitpairs[key]
			juxtas[label].update({'read': read, 'ct': ct,'q1R': hitpairs['mob1off'],'q2L': hitpairs['mob2'],
				'hit1dna':dna1,'hit1id':id1,'hit1scoords':(sL1,sR1),'hit2dna':dna2,'hit2id':id2,
				'hit2scoords':(sL2,sR2)}) 
			if 'cycle' in juxtas[label]:
				if juxtas[label]['cycle']==cycle or (cycle==0 and juxtas[label]['cycle']==1):
					return
			if intfilter:
				IntID=makelink(juxtas[label],label)
				if IntID:
					juxtas[label]['IntID']=IntID
					if cycle==0:
						juxtas[label]['cycle']=1
					else:
						juxtas[label]['cycle']=cycle
				else:
					juxtas[label]['cycle']=0
		else:
			if not intfilter or (intfilter and cycle==3):
				juxtas[label]['ct'] += 1
			if intfilter:
				if juxtas[label]['cycle']==cycle or (cycle==0 and juxtas[label]['cycle']==1):
					return
				juxtas[label]['cycle']=cycle
				IntID=makelink(juxtas[label],label) #if new cycle check for new connections
				if IntID and not 'IntID' in juxtas[label]:
					juxtas[label]['IntID']=IntID
				elif IntID:
					for newmob in IntID:
						if not newmob in juxtas[label]['IntID']:
							juxtas[label]['IntID'][newmob]=IntID[newmob]
def connect_juxtas(oldjuxtas,cycle):
	newjuxtas={}
	for lab in oldjuxtas:
		newjuxtas[lab]=oldjuxtas[lab]
		IntID=makelink(newjuxtas[lab],lab)
		if IntID and not 'IntID' in newjuxtas[lab]:
			newjuxtas[lab]['IntID']=IntID
			newjuxtas[lab]['cycle']=cycle
		elif IntID:
			for newmob in IntID:
				if not newmob in newjuxtas[lab]['IntID']:
					newjuxtas[lab]['IntID'][newmob]=IntID[newmob]
					newjuxtas[lab]['cycle']=cycle if not newjuxtas[lab]['cycle'] else newjuxtas[lab]['cycle']
	return newjuxtas

def makelink(juxline,juxlink): #,mobSeg,tempmobSeg,oriswitch):
	juxlinkdict=juxtalabel2dict(juxlink)
	intinfo=defaultdict(list)
	hits=[[juxline[key] for key in f"hit{n}dna hit{n}scoords hit{n}".split()] for n in range(1,3)]
	for i,(dna,(left,right),guide) in enumerate(hits):
		readside,dir=guide[0],guide[-1]
		intgenes=findnearbyInt(dna,left,right) #,readside,dir)
		if intgenes:
			linkerdna,(linkerL,linkerR),linkerGuide=hits[0] if i==1 else hits[1]
			for intdna,pos_i,intToTheRight,gapsize in sorted(intgenes,key=lambda x:x[-1],reverse=True): # sort mobgenes by closest proximity
				mobdna,mob_id,L,R,mobdir,midpoint,history=[mobSeg[intdna][pos_i][key] for key in keyorder]
				mobside='B' if mobdir=="+" else 'E'
				#if 'IntID' in juxline:
				#	if mob_id in juxline['IntID']:
				#		continue # skip hit if the mobility gene already has a hit
				#if primary then mob code=> hit side, dir on genome, intside
				#if secondary then mob code=> linkage side oriented to hit, dir on genome,intside
				if history: #secondary
					primary=False
					last_mobcode=mobSeg[intdna][pos_i].get("last_mobcode")
					last_side_of_uniq,mobside=last_mobcode[0::2]
					lastlinkcoord=mobSeg[intdna][pos_i].get("linking_coord") #the link from primary hit
					lastoverlap=mobSeg[intdna][pos_i].get("last_overlap") #overlap in primary connection
					hitcoord=juxlinkdict.get(dna) #hit in jux that is near the lastlinkcoord
					hitdir=guide[-1]
					side_of_uniq="R" if hitcoord==left else "L" #genomic side of uniq seq
					xlen=abs(hitcoord - lastlinkcoord )+1 # putative cross over length
					if xlen  >= juxline['overlap'] or xlen >= lastoverlap: 
					# if lastlinkcoord at or beyond overlap distance to hitcoord then the side_of_uniq switches
						mobcode=f"{oriswitch.get(last_side_of_uniq)}" 
					else:
						mobcode=last_side_of_uniq
					if mobcode!=side_of_uniq: # based on co-orientation of overlapping sequence
						hitdir=oriswitch.get(hitdir) 
					mobcode += "0" if hitdir==1 else "1"
				else: #primary
					primary=True
					intside=["L","R"][int(intToTheRight)]
					hitside,hitdir=guide[0],guide[-1]
					linkside,linkdir=linkerGuide[0],linkerGuide[-1]
					side_of_uniq ="R" if juxlinkdict.get(dna)==left else "L"
					if intside!='R': # coorient so that integrase is on the right
						mobdir,hitdir,hitside,linkdir,linkside,side_of_uniq,mobside=[oriswitch.get(val) for val in [mobdir,hitdir,hitside,linkdir,linkside,side_of_uniq,mobside]]
					mobcode=f"{side_of_uniq}{0 if hitdir==1 else 1}"
					linking_coord=juxlinkdict.get(linkerdna)

				try:
					mobcode+=mobside
					new_history=[{juxlink:mobcode}]+history
					intinfo[mob_id].append(new_history)
				except Exception as e:
					print("huh",cycle,e,type(history),history,juxline['read'])
					sys.exit()
				if primary or runcyc>2:
					remainsLeft,remainsRight=[gapsize+juxline['overlap']]*2
					tempmobSeg[linkerdna].append({'dna':linkerdna,'id':mob_id,'L': linkerL, 
						'R': linkerR,'dir':mobdir,'midpoint':int((linkerL+linkerR)/2),
						'history':new_history,'last_mobcode':mobcode,
						'linking_coord':linking_coord,"last_overlap":juxline['overlap']})
	#there can be more than one link to the same mobile gene, choose one with a shared sequence if possible:
		#given a sample assembly with  species 1 and 2: attB1, attL2, attR2
		#The following connections of attR will happen in cycle 2, if attP recombinant read in readset:
			#1. attPL-attBR could be linked by primary attPL :attPR (attPL shared sequence) 
			#2. attPL-attBR could be linked by primary attBR2:attPR (no shared sequence, indirect connection via attB)
		#1 will be chosen over #2
		#It attP was not present, #1 would not be possible, so #2 would be kept

	info={}
	for mob in intinfo:
		if len(intinfo[mob])==1:
			info[mob]=intinfo[mob][0]
		elif len(intinfo[mob]) >1:
			for mobconnect in intinfo[mob]:
				side_of_uniq=set([link[key][0] for link in mobconnect for key in link])
				if len(side_of_uniq)==2:
					info[mob]=mobconnect
					break
			if not mob in info:
				info[mob]=intinfo[mob][0]
	return info

def write_juxtas_bed():
	bed_file = os.path.join(outDir, "juxtas.bed")
	with open(bed_file, 'w') as out:
		overlaps = {}
		for label, line in juxtas.items():
			if not re.search(r'circle|scar',line['call'],re.I):
				continue
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

def minimize_overlap(pos,overlap,btop):
	goff=roff=last=0
	newbtop=btop
	while btop:
		m = re.search(mo_patterns[pos], btop)
		if m:
			last = int(m.group(1))
			btop = btop[: -len(m.group(1))] if pos==0 else btop[len(m.group(1)):]
			if last + roff > overlap:
				break
		if len(btop) < 2:
			break
		pair_chars = btop[-2:] if pos==0 else btop[:2]
		btop = btop[:-2] if pos==0 else btop[2:]
		roff += last + 1
		goff += last + 1
		last = 0
		if pair_chars[0] == '-':
			roff -= 1
		elif pair_chars[1] == '-':
			goff -= 1
		newbtop=btop
		if roff > overlap:
			break
	return goff,roff,newbtop
def collapsejuxtas():
	global juxtas
	hotspots=defaultdict(list)
	for label in juxtas:
		dna1,c1,dna2,c2,config=label.split('/')
		b1,b2=[genomebins.get(c) for c in [c1,c2]]
		hotspots[f"{dna1}/{b1}/{dna2}/{b2}/{config}"].append((label,float(juxtas[label]['hit1id'])+float(juxtas[label]['hit2id'])))
	for hotspot in hotspots:
		byid=sorted(hotspots[hotspot],key=lambda x: x[1],reverse=True)
		bestlab,bestid=byid[0]
		if len(byid)>1:
			dna1,bc1,dna2,bc2,config=bestlab.split('/')
			bc1,bc2=int(bc1),int(bc2)
			delete=[]
			for rank,(label,id) in enumerate(byid):
				if label!=bestlab:
					dna1,c1,dna2,c2,config=label.split('/')
					c1,c2=int(c1),int(c2)
					if all([i<min([5,int(juxtas[bestlab]['overlap'])]) for i in [abs(bc1-c1),abs(bc2-c2)]]):
						delete.append((label,id))
			for label in delete:
				if 'IntID' in juxtas[label[0]]:
					intid=juxtas[label[0]]['IntID']
					replaceme={}
					for mob in intid:
						for link in intid[mob]:
							if label[0] in link:
								replaceme[mob]={bestlab:link[label[0]]}
					for line in juxtas:
						if line==label[0]:
							continue
						for key in replaceme:
							if key in juxtas[line].get("IntID",[]):
								for reind,lab in enumerate(juxtas[line]["IntID"][key]):
									if not label[0] in lab:
										continue
									juxtas[line]["IntID"][key][reind]=replaceme[mob]
				juxtas[bestlab]['ct']+=juxtas[label[0]]['ct']
				del juxtas[label[0]]
				try:
					ind= byid.index(label)
					del byid[ind]
				except:
					print("wuh",delete)
					print(label)
					print(byid)
					print(bestlab,bestid)
		for rank,(label,id) in enumerate(byid):
			juxtas[label]["best_local"]=bestlab
			juxtas[label]["local_count"]=len(byid)
			juxtas[label]["local_rank"]=rank+1


def print_counts():
	for key in ['nonstandard', 'genomic', 'multi', 'double', 'shift', 'circleJxn', 'scar', 'palindrome', 'transpose', 'tn']:
		print(f"{key}={cts[key]}", end=' ')
	print()
	if 'overlap' in cts:
		for k in sorted(list(cts['overlap'].keys())):
			print(f"overlap {k}={cts['overlap'][k]}")
# Main execution
if __name__ == "__main__":
	keyorder=['dna','id', 'L', 'R', 'dir','midpoint', 'history']
	load_ref_files()
	tempmobSeg=defaultdict(list)
	per_read()
	if intfilter:
		done=False
		for cycle in range(2,runcyc+1):
			mobSeg=tempmobSeg
			tempmobSeg=defaultdict(list)
			if cycle==runcyc:
				for key in cts:
					if not key in ['nonstandard','mobs','tn']:
						cts[key]=0
				cts['overlap']=defaultdict(int)
			per_read()
		collapsejuxtas()
	else:
		collapsejuxtas()
		juxtas=connect_juxtas(juxtas,cycle)
		cycle+=1
		mobSeg=tempmobSeg
		tempmobSeg=defaultdict(list)
		juxtas=connect_juxtas(juxtas,cycle)
	with open(juxtasout,'wb') as p:
		pickle.dump(juxtas,p)
	print_counts()

