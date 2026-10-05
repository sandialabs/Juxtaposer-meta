#!/usr/bin/env python3

import sys,os,re,subprocess,pickle,math,json
from collections import defaultdict, Counter
import json
# Usage check
if len(sys.argv) < 3:
	print(f"Usage: python {sys.argv[0]} juxpickle readmax")
	sys.exit(1)

jux_file,readmax = sys.argv[1:]
readmax=int(readmax)

if not os.path.isfile(jux_file):
	sys.exit(f"Could not open {jux_file}")

base,fn=os.path.split(jux_file)
basedir=os.listdir(base)
for i in basedir:
	if i.endswith(".mobile.bed"):
		prefix=i.replace(".mobile.bed","")
		break
reffiles=[os.path.join(base,filename) for filename in f"{prefix}.lens {prefix}.nin {prefix}.mobile.bed".split()]
for file in reffiles:
	if not os.path.exists(file):
		sys.exit() #add print message
lens_file,blastdbs,mob_file=reffiles
blastdbs=blastdbs.rstrip(".nin")
tempfa,tempout=[os.path.join(base,f"clusteroverlaps.{suf}") for suf in "fa blast".split()]
if os.path.exists(tempfa):
	os.remove(tempfa)

juxtasout=jux_file.replace('.pkl','.txt')

networkout=os.path.join(base,"networks.json")
cutoff=70
class jux:
	def __init__(self, **kwargs):
		for key, val in kwargs.items():
			setattr(self, key, val)
	@classmethod
	def from_dict(cls, data_dict):
		return cls(**data_dict)
	def __repr__(self):
		keys = sorted(self.__dict__.keys())
		attrs = ', '.join(f"{k}={repr(getattr(self, k))}" for k in keys)
		return f"{self.__class__.__name__}({attrs})"


def readfile(file,obj_type=open):
	with obj_type(file,'r') as f:
		for line in f:
			yield line.strip()
def revcomp(seq):
	complement = str.maketrans('ACGTacgt', 'TGCAtgca')
	return seq.translate(complement)[::-1]

def comseq(seq1,seq2):
	rseq1,rseq2=revcomp(seq1),revcomp(seq2)
	if len(seq1)!=len(seq2):
		shorter,longer=int(len(seq1)>len(seq2)),int(len(seq1)<len(seq2))
		shorterseq,longerseq=[[seq1,seq2][i] for i in [shorter,longer]]
		rshorter=rseq1 if not shorter else rseq2
		if shorterseq in longerseq:
			return (0,0), shorterseq,[(shorter,0),(longer,longerseq.find(shorterseq))]
		elif rshorter in longerseq:
			return (1,0), rshorter,[(shorter,0),(longer,longerseq.find(rshorter))]
	kmers={(i,0):0 for i in [0,1]}
	kmer=min([max([2,int(min([len(seq1),len(seq2)])/3)]),5])
	for i,j in kmers:
		s=[seq1 if not i else rseq1, seq2 if not j else rseq2]
		store=defaultdict(set)
		for si in range(2):
			for pos in range(0,len(s[si])-kmer):
				store[ s[si][pos:pos+kmer] ].add(si)
		kmers[(i,j)]=sum([1 for subseq in store if len(store[subseq])==2])
	seqOris=sorted(list(kmers.items()),key=lambda x:x[1])[-1][0]
	seq1,seq2=[seq1 if not seqOris[0] else rseq1, seq2 if not seqOris[1] else rseq2]
	if seq1 in seq2:
		return seqOris, seq1,[(0,0),(1,seq2.find(seq1))]
	if seq2 in seq1:
		return seqOris, seq2,[(1,0),(0,seq1.find(seq2))]
	subseqs=defaultdict(list)
	mink=int(min([len(seq1),len(seq2),6])*0.5)
	for k in range(min([len(seq1),len(seq2)]),mink-1,-1):
		for index,seq in enumerate([seq1,seq2]):
			for i in range(0,len(seq)-k+1):
				subseq=seq[i:i+k]
				subseqs[subseq].append((index,i))
		for sub,records in subseqs.items():
			if len(records)>1:
				if len(set([record[0] for record in records]))>1:
					return seqOris,sub,records
		subseqs=defaultdict(list)
	return None,None,None

def load_ref_files():
	global phgint,tiles
	biggest=0
	for line in readfile(lens_file):
		line = line.rstrip('\n')
		m = re.match(r'^(\S+).*\t(\d+)\s*', line)
		if m:
			key, size= m.groups()
			dnas[key] = int(size)
			if biggest<int(size):
				biggest=int(size)
	for line in readfile(mob_file):
		contig,L,R,id_,score,dir=line.split('\t')
		mobsLoc[id_]={contig:{'L':int(L),'R':int(R)}}
def iterthru(temp,target):
	for key,val in temp.items():
		if isinstance(val,str):
			if key in target: #and target[key]!='':
				if not val or target[key]==val: #nothing to add
					continue
				else:
					target[key].update({val:""}) # typo ,->: and changed .add method to .update
			else:
				target[key]={val:""} if val else val
		if isinstance(val,dict):
			if key in target and not isinstance(target[key],str):
				target[key]=iterthru(val,target[key])
			else:
				target[key]=val
	return target
def getstats(data):
	outliers=[]
	if not len(data):
		return 
	sorted_data=sorted(data)
	mindata,maxdata,sumdata,lendata=[f(data) for f in [min,max,sum,len]]
	mu=sumdata/lendata
	sumsqdifs=sum([(mu-float(ex))**2 for ex in data])
	variance=sumsqdifs/lendata
	stdev=variance**0.5
	rangedata=maxdata-mindata
	if lendata<3:
		median=data
		if mindata*2<maxdata:
			outliers=data
	else:
		if lendata % 2 == 0:
			median=[sorted_data[int((lendata+i)/2)] for i in [-1,1]]
		else:
			median=[sorted_data[int((lendata/2)+0.5)]]

		l,u=[int(mu+(i*stdev)) for i in [-1,1]]
		if rangedata>10:
			for point in data:
				if not (point>=l and point<=u):
					outliers.append(point)

	return {"len":lendata,"sum":sumdata,"average":mu,"variance":variance,"stdev":variance**0.5,"range":rangedata,'median':median,'data':data,'outliers':outliers}


def getjuxdict(label):
	parts=label.split('/')
	if len(parts)==5 or (len(parts)==4 and parts[-1].isnumeric):
		return {parts[i]:parts[i+1] for i in [0,2]}
	return dict()
def read_juxtas():
	pruned=0
	tot=0
	for label in juxtas:
		juxobj=jux.from_dict(juxtas[label])
		replacelabs[label]=juxobj.best_local
		rterm,dterm=juxobj.read_termini,juxobj.dna_termini
		dnatermd[label]=dterm
		seqs[label]=(juxobj.jxn,juxobj.config)
		overlaps[label]=juxobj.overlap
		coverages[label]+=juxobj.ct
		sumpid[label]=juxobj.hit1id+juxobj.hit2id
		multiplicity[label]=juxobj.local_count
		losingsites[label]=juxobj.losingsites
		if not juxobj.overlap:
			questionable[label]="overlap"
			continue
		if not "IntID" in juxtas[label]:
			for dna,posx in [(juxobj.hit1dna,juxobj.coord),(juxobj.hit2dna,juxobj.coordp)]:
				if not posx in negatives[dna]:
					negatives[dna][posx]=[]
				negatives[dna][posx].append(label)
			continue
		if label!=juxobj.best_local:
			coverages[juxobj.best_local]+=juxobj.ct
		if rterm!="1,1":
			rterms,dterms=[t.split(',') for t in [rterm,dterm]]
			questionable[label]="read-internal"
			continue
		if len(juxobj.losingsites)>5:
			questionable[label]=f"{len(juxobj.losingsites)} losers"
			continue
		if juxobj.call!=".":
			calls[label]=juxobj.call
		mobs=juxobj.IntID
		for mob in mobs:
			if not mob in director:
				director[mob]=defaultdict(list)
			links={replacelabs.get(j,j):i[j] for i in mobs[mob] for j in i}
			baselinks=list(links.keys())
			mylinks=defaultdict(list)
			temp,trim=None,None
			if len(links)==2:
				for i,n in enumerate(baselinks):
					c1,dna1,c2,dna2,cf=n.split('/')
					mylinks[c1].append((i,dna1))
					mylinks[c2].append((i,dna2))
				for n in mylinks:
					if len(mylinks[n])==2:
						dif=abs(int(mylinks[n][0][-1])-int(mylinks[n][1][-1]))
						if dif>readmax:
							trim=baselinks[0]
			for link in links:
				if trim:
					if link==trim:
						break
				director[mob][link].append(links[link])
				if temp:
					temp={link:temp}
				else:
					temp={link:""}
			if not mob in story and temp:
				story[mob]=temp
			elif not temp:
				continue
			else:
				iterthru(temp,story[mob])
	tempstory={}
	for mob,links in story.items():
		nodecount=sum([len(links[p])+1 for p in links])
		if nodecount<2:
			continue
		tempstory[mob]=defaultdict(dict)
		pp,ps=0,0
		try:
			for p in story[mob]:
				bestp=replacelabs.get(p,p)
				if bestp in questionable:
					bestp=p
				if not bestp in tempstory[mob]:
					if story[mob][p]:
						tempstory[mob][bestp]={i:{} for i in story[mob][p]}
					else:
						tempstory[mob][bestp]={}
				if not len(story[mob][p]):
					continue
				for s in story[mob][p]:
					bests=replacelabs.get(s,s)
					if bests in questionable:
						bests=s
					tempstory[mob][bestp][bests]={}
		except Exception as e:
			sys.stderr.write(f"{e}, line 255\n")
	clusternets(tempstory)
	return

def clusternets(story):
	for mob in story:
		story_id={side:defaultdict(dict) for side in "BE"}
		tertiaries={}
		seqdict={side:defaultdict(list) for side in "BE"}
		mobcodes=director.get(mob)
		for p in story[mob]:
			pmc=mobcodes.get(p)[0][-1]
			pcross=seqs.get(p)[0].split('/')[1].lower()
			if len(pcross)<=3:
				if revcomp(pcross) in story_id[pmc]:
					pcross=revcomp(pcross)
				elif not pcross in story_id[pmc]:
					story_id[pmc][pcross]={}
				story_id[pmc][pcross][p]={}
			else:
				pfullseq=seqs.get(p)[0].replace('/','').lower()
				seqdict[pmc][pfullseq].append(p)
			for s in story[mob][p]:
				scross=seqs.get(s)[0].split('/')[1].lower()
				if len(pcross)<=3 and len(scross)<=3:
					story_id[pmc][pcross][p][s]={}
				else:
					sfullseq=seqs.get(s)[0].replace('/','').lower()
					seqdict[pmc][sfullseq].append(s)
				juxdict_s=[(key,val) for key,val in getjuxdict(s).items() if not key in getjuxdict(p)]
				if juxdict_s:
					cont,posx=juxdict_s[0]
					if cont in negatives:
						for psx,lab in negatives[cont].items(): # these are juxtapositions not linked to mobile gene after two cycles
							if abs(int(posx)-int(psx))<overlaps.get(s)+5: #has to be within id-block+ wiggle room
								tertiaries[s]=sorted(lab,key=lambda x: sumpid[x],reverse=True)[0]
		for pmc in seqdict:
			if not len(seqdict[pmc]):
				continue
			connectionswithseq=seqdict[pmc]
			for i in range(2):
				log=set()
				shared=defaultdict(set)
				sorted_seqs=sorted(list(connectionswithseq.items()),key=lambda x: len(x[0]),reverse=True)
				for seq,connections in sorted_seqs:
					for seq2,connections2 in sorted_seqs:
						if seq==seq2:
							continue
						rcstatus,shared_seq,positions=comseq(seq,seq2)
						if not shared_seq or len(shared_seq)<int(min([len(seq),len(seq2)])/2):
							continue
						if revcomp(shared_seq) in shared:
							shared_seq=revcomp(shared_seq)
						shared[shared_seq].update(set(connections+connections2))
						log.update(shared[shared_seq])
				addlog={connection:key for key,connections in sorted_seqs for connection in connections if not connection in log}
				if shared:
					replacethese=defaultdict(set)
					replaceiter={}
					subnets=set()
					for seq,con in sorted(list(shared.items()),key=lambda x: len(set(x[1])),reverse=True):
						if all([c in addlog for c in con]):
							continue
						for c in con:
							if c in addlog:
								if len(seq)<len(addlog[c]): 
									replacethese[seq].update(con)
									break
								else:
									replacethese[addlog[c]].update(con)
							addlog[c]=seq

					for seq in replacethese:
						for c in replacethese[seq]:
							addlog[c]=seq

					subnets=set([seq for c,seq in addlog.items()])
					if len(subnets)>1:
						for seq1 in subnets:
							for seq2 in subnets:
								if seq1==seq2 or len(seq1)==len(seq2):
									continue
								bothseqs=[seq1,seq2]
								long,short=bothseqs[int(len(seq1)<len(seq2))],bothseqs[int(len(seq1)>len(seq2))]
								if short in long or revcomp(short) in long:
									if not long in replaceiter or len(short)<len(replaceiter[long]):
										replaceiter[long]=short

					connectionswithseq=defaultdict(list)
					for c in addlog:
						if addlog[c] in replaceiter:
							addlog[c]=replaceiter[addlog[c]]
						connectionswithseq[addlog[c]].append(c)
			for p in story[mob]:
				if not p in addlog:
					continue
				seq=addlog[p]
				if not seq in story_id[pmc]:
					story_id[pmc][seq]={}
				story_id[pmc][seq][p]={}
				for s in story[mob][p]:
					if not s in addlog or seq!=addlog[s]:
						continue
					story_id[pmc][seq][p][s]={}
		analyzeclusters(story_id,mob,tertiaries)
	return

def analyzeclusters(story_id,mob,tertiaries):
	if story_id:
		overs=[]
		oversid=defaultdict(list)
	# network cluster exploration
		clusters={}
		dirs=director[mob]
		count_cross=0
		searchtype=1
		for side,crossover in story_id.items():
			cnodes=0
			for seq,p in crossover.items():
				cnodes+=len(p)
				for pri,s in p.items():
					cnodes+=len(s)
			if cnodes<2 or cnodes>18: #maxcount for hgt is 6 unique, assume multiplicity of 3
				continue
			for seq,pris in crossover.items():
				overlapdata=[]
				notes=[]
				uniqcontigs=Counter()
				log,contcoord=defaultdict(list),defaultdict(list)
				for p,secs in pris.items():
					overlaplen=overlaps.get(p,0)
					overlapdata.append(overlaplen)
					oversid[overlaplen].append(p)
					c1,p1,c2,p2=[p.split('/')[i] for i in range(4)]
					for cont,pos in [(c1,p1),(c2,p2)]:
						uniqcontigs[cont]+=1
						log[cont].append(int(pos))
						contcoord[f"{cont}/{pos}"].append(p)
					for s in secs:
						overlaplen=overlaps.get(s,0)
						overlapdata.append(overlaplen)
						oversid[overlaplen].append(s)
						c1,p1,c2,p2=[s.split('/')[i] for i in range(4)]
						for cont,pos in [(c1,p1),(c2,p2)]:
							uniqcontigs[cont]+=1
							log[cont].append(int(pos))
							contcoord[f"{cont}/{pos}"].append(s)
	#network analysis
				if len(uniqcontigs)>6:
					continue
				ostat=getstats(overlapdata)
				ignoreoverlaps={node:i for i in ostat['outliers'] for node in oversid[i]}
				if len(overlapdata)-len(ignoreoverlaps)<=1:
					continue
				splits=defaultdict(list)
				for contig in log:
					if len(log[contig])>1:
						splitcounts=Counter(log[contig])
						for split,count in splitcounts.items():
							splits[contig].append((split,count))
				nearterminal=defaultdict(list)
				genomicatt=[]
				complete=[]
				allbroken=True
				pred_gsites=defaultdict(list)
				pred_recsites={}
				for n in splits:
					for i in splits[n]:
						for j in splits[n]:
							if i==j:
								continue
							if (i,j) in complete or (j,i) in complete:
								continue
							complete.append((i,j))
							gensite=abs(j[0]-i[0])
							inrange=False
							if gensite in overlapdata:
								inrange=True
							else:
								buffer=gensite-1,gensite+1
								if (buffer[0]>=min(overlapdata) and buffer[0]<=max(overlapdata)) or (buffer[1]>=min(overlapdata) and buffer[1]<=max(overlapdata)):
									inrange=True
							if inrange:
								allbroken=False
								sitetype=None
								gen1,gen2=contcoord[f"{n}/{j[0]}"],contcoord[f"{n}/{i[0]}"]
								for p in pris:
									if p in gen1 or p in gen2:
										if sitetype==2 or sitetype==None:
											sitetype=1
											genomicside=gen1 if p in gen1 else gen2
									for s in pris[p]:
										if s in gen1 or gen2:
											if sitetype==None:
												sitetype=2
												genomicside=gen1 if s in gen1 else gen2
								if sitetype==1 or sitetype==2:
									pred_gsites[n].append((min([j[0],i[0]]),max([j[0],i[0]])))
				if pred_gsites:
					losers=defaultdict(set)
					for contig,sites in pred_gsites.items():
						if len(sites)<2:
							continue
						for site in sites:
							mp=abs(int((site[1]-site[0])/2))
							for site2 in sites:
								if site==site2:
									continue
								mp2=abs(int((site[1]-site[0])/2))
								if abs(mp-mp2)<readmax:
									difs=[abs(i[0]-i[1]) for i in [site,site2]]
									if difs[0]<difs[1] and (difs[1]-difs[0])>min(overlapdata):
										losers[contig].add(site2)
										for p in site2:
											if not p in nearterminal[contig]:
												nearterminal[contig].append(p)

										if difs[1] in oversid:
											for juxid in oversid[difs[1]]:
												if "1,1" in dnatermd.get(juxid):
													genomicatt.append(juxid)
													ignoreoverlaps[juxid]=difs[1]

					for contig,loser in losers.items():
						pred_gsites[contig]=[i for i in pred_gsites[contig] if i!=loser]
						if not pred_gsites[contig]:
							del pred_gsites[contig]
				elif allbroken:
					islandconts=[]
					for cont in uniqcontigs:
						if uniqcontigs[cont]==3:
							islandconts.append(cont)
					sites={}
					orgs=[]
					if not islandconts:
						notes.append(f"attP site could not be validated. Tertiary cycle required... {len(tertiaries)} tertiary sites found")
						searchtype=0
						uniqconts={k:uniqcontigs[k] for k in uniqcontigs}
						attp=Counter()
						cyclelog={1:[],3:[]}
						for key,val in tertiaries.items():
							tertjux=[(val.split('/')[i],val.split('/')[i+1]) for i in [0,2]]
							for t,ctrt in tertjux:
								attp[f"{t}/{ctrt}"]+=1
								if not t in uniqconts:
									uniqconts[t]=0
								uniqconts[t]+=1
								cyclelog[3].append(f"{t}/{ctrt}")
						for p in pris:
							c1,p1,c2,p2,C=p.split('/')
							attp[f"{c1}/{p1}"]+=1
							attp[f"{c2}/{p2}"]+=1
							cyclelog[1]+=[f"{c1}/{p1}",f"{c2}/{p2}"]
						remove=[cont for cont in attp if attp[cont]<2]
						for cont in remove:
							del attp[cont]
						islandconts=[i.split('/')[0] for i in attp]
						if len(islandconts)!=2:
							continue
							notes.append("island contigs could not be resolved, attP prediction is incomplete")
					for p in pris:
						test=[cont for cont in islandconts if cont in p]
						if len(test)==2:
							pred_recsites[p]="attP"
						elif len(test)==1:
							pconts=[p.split('/')[i] for i in [0,2]]
							orgcont=pconts[0] if pconts[1] in test else pconts[1]
							if not [orgcont] in orgs:
								orgs.append([orgcont])
							pred_recsites[p]=f"attL{orgs.index([orgcont])+1}"
					for p in pris:
						if not p in pred_recsites:
							continue

						for s in pris[p]:
							if pred_recsites[p]=="attP":
								stype="R"
								pred_recsites[s]="attR"
								c1,c2=[s.split('/')[i] for i in [0,2]]
							elif pred_recsites[p].startswith("attL"):
								stype="B"
								c1,c2=[s.split('/')[i] for i in [0,2]]
								test1,test2=[[i for i in orgs if c in i] for c in [c1,c2]]
								if test1:
									ind=orgs.index(test1[0])
									pred_recsites[s]=f"attB{ind+1}"
									orgs[ind].append(c2)
								elif test2:
									ind=orgs.index(test2[0])
									pred_recsites[s]=f"attB{ind+1}"
									orgs[ind].append(c1)
							if not searchtype and s in tertiaries:
								t=tertiaries.get(s)
								c1,c2=[t.split('/')[i] for i in [0,2]]
								test1,test2=[[i for i in orgs if c in i] for c in [c1,c2]]
								if test1:
									ind=orgs.index(test1[0])
									orgs[ind].append(c2)
									pred_recsites[t]=f"attR{ind+1}"
								elif test2:
									ind=orgs.index(test2[0])
									orgs[ind].append(c1)
									pred_recsites[t]=f"attR{ind+1}"
					if not searchtype:
						pred_recsites['/'.join([key for key in attp])]="attP"

					addlink,remove={},[]
					for site in pred_recsites:
						if pred_recsites[site]=="attR":
							c1,p1,c2,p2=[site.split('/')[i] for i in range(4)]
							orgcont,orgpos=[c2,p2] if c1 in islandconts else [c1,p1] if c2 in islandconts else [None,None]
							if orgcont==None:
								 continue #del pred_recsites[site]
							try:
								orgid=[orgs.index(i)+1 for i in orgs if orgcont in i]
								if not orgid:
									if not site in tertiaries:
										remove.append(site)
										continue
									addlink[tertiaries[site]]="attR"
								pred_recsites[site]=f"attR{orgid[0]}"
							except Exception as e:
								sys.stderr.write(f"{e}, line 591\n")
					for site in remove:
						del pred_recsites[site]
					for site in addlink:
						pred_recsites[site]=addlink[site]
					if len(pred_recsites)>7 or len(pred_recsites)<2:
						continue
				else:
					continue

				for p in pris:
					o,d=overlaps.get(p),dnatermd.get(p)
					if "1" in d:
						if o-min(overlapdata)>int(sum(overlapdata)/len(overlapdata)):
							ignoreoverlaps[p]=o
							for s in pris[p]:
								o,d=overlaps.get(s),dnatermd.get(s)
								if o-min(overlapdata)>int(sum(overlapdata)/len(overlapdata)):
									ignoreoverlaps[s]=o

				totalfirst=totalseconds=scar_cj=0
				localcountsum,overlapdata,covsum=[],[],[]
				uniqcontigs=Counter()
				finalnet={}
				totalthirds=0
				directed_net={}
				for p,secs in pris.items():
					if not p in ignoreoverlaps:
						overlaplen=overlaps.get(p,0)
						overlapdata.append(overlaplen)

					totalfirst+=1
					localcountsum.append(multiplicity.get(p,0))
					covsum.append(coverages.get(p,0))
					if p in calls:
						scar_cj+=1
					c1,c2=[p.split('/')[i] for i in [0,2]]
					uniqcontigs[c1]+=1
					uniqcontigs[c2]+=1
					finalnet[p]={}
					directed_net[p]=list(set(dirs.get(p)))
					for s in secs:
						if not s in ignoreoverlaps:
							overlaplen=overlaps.get(s,0)
							overlapdata.append(overlaplen)
						totalseconds+=1
						localcountsum.append(multiplicity.get(s,0))
						covsum.append(coverages.get(s,0))
						if s in calls:
							scar_cj+=1
						c1,c2=[s.split('/')[i] for i in [0,2]]
						uniqcontigs[c1]+=1
						uniqcontigs[c2]+=1
						finalnet[p][s]=str()
						directed_net[s]=list(set(dirs.get(s)))
						if not searchtype and s in tertiaries:
							t=tertiaries[s]
							overlaplen=overlaps.get(t,0)
							overlapdata.append(overlaplen)
							totalthirds+=1
							localcountsum.append(multiplicity.get(t,0))
							covsum.append(coverages.get(t,0))
							if t in calls:
								scar_cj+=1
							c1,c2=[t.split('/')[i] for i in [0,2]]
							uniqcontigs[c1]+=1
							uniqcontigs[c2]+=1
							finalnet[p][s]=t

				X1,X2,X3=totalfirst,totalseconds,len(nearterminal)+len(pred_gsites)
				leaf=None
				if pred_gsites:
					tests=[X2<=4, X1<=3, X3<=4, X3<=9]
					if tests[0] or not tests[1] or (tests[1] and not tests[2] and tests[3]):
						pos=True
					else:
						pos=False
				elif allbroken:
					tests=[X1>2,X1<=1,X2<=1,X2<=4]
					if tests[0] or (not tests[1] and not tests[2] and tests[3]):
						pos=True
					else:
						pos=False
				overlapstat=getstats(overlapdata)
				if not overlapstat:
					overlapstat=ostat
				if not overlapstat:
					continue
				if ignoreoverlaps:
					overlapstat['outliers']=ignoreoverlaps
				else:
					overlapstat['outliers']={node:i for i in overlapstat['outliers'] for node in oversid[i]}
				covstat=getstats(covsum)
				localstat=getstats(localcountsum)
				crossandnt=defaultdict(list)
				for contig in nearterminal:
					if contig in pred_gsites:
						ntmp=sum(nearterminal[contig])/2
						for coords in pred_gsites[contig]:
							xmp=sum(coords)/2
							if abs(xmp-ntmp)<readmax:
								crossandnt[contig].append(coords)
				accounted_overlap_outliers=[jux for jux in overlapstat['outliers'] if any([halflab in jux for halflab in [f"{key}/{coord}" for key in nearterminal for coord in nearterminal[key]]])]
				if pred_gsites:
					modeltype=1
					unaccounted_outliers=len(overlapstat['outliers'])-len(accounted_overlap_outliers)
#					c1,c2,c3,b=(3.6259904732490797,-2.1790068918313783,-0.7276079610398841,0.6646529138187754)
					c1,c2,c3,b=(2.612901701004964,-0.9368229945890589,-0.7538679871906498,-2.2320736221579254)
					x1,x2,x3=[int(i) for i in [totalfirst,totalseconds,unaccounted_outliers]]
				elif pred_recsites:
					modeltype=2
#					c1,c2,c3,b=(0.35526901602611805,0.4014174857330322,1.1119555177852682,-5.4899196950553995)
					c1,c2,c3,b=(0.017788679453884324,1.0441824352224296,1.2052296862337888,-7.5356735633569505)
					x1,x2,x3=[int(i) for i in [totalfirst,totalseconds,len(pred_recsites)]]
				else:
					continue
				count_cross+=1
				z=c1*x1+c2*x2+c3*x3+b
				percentage=int(round(1/(1+math.exp(-z))*100))
				serialid=f"{mob}_{side}{count_cross}{modeltype}-{percentage}"
				for p in finalnet:
					if not p in elfin:
						elfin[p].append(serialid)
					for s in finalnet[p]:
						if not s in elfin: 
							elfin[s].append(serialid)
						for t in finalnet[p][s]:
							if not t in elfin:
								elfin[t].append(serialid)
				clusterednetworks[serialid]={"notes":notes}
				networkdata={key:val for key,val in [("network",finalnet),("co-orientations",directed_net),("p_crosssites",pred_gsites),("p_recsites",pred_recsites),("all_broken",int(allbroken)),("p_nearterminal",nearterminal),("p_crossandnt",crossandnt),("cross_seq",seq)]}
				statsdata={key:val for key,val in [("overlaps",overlapstat),("coverages",covstat),("local",localstat),("primary_ct",totalfirst),("secondary_ct",totalseconds),("tertiary_ct",totalthirds),("uniq_contigs",log),("accounted_overoutliers",accounted_overlap_outliers)]}
				clusterednetworks[serialid]["network"]=networkdata
				clusterednetworks[serialid]["stats"]=statsdata
				clusterednetworks[serialid]["model"]={"decision_tree":pos}
				clusterednetworks[serialid]["model"]["prediction"]=pos and percentage>50
				clusterednetworks[serialid]["model"]["probability"]=percentage
				clusterednetworks[serialid]["model"]["modeltype"]=modeltype
				clusterednetworks[serialid]["network"]["networkid"]=serialid
	return
def write_juxtas_txt(elfin):
	with open(networkout,'w') as f:
		json.dump(clusterednetworks, f, indent=4)
	with open(juxtasout, 'w') as out:
		header = ['dnaLo', 'CoordLo', 'dnaHi', 'CoordHi', 'config', 'call', 'overlap', 'origin',
				  'idLo', 'idHi', 'read_termini', 'dna_termini', 'local_count', 'local_rank', 'sampleRead',
				  'overlapSeq', 'readseq', 'seqLo', 'seqHi', 'count', 'sites', 'best_local','Mob', 'network']
		out.write('#' + '\t'.join(header) + '\n')
		for label in sorted(juxtas.keys(), key=lambda x: juxtas[x]['ct'], reverse=True):
			line = juxtas[label]
			line["mobile"],line["network"]=".","."
			if label in elfin:
				line["network"]=','.join(elfin[label])
			elif 'IntID' in line:
				tmob=defaultdict(set)
				for mob in line['IntID']:
					if len(line['IntID'][mob])<=1:
						tmob[mob]=[]
					else:
						tmob['='.join([f"{key}.{val}" for d in line['IntID'][mob] for key,val in d.items()][1:])].add(mob)
				proxmobs=[f"{key}={','.join(list(vals))}" if vals else f"{key}" for key,vals in tmob.items()]
				line["mobile"] = f"{';'.join(proxmobs)}"
			fields = [
				line['hit1dna'], line['coord'], line['hit2dna'], line['coordp'], line['config'], line['call'],
				line['overlap'], line['origin'], line['hit1id'], line['hit2id'], line['read_termini'], line['dna_termini'],
				line['local_count'], line['local_rank'], line['read'], line['jxn'], line['readseq'], line['hitseq1'],
				line['hitseq2'], line['ct'], line['sites'], line['best_local'],line["mobile"],
				line["network"]
				]
			out.write('\t'.join(map(str, fields)) + '\n')
		out.write("# ok\n")
	return


if __name__ == "__main__":
	elfin=defaultdict(list)
	cts={"PhageIntegrase":Counter(),"Recombinase":Counter(),"Resolvase":Counter(),"non-mob":Counter()} 
	negatives,bycross=defaultdict(dict),defaultdict(dict)
	questionable=defaultdict(list)
	coverages=Counter()
	tiles,multiplicity,overlaps,seqs,dnas,mobsLoc,story,calls,clusterednetworks,mapback,director,replacelabs,losingsites,sumpid,dnatermd={},{},{},{},{},{},{},{},{},{},{},{},{},{},{}
	attlines=[]
	pruned=tot=0
	load_ref_files()
	with open(jux_file,'rb') as f:
		juxtas=pickle.load(f)
	read_juxtas()
	write_juxtas_txt(elfin)



