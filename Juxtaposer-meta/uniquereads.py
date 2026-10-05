import os,sys,re
from collections import defaultdict

inputs=sys.argv
usage=f"python {inputs[0]} ReadsFile.fasta"

if len(inputs)!=2:
	print(f"Usage:{usage}")
	sys.exit()
else:
	fasta=inputs[-1]
dir,file=os.path.split(fasta)
uniqreads=defaultdict(list)

def readfile(file,obj_type=open):
	with obj_type(file,'r') as f:
		for line in f:
			yield line.strip()

def uniquer():
	count=0
	reads=readfile(fasta)
	defline=next(reads)
	while len(defline)>0:
		try:
			count+=1
			seq=next(reads)
			uniqreads[seq].append(defline)
			defline=next(reads)
		except StopIteration:
			break
	return count
count=uniquer()
ufasta=os.path.join(dir,"ucandidates.fa")
results=os.path.join(dir,"ucandidates.results")

with open(ufasta,'w') as f1, open(results,'w') as f2:
	_ = f1.write('\n'.join( (f"{d[0]}\n{seq}" for seq,d in uniqreads.items()) ) )
	_ = f2.write(f"#Found {len(uniqreads)} unique reads from the original {count} reads\n") 
	_ = f2.write(f"#unique_id\tidentical_reads\n")
	_ = f2.write('\n'.join((f"{d[0]}\t"+','.join(d[1:]) for d in uniqreads.values() if len(d)>1)))
