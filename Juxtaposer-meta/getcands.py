import os,sys,re,gzip,subprocess,concurrent.futures
from collections import defaultdict,Counter
import pysam
import time
from datetime import datetime
#import tracemalloc

#tracemalloc.start()
inputs = sys.argv
if len(inputs)==5:
	aln_file,fqreads,outdir,readlen=inputs[1:]
	readlen=int(readlen)
else:
	print(f"Usage: python {inputs[0]} alignment_file[BAM|SAM] reads(comma separated) output_directory read_length")
	sys.exit()
towrite,files,counts=[],{},Counter()
seqs={}
pre="candidates"
faout,outsum=[os.path.join(outdir,name) for name in f"{pre}.fa {pre}.results".split()]

for file in fqreads.split(',')+[aln_file]:
	test=subprocess.run(["file",file], capture_output=True)
	ftype=True if "gzip" in str(test.stdout) else False
	files[file]=ftype

def readfile(file,zipped,obj_type=open):
	if obj_type==open and zipped:
		obj_type,readtype=gzip.open,'rt'
	else:
		readtype='rb' if zipped else 'r'
	with obj_type(file,mode=readtype) as f:
		try:
			for line in f:
				yield line if obj_type==pysam.AlignmentFile else line.strip()
		except Exception as e:
			print(f"Could not read file {file} due to the following error: {e}")

def processFQ(fq,zipped):
	global s
	print(f"Reading FASTQ file")
	fq_obj=readfile(fq,zipped)
	defline=next(fq_obj)
	utseq,tseqs={},{}
	while len(defline)>0:
		try:
			seq=next(fq_obj)
			blank=next(fq_obj)
			qual=next(fq_obj)
			read_id=defline.replace("@","")
			if seq in utseq:
				utseq[seq].append(read_id)
			else:
				utseq[seq]=[read_id]
				tseqs[read_id]=seq
			defline=next(fq_obj)
		except StopIteration:
			utseq={read[0]:(tseqs[read[0]],None if len(read)==1 else read[1:]) for read in utseq.values()}
			print(f"{time.time()-s} Finished reading FASTQ file")
			return utseq
#def categorize(read,line,read_len):
	#try:
	#	read,cont=[line.query_name,line.reference_name]
	#except Exception  as e:
#		print(e)
#		sys.exit() #raise?
#	if line.is_paired:
#		read=f"{read}/{[1,2][int(line.is_read2)]}"
def categorize(mydict):
	results={}
	m,u="well-mapped","unmapped"
	for read,line in mydict.items():
		try:
			cont = line.reference_name
			if not cont:
				results[read]=u
				continue
				#return "unmapped", read
			asco, bdtype = line.get_tag("AS", True)
		except Exception:
			results[read]=u
			continue
			#raise "Rerun bowtie recommended"
			#return "unmapped", read
		alen=line.query_alignment_length
		if  asco<=((-0.073*alen)-16.6):
			results[read]=u
			continue
			#return "unmapped",read
		results[read]=m
		#return "well-mapped",read
	return results
def catreads(sam,zipped,read_len,n=None):
	counts=Counter()
	global s
	available=os.cpu_count()
	MAX_WORKERS=min([32,(os.cpu_count() or 1)+4])
	with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_WORKERS,thread_name_prefix='CtgrzRds') as executor, open(outsum,'w') as f1, open(faout,'w') as f:
		lastread = None
		tempstore={}
		futures = set()
		cycle=0
		for samline in readfile(sam,zipped,pysam.AlignmentFile):
			if not isinstance(samline,pysam.libcalignedsegment.AlignedSegment):
				print(samline, type(samline))
				continue
				sys.exit() #raise?
			if samline.query_name==lastread:
				continue
			lastread = samline.query_name
			read_id=lastread if not samline.is_paired else f"{lastread}/{[1,2][int(samline.is_read2)]}"
			if not read_id in seqs:
				continue
			#try:

			#	cont,alen,asco=samline.reference_name,samline.query_alignment_length,samline.get_tag("AS",True)[0]
			#except Exception as e:
			#	print(f"Could not process alignment entry due to the following issue.\n{e}")
			#	continue
			tempstore[read_id]=samline
			if len(tempstore)>MAX_WORKERS:
				futures.add(executor.submit(categorize,tempstore))
				tempstore={}
#			futures.add(executor.submit(categorize,read_id,samline,read_len))
			if len(futures) > MAX_WORKERS:
				cycle+=1
				finished, futures = concurrent.futures.wait(futures, return_when=concurrent.futures.FIRST_COMPLETED)
				for completed_future in finished:
#					mapping,read = completed_future.result()
					cats = completed_future.result()
					for read,mapping in cats.items():
						if mapping!="well-mapped":
							#towrite.append(read)
							if read in seqs:
								seq,reads=seqs[read]
								identical=""
								if not reads==None:
#								identical="" if len(shared)<=1 else " "+' '.join(shared)
									identical=f" {' '.join(reads)}"
								_ = f.write(f">{read}{identical}\n{seq}\n")
						counts[mapping]+=1
				#print(f"{time.time()-s} finished cycle {cycle}")
		for future in concurrent.futures.as_completed(futures):		#leftover futures if any
#			mapping,read = future.result()
			cats = future.result()
			for read,mapping in cats.items():
				if mapping!="well-mapped":
					if read in seqs:
					#towrite.append(read)
						seq,reads=seqs[read]
						identical=""
						if not reads==None:
#						identical="" if len(shared)<=1 else " "+' '.join(shared)
							identical=f" {' '.join(reads)}"
						#shared=uniqseq.get(read)
						#identical="" if not shared else " "+' '.join(shared)
						_ = f.write(f">{read}{identical}\n{seq}\n")
				counts[mapping]+=1
		_ = f1.write("\n".join([f"{key}\t{counts[key]}" for key in "unmapped well-mapped".split()+["No data"]]))
		del counts
def processFQs(readsets):
	global seqs
	with concurrent.futures.ProcessPoolExecutor() as executor:
		futures = {executor.submit(processFQ,filename,files[filename]): filename for filename in readsets}
		for future in concurrent.futures.as_completed(futures):
			filename = futures[future]
			try:
				result = future.result()
				seqs.update(result)
			except Exception as e:
				print(f'File {filename} generated an exception: {e}')

if __name__ == '__main__':
	#if cpu==1:
	#	noparallel(aln_file,ftype,readlen)
	#else:
	s=time.time()

	processFQs(fqreads.split(','))
#	seqs=processFQ(fqreads,files[fqreads])
	catreads(aln_file,files[aln_file],readlen)
	print(f"Finished in {time.time()-s}")
#	snapshot=tracemalloc.take_snapshot()
#	top_stats = snapshot.statistics('lineno')
#	for stat in top_stats[:10]:
#		print(stat)
