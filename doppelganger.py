#!/usr/bin/env python3

from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import glob
import os
import sys
from collections import OrderedDict
import shutil
wd = os.getcwd() #TODO path problem

def sam_parser(bed_file, out_file, bam_file, threads):
	global sv
	sv = []
	with open(bed_file) as ji:
		for j in ji:
			j = j.rstrip().split('\t')
			a = int(j[1])-3
			b = int(j[2])+3
			bed = j[0]+":"+str(a)+"-"+str(b)
			out = out_file+'.sam'
			os.system('samtools view {} {} -@ {} >> {}.sam'.format(bam_file, bed, threads, out_file))
			sv.append(j[0]+"\t"+j[1])

def sam_regions(out_file, region): #TODO pysam 
	with open(wd+"/"+out_file+".sam") as ij:
		for i in ij:
			i = i.rstrip().split('\t')
			for s in sv:
				s = s.rstrip().split('\t')
				if len(i[9]) > 1 and (int(s[1])-int(region)) > int(i[3]):
					s1 = (int(s[1])-int(region))-int(i[3])
					s2 = (int(s[1])+int(region))-int(i[3])
					with open(wd+"/"+out_file+"_pre-tsd-sequnces.fa", 'a') as k:
						k.write(">"+s[0]+":"+str(s[1])+"-"+str(s[1])+"\n"+i[9][s1:s2]+'\n')					

def composite(out_file):
	global tmp_path_mv 
	for fa in SeqIO.parse(wd+"/"+out_file+"_pre-tsd-sequnces.fa", "fasta"):
		name = fa.id; seq = fa.seq
		if seq != '':
			with open(wd+"/"+out_file+"_pre-tsd-sequnces-composite.fa", 'a') as k:
				k.write(">{}\n{}\n".format(str(name),str(seq)))

	for record in SeqIO.parse(wd+"/"+out_file+"_pre-tsd-sequnces-composite.fa", "fasta"):
		ids = record.id.split()[0]
		ids_split = ids.split(':')
		with open(f'{ids_split[0]+"-"+ids_split[-1]}'+out_file+'.pre.fa', 'a') as f:
			SeqIO.write(record,f,"fasta")
	
	os.remove(out_file+"_pre-tsd-sequnces.fa")
	pre_cons = out_file+"_pre-tsd-sequnces-composite.fa"
	tmp_path = os.system('mkdir {}.tmp'.format(out_file)); tmp_path_mv = str(wd)+"/"+out_file+".tmp"

	os.system('mv *'+out_file+'.pre.fa {}'.format(tmp_path_mv))
	os.remove(out_file+".sam")
	os.remove(out_file+"_pre-tsd-sequnces-composite.fa")

def ins_fasta(out_file, region, fasta):
	global ins_path_mv
	trunked = glob.glob(tmp_path_mv+"/*.fa")
	for t in trunked:
		for fa in SeqIO.parse(t, "fasta"):
			name = fa.id; seq = fa.seq; name1 = name.split("-")
			seq0 = seq[:int(region)-1]; seq1 = seq[int(region)+1:]
			try:
				alignments = pairwise2.align.globalxx(seq0, seq1)
				for a in alignments:
					if int(a[2]) >= int(region)*0.8 or int(a[2]) >= float(region)*0.8: #80% similarity 
						with open(wd+"/"+out_file+"_TSD-prediction.tmp.txt",'a') as pt:
							pt.write(str(name1[0])+"|"+str(seq0) +"|"+str(seq1)+"|"+str(int(a[2]))+'\n')
						
						for ins in SeqIO.parse(fasta, "fasta"): 
							ins_name = ins.id; seq = ins.seq
							ins_id = ins.id.split('_')
							ins_name = ins.id; seq = ins.seq; ins_ids = ins.id.split()[0]; ins_split = ins_ids.split(':')
							if fa.id == ins_id[0]:
								with open(f'{ins_split[0]+"-"+ins_split[-1]}'+out_file+'.INS.fa', 'w') as ins_cond:
									ins_cond.write(">"+str(ins_name)+"\n"+str(ins.seq)+"\n")

			except ValueError:
				print("Continuing... ")
				pass

	ins_path = os.system('mkdir {}.INS'.format(out_file)); ins_path_mv = str(wd)+"/"+out_file+".INS"
	os.system('mv *'+out_file+'.INS.fa {}'.format(ins_path_mv))
	shutil.rmtree(out_file+'.tmp')

# blast sequences of ORF1/2
def blastn_results(out_file, threads, org):
	global blast_out_mv
	ins_condidates = glob.glob(ins_path_mv+"/*.fa")
	for ins_file in ins_condidates:
		iname = ins_file.split()
		iname1 = iname[0].split("/")
		iname2 = iname1[-1].split(".")

		os.system('blastn -db {} -query {} -out {} -evalue 20 -word_size 5 -window_size 0 -num_threads {} -outfmt "6 qseqid sseqid qlen slen qseq sseq sstart send qstart qend evalue length pident nident mismatch gaps"'.format(wd+"/repeats/polyA_tail.txt", ins_file, iname2[0]+'_PA.'+out_file+'.out', threads))

		if org == 'mouse':
			os.system('blastn -db {} -query {} -out {} -evalue 20 -word_size 5 -window_size 0 -num_threads 6 -outfmt "6 qseqid sseqid qlen slen qseq sseq sstart send qstart qend evalue length pident nident mismatch gaps"'.format(wd+"/repeats/mouse.Line1.ORF1", ins_file, iname2[0]+'_ORF1.'+out_file+'.out', threads))
			os.system('blastn -db {} -query {} -out {} -evalue 20 -word_size 5 -window_size 0 -num_threads 6 -outfmt "6 qseqid sseqid qlen slen qseq sseq sstart send qstart qend evalue length pident nident mismatch gaps"'.format(wd+"/repeats/mouse.Line1.ORF2", ins_file, iname2[0]+'_ORF2.'+out_file+'.out', threads))
			os.system('blastn -db {} -query {} -out {} -evalue 20 -word_size 5 -window_size 0 -num_threads 6 -outfmt "6 qseqid sseqid qlen slen qseq sseq sstart send qstart qend evalue length pident nident mismatch gaps"'.format(wd+"/repeats/mouse.Line1.5UTR", ins_file, iname2[0]+'_5UTR.'+out_file+'.out', threads))
			os.system('blastn -db {} -query {} -out {} -evalue 20 -word_size 5 -window_size 0 -num_threads 6 -outfmt "6 qseqid sseqid qlen slen qseq sseq sstart send qstart qend evalue length pident nident mismatch gaps"'.format(wd+"/repeats/mouse.Line1.3UTR", ins_file, iname2[0]+'_3UTR.'+out_file+'.out', threads))
		if org == 'human':
			os.system('blastn -db {} -query {} -out {} -evalue 20 -word_size 5 -window_size 0 -num_threads 6 -outfmt "6 qseqid sseqid qlen slen qseq sseq sstart send qstart qend evalue length pident nident mismatch gaps"'.format(wd+"/repeats/human.LINE1.ORF1", ins_file, iname2[0]+'_ORF1.'+out_file+'.out', threads))
			os.system('blastn -db {} -query {} -out {} -evalue 20 -word_size 5 -window_size 0 -num_threads 6 -outfmt "6 qseqid sseqid qlen slen qseq sseq sstart send qstart qend evalue length pident nident mismatch gaps"'.format(wd+"/repeats/human.LINE1.ORF2", ins_file, iname2[0]+'_ORF2.'+out_file+'.out', threads))
			os.system('blastn -db {} -query {} -out {} -evalue 20 -word_size 5 -window_size 0 -num_threads 6 -outfmt "6 qseqid sseqid qlen slen qseq sseq sstart send qstart qend evalue length pident nident mismatch gaps"'.format(wd+"/repeats/human.Line1.5UTR", ins_file, iname2[0]+'_5UTR.'+out_file+'.out', threads))
			os.system('blastn -db {} -query {} -out {} -evalue 20 -word_size 5 -window_size 0 -num_threads 6 -outfmt "6 qseqid sseqid qlen slen qseq sseq sstart send qstart qend evalue length pident nident mismatch gaps"'.format(wd+"/repeats/human.Line1.3UTR", ins_file, iname2[0]+'_3UTR.'+out_file+'.out', threads))
	
	os.system("find {}/*out -size 0 -delete".format(wd))
	blast_out = os.system('mkdir {}.blast'.format(out_file)); blast_out_mv = str(wd)+"/"+out_file+".blast"
	os.system('mv *'+out_file+'.out {}'.format(blast_out_mv))
	shutil.rmtree(out_file+'.INS')

	os.system(""" for i in """ +blast_out_mv+"""/*out; do awk '{split($1,a,"-") split(a[1], b, ":"); print a[1]"\t"b[1]":"b[2]+$9"-"b[2]+$10"\t"int(($12/$4)*100)"%""\t"$12"\t"$5"\t"$2}' $i | sort | uniq > ${i%.out*}'.sum'; done """)	
	
	with open(wd+"/"+out_file+"_TSD-prediction.tmp.txt") as tmpfile:
		for f in OrderedDict.fromkeys(tmpfile):
	 		with open(wd+"/"+out_file+"_TSD-prediction.txt", 'a') as tsfile:
	 			tsfile.write(f)
	
	os.remove(wd+"/"+out_file+"_TSD-prediction.tmp.txt")

def blast_fmt(out_file):
	
	blast_condidates = glob.glob(blast_out_mv+"/*.sum")
	for blastn in blast_condidates:
		with open(blastn) as blasts:
			for blast in blasts:
				bl = blast.rstrip().split('\t')
				with open(wd+"/"+out_file+"_TSD-prediction.txt") as tsdp:
					for tsd in tsdp:
						ts = tsd.rstrip().split("|")
						if str(ts[0]) == str(bl[0]):
							#print(ts[0])
							with open(wd+"/"+out_file+"_TSD-PA-L1-predictions.tmp.txt", 'a') as tpl:
								tpl.write("|".join((ts))+"|"+str(bl[1])+"|"+str(bl[3])+"|"+str(bl[4])+"|"+str(bl[5])+"\n")

	os.remove(wd+"/"+out_file+"_TSD-prediction.txt") 
	shutil.rmtree(out_file+".blast")
	
	final_file = str(wd)+"/"+out_file+"_TSD-PA-L1-predictions.tmp.txt"
	os.system('sort {} | uniq > {}"_TSD-PA-L1-predictions.txt"'.format(final_file, out_file))
	os.remove(wd+"/"+out_file+"_TSD-PA-L1-predictions.tmp.txt") 
	
	#filter the final file to produce human readable file with SV position showing TSD, internal structure of inserted seq
	final_file_sum = str(wd)+"/"+out_file+"_TSD-PA-L1-predictions.txt"
	os.system(""" awk '{split($0, a, "|"); print a[1]"\t"a[8]}' """ +final_file_sum+ """ | sort | uniq | awk -v OFS='\t' '{ line[$1] = (line [$1] ? line [$1] ";" $2 : $0)} END { for (l in line) print line[l]}' > """ +out_file+"""'_TSD-PA-L1-predictions.summary'""")

	print("Step - 6: Final predictions of TSD, poly(A)tail and L1 repeats are completed for "+out_file)
	print("For more information see "+out_file+"_TSD-PA-L1-predictions.txt and "+out_file+"_TSD-PA-L1-predictions.summary"+"\n")

def doppelganger(bam_file, bed_file, out_file, region, org, fasta, threads):

	print("\n"+"\033[3mdoppelganger\033[0m:"+"\n"+"\t"+"The pipeline can predict target site duplication (TSD), "+"\n"+"poly(A)tail as well as non-reference L1 elements from whole genome sequencing;")
	print("Created and maintained by Ahmed Arslan - aarslan@sbpdiscovery.org"+"\n")

	try:
		sam_parser(bed_file, out_file, bam_file, threads)
		print("Step - 1: Parsing input bam file and selecting reads... ")

	except (IOError, IndexError):
		print("please provide bam file...")

	try:
		sam_regions(out_file, region)
		print("Step - 2: Formating reads data from previous step... ")
	except (IOError, IndexError):
		pass

	try:
		composite(out_file)
		print("Step - 3: Finding TSD in reads... ")
	except (IOError, IndexError):
		pass

	try:
		ins_fasta(out_file, region, fasta)
		print("Step - 4: Finding Poly(A)tail & L1 signals in TSD containing reads... ")
	except (IOError, IndexError):
		pass

	try:
		blastn_results(out_file, threads, org)
		print("Step - 5: Finding L1 genes (ORF1/ORF2) and UTRs for TSD containing reads... ")
	except (IOError, IndexError):
		pass

	try:
		blast_fmt(out_file)
	except (IOError, IndexError):
		pass		


# split input bed
# grep -v 'GL\|JH' TSM.filter.INS.size.bed > TSM.filter.INS.size.GL.bed
# sh ./bedSplit.sh /Users/fortunecookie/Desktop/data/doppelganger/TSM.filter.INS.size.GL.bed
# split fasta
# sh fastSplit.sh /Users/fortunecookie/Desktop/data/doppelganger/TSM.filter.INS.fa
# bam split
# bamtools split -in TSM_all_HiFiCCS.sorted.bam -reference
# rm *JH* && rm *GL* && rm *un*
# for i in *bam; do samtools index $i;done


if __name__ == '__main__':

	import argparse
	from multiprocessing import Pool

	#sys.setrecursionlimit(2000)

	parser = argparse.ArgumentParser(description='doppelganger: use multiprocessing to analyze long-reads and predict non-ref. TSD, PAT, and L1.')
	parser.add_argument('--bam', help='sorted bam file of ref. aligned long-reads')
	parser.add_argument('--bed', help='bed file containing insertion positions [chr\tstart\tend]')
	parser.add_argument('--out', help='name of the output file')
	parser.add_argument('--region', type=int, help='specify region around an insertion in bp [default: 20]')
	parser.add_argument('--fasta', help='fasta file containing the insertion sequences')
	parser.add_argument('--org', help='name of the organism [use: human or mouse]')
	parser.add_argument('-v', '--version', action='version', version='doppelganger v1.0.0')
	parser.add_argument('--threads', type=int, help='specify number of threads')

	args = parser.parse_args()
	p = Pool(int(args.threads))
    
	try:

		p.apply_async(doppelganger, args=(args.bam, args.bed, args.out, args.region, args.org, args.fasta, args.threads))
		p.close()
		p.join()

	except (IOError, IndexError):
		print("To run program properly, follow instructions... ")

	