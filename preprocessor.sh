#!/usr/bin/env python3
#this simple pipeline takes raw pacbio reads from SQII and produce gencDNA signal data
#PS: change tool according to the machine

#conda create -n doppelganger -y -c bioconda lima, ccs, pdindex, bam2fastq, minimap2, samtools, sniffles, bcftools
#conda install -c bioconda lima, ccs, pdindex, bam2fastq, minimap2, samtools, sniffles, bcftools
#conda activate doppelganger

raw=$1
barcodes=$2
threads=$3
out=$4
reference=$5

echo "starting pipeline with removing the barcodes... "
lima --same --split-bam-named $raw $barcodes $out'.bam' -j $threads
echo "producing ccs-hifi reads... "
ccs m64098_220807_202538.subreads.bc2049.bc2049--bc2049.bam $out'.ccs.bam' -j $threads
echo "hifi reads are produced... "
pbindex $out'.ccs.bam'
echo "converting bam2fastq for alignment... "
bam2fastq -o $out $out'.ccs.bam'
echo "starting double alignment with vulcan... "
minimap2 -ax map-hifi -MD $reference $out'.fastq.gz' > $out'.sam' -t $threads
echo "alignment is done... "
samtools view -S -b $out'.sam' -@ $threads > $out'.bam' 
samtools sort -o $out'.sorted.bam' -O bam -@ $threads $out'.bam'
samtools index $out'.sorted.bam' 
echo "starting structural variant calling... "
sniffles --input $out'.sorted.bam' --vcf $out'.vcf' --threads $threads --reference $reference --non-germline --minsupport 1
echo "postprocessing vcf..."
bcftools view -i 'GT="alt"' -f PASS -c 1 $out'.vcf'  -o $out'.filter.vcf' 
bcftools query -f "%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%INFO/SVLEN\t%REF\t%ALT\t%STRAND\n" $out'.filter.vcf' > $out'.filter.strand.bed' 
grep 'INS' $out'.filter.strand.bed' > $out'.filter.strand.INS.bed' 
awk '{print ">"$1":"$2"-"$3"_"$4"_"$5"_"$6"\n"$8}' $out'.filter.strand.INS.bed' > $out'.fa'
awk '{print ">"$1":"$2"-"$3"_"$4"_"$5"_"$6"}' $out'.filter.strand.INS.bed' > $out'.bed'
rm $out'.filter.strand.INS.bed' $out'.filter.strand.bed' $out'.bam' $out'.filter.vcf' 
echo 'data files for doppelganger are ready... '


# preprocessor.sh [fastq] [barcode.fa] [threads] [out_name] [reference]
