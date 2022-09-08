# doppelganger
Target site duplication assessment from alignment file and genomic location of insertion sites. 


Instructions to produce input files:
  
  1 - bam file (sorted)
  2 - bed file (containing insertion locations)
  3 - fasta file (containing insertion sequences)
  
  
1 - Generating (sorted) BAM file:
    
    Softwares needed:[minimap2](https://github.com/lh3/minimap2),[samtools](http://www.htslib.org/),[reference](https://www.ncbi.nlm.nih.gov/grc/human) 

    Alignment:
    Command: minimap2 -ax map-hifi -MD GRCh38_genomic.fa HiFiCCS.fastq.gz > HiFiCCS.sam -t 10

      A - change map-hifi to map-ont, if using nanopore reads.
      B - change -t to increase or decrease number of threads. 
      C - change reference genome (GRCh38_genomic.fa) to your requirement.

    Convert sam to bam file format:
    Command: samtools view -Sb HiFiCCS.sam -@ 10 > HiFiCCS.bam 
      A - change -@ to increase or decrease number of threads.
    
    Sort bam file:
    Command: samtools sort -o HiFiCCS.sorted.bam -O bam -@ 10 HiFiCCS.bam 
    Index sorted bam file:
    Command: samtools index HiFiCCS.sorted.bam
  
  
2 - Generating bed and fasta from bed file
    
    Softwares needed:[Sniffles](https://github.com/fritzsedlazeck/Sniffles), [bcftools](https://samtools.github.io/bcftools/howtos/install.html)
    
    Step - Insertion variant calling:
    Command: sniffles --input HiFiCCS.sorted.bam --vcf HiFiCCS.sorted.vcf --threads 10 --reference GRCh38_genomic.fa --non-germline --minsupport 1
    
      A - change reference genome (GRCh38_genomic.fa) to your requirement.
      B - change to germline variant prediction by removing option --non-germline
      C - increase number of reads required to predict an insertion by changing option --minsupport [int]
      
    Step - vcf to bed/fasta:
    Command: bcftools view -i 'GT="alt"' -f PASS -c 1 HiFiCCS.sorted.vcf -o HiFiCCS.filter.vcf
    Command: bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%INFO/SVLEN\t%REF\t%ALT\t%STRAND\n' HiFiCCS.filter.vcf > HiFiCCS.filter.bed
    Command: grep 'INS' HiFiCCS.filter.bed > HiFiCCS.filter.INS.bed
    Command: awk '{print ">"$1":"$2"-"$3"_"$4"_"$5"_"$6"\n"$8}' HiFiCCS.filter.INS.bed > HiFiCCS.fa
    Command: awk '{print ">"$1":"$2"-"$3"_"$4"_"$5"_"$6}' HiFiCCS.filter.INS.bed > HiFiCCS.bed
    
        A - HiFiCCS.fa, HiFiCCS.bed and HiFiCCS.sorted.bam are the input files for doppelganger. 
        

    Remove extra/tmp files:
    rm HiFiCCS.sam HiFiCCS.bam HiFiCCS.sorted.vcf HiFiCCS.filter.bed HiFiCCS.filter.INS.bed HiFiCCS.filter.vcf