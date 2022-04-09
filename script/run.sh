#!/bin/bash

while true
do
    echo "$(date) - Start sequencing data analysis"
    fastqc --no-extract -t 4 s1_y_1.fq.gz s1_y_2.fq.gz s2_y_1.fq.gz s2_y_2.fq.gz

    echo "$(date) - Trimming low quality reads"
    trim_galore --paired --stringency 6 --gzip s1_y_1.fq.gz s1_y_2.fq.gz
    trim_galore --paired --stringency 6 --gzip s2_y_1.fq.gz s2_y_2.fq.gz
   
    echo "$(date) - Runnning QC"
    multiqc .

    echo "$(date) - Building Bwa index"
    bwa index yeast.fa

    echo "$(date) - Mapping using BWA-MEM"
    bwa mem yeast.fa s1_y_1_val_1.fq.gz s1_y_2_val_2.fq.gz > s1_bwa.sam
    bwa mem yeast.fa s2_y_1_val_1.fq.gz s2_y_2_val_2.fq.gz > s2_bwa.sam

    echo "$(date) - Building HISAT2 index"
    hisat2-build yeast.fa yeast.fa

    echo "$(date) - HISAT2 alignment"
    hisat2 --dta -x yeast.fa -1 s1_y_1_val_1.fq.gz -2 s1_y_2_val_2.fq.gz > s1_hisat.sam
    hisat2 --dta -x yeast.fa -1 s2_y_1_val_1.fq.gz -2 s2_y_2_val_2.fq.gz > s2_hisat.sam
    
    echo "$(date) - convert sam to bam format"
    samtools view -bS s1_hisat.sam > s1_hisat.unsorted.bam
    samtools view -bS s2_hisat.sam > s2_hisat.unsorted.bam

    echo "$(date) - sort alignment results"
    samtools sort s1_hisat.unsorted.bam -o s1_hisat.sorted.bam
    samtools sort s2_hisat.unsorted.bam -o s2_hisat.sorted.bam

    echo "$(date) - indexing sorted bam"
    samtools index s1_hisat.sorted.bam
    samtools index s2_hisat.sorted.bam

    echo "$(date) - transcript quantification"
    stringtie -e -G yeast.gff -o s1_out.gtf -A s1_genes.list s1_hisat.sorted.bam
    stringtie -e -G yeast.gff -o s2_out.gtf -A s2_genes.list s2_hisat.sorted.bam

    echo "$(date) - Pseudo alignment"
    /histor/public/software/kallisto/kallisto index -i yeast_transcriptome.idx yeast_transcriptome.fa
    /histor/public/software/kallisto/kallisto quant -i yeast_transcriptome.idx -o ./kallisto_s1 s1_y_1.fq.gz s1_y_2.fq.gz
    /histor/public/software/kallisto/kallisto quant -i yeast_transcriptome.idx -o ./kallisto_s2 s2_y_1.fq.gz s2_y_2.fq.gz

    echo $(date)
    sleep 2
done
