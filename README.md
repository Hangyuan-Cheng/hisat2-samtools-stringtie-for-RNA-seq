# hisat2-samtools-stringtie-for-RNA-seq
This repository will guide you complete the analysis of RNA-seq offline using the following software:hisat2, samtools, and stringtie. Double-terminal sequencing data has been filtered by FastQC, and the specific filtering process is not involved in this paper. The sample data is for reference only. For the guided in Chinese version, you can visit my CSDN: https://blog.csdn.net/m0_64240043/article/details/122450444?spm=1001.2014.3001.5501

Let's look at the data to be processed as a whole
```shell
@biocloud:~/1223/NGS2022$ tree
.
├── gene_data.csv # final result!!
├── genome
│   ├── yeast.fa# reference genome
│   ├── yeast.gff# genome annotation for reference
│   └── yeast_transcriptome.fa # The transcriptome index file required by Kallisto, which is used to achieve transcription quantification based on pseudo Alignment, is not covered in this article
├── reads	
│   ├── s1_y_1.fq.gz
│   ├── s1_y_2.fq.gz
│   ├── s2_y_1.fq.gz
│   └── s2_y_2.fq.gz
├── script	 # Script file. If necessary, add an absolute path to reference the script.
│   ├── edgeR.R
│   ├── prepDE.py
│   ├── prepDE.py3
│   └── run.sh
└── src		# A software installation package that may be used. If the software has been installed on the server, you do not need to install it again.
    ├── fastqc_v0.11.9.zip
    ├── kallisto_linux-v0.46.1.tar.gz
    ├── samtools-1.14.tar.bz2
    ├── stringtie-2.2.0.Linux_x86_64.tar.gz
```
--------
## 1.hisat2 builds a reference genome index
```shell
$ cd ./genome 
$ hisat2-build yeast.fa genome.fa 
```
## 2.hisat2 compared double-ended sequencing reads under the same treatment to the reference genome
```shell
$ cd ..
$ mkdir alignment
$ cd alignment

$ hisat2 -p 6 -x ../genome/genome.fa -1 ../reads/s1_y_1.fq.gz -2 ../reads/s1_y_2.fq.gz -S ./s1.sam
```
## 3.samtools：binary transforms, sorts, and indexes Sam files.
```shell
$ samtools view -bS s1.sam -o s1.bam 
$ samtools sort s1.bam s1.sorted 
$ samtools index s1.sorted.bam 
```
## 4.stringtie estimates transcript expression based on GTF annotation and bam alignment after sequencing
```shell
$ stringtie ./s1.sorted.bam -G ../genome/yeast.gff -e -p 2 -o ./s1_out.gtf -A ./s1_genes.list 
```
The comparison results of s2_out. GTF and S2_genes. List were obtained by the same method. The steps are exactly the same as above and will not be repeated here. A pair of double-ended sequencing results under each process produces an S1_out.gtf file. After following the steps above to process another pair of double-endsequencing, you should have two GTF files: s1_out.gtf and s2_out.gtf.

## 5.The transcriptome expression levels under the two treatments were summarized
```shell
$ cd ..
$ mkdir differential_expression
$ cd differential_expression

$ vim sample_list.txt 
```

## 6.The prepde.py3 script for Stringtie is used to generate a list of differentially expressed genes.
prepDE.py3 could be download from github（[https://github.com/gpertea/stringtie](https://github.com/gpertea/stringtie)）
```shell
$ python prepDE.py3 -i sample_list.txt -g gene_count.csv -t transcript.csv
```
The gene_count.csv obtained after the above steps is the list of differentially expressed genes.
