# Analysis of differential gene expression in woodrats

![Gene Expression](/gene_plus_RNA.png "Annotation/Alignment Image")

Software required:

* [GNU Parallel](https://www.gnu.org/software/parallel/)
* [Trim Galore!](https://github.com/FelixKrueger/TrimGalore)
* [HISAT2](http://daehwankimlab.github.io/hisat2/manual/)
* [samtools](http://www.htslib.org/)
* [HTSeq](https://htseq.readthedocs.io/en/master/)
* R
* RStudio
* [DESeq2](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)


## Contents

* [Build Directory Structure](#build-directory-structure)
* [Read Trimming and Filtering](#read-trimming-and-filtering)


### Build Directory Structure 
Set up working directory. Place raw fastq files containing RNA-seq paireed in reads in a subfolder called "raw"

```
mkdir Neotoma_transcriptomics
cd Neotoma_transcriptomics
mkdir raw
mkdir trimmed
mkdir alignments
mkdir htseq
```

### Read Trimming and Filtering
We will impose these filters to trim reads:

* Remove 5' end bases if phred quality is below 20
* Remove 3' end bases if phred quality is below 20
* Remove Illumina adapter strings
* Minimum read length of each read in a pair = 20bp
```
cat sample_names.txt | parallel -j 6 "~/software/TrimGalore-0.6.5/trim_galore \
--cores 4 -o trimmed \
--paired raw/{}_R1_001.fastq.gz raw/{}_R2_001.fastq.gz &> trimmed/{}_tg.log &"
```


### HISAT2 Alignment
Build a HISAT2 database from the Neotoma bryanti genome

```
cd genomes
hisat2-build -p 23 Neotoma_bryanti.Arrow.fasta Nbryanti_db &
```



