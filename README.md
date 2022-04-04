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
mkdir Neotoma_Transcriptomics
cd Neotoma_Transcriptomics
mkdir raw
mkdir trimmed
mkdir alignments
mkdir htseq
```

### Read Trimming and Filtering
```
cat list.txt | parallel -j 12 "cd {}; trim_galore --cores 4 --paired {}_1.fq.gz {}_2.fq.gz &> {}_tg.log &"
```
