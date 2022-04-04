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
cd ../
```

Align paired reads with splice-site-aware HISAT2 aligner

Running /Shell_scripts/rna_submitarray.sh on a Slurm system will create one job with one 8 core cpu for each sample in sample_names.txt

```
bash rna_submitarray.sh
```
The above shell scripts will execute HISAT2 alignments and subsequent conversion to a sorted and indexed bam file (using samtools) as follows:

```
cd alignments
mkdir ${FILENAME}_align
cd ${FILENAME}_align
hisat2 -q --phred33 --no-temp-splicesite --no-mixed --no-discordant --max-intronlen 150000 \
--rna-strandness RF --no-unal -p 8 \
-x /data/gpfs/assoc/matocqlab/Neotoma_transcriptomics/genome/Nbryanti_db \
-1 /data/gpfs/assoc/matocqlab/Neotoma_transcriptomics/trimmed/${FILENAME}_R1_001_val_1.fq.gz \
-2 /data/gpfs/assoc/matocqlab/Neotoma_transcriptomics/trimmed/${FILENAME}_R2_001_val_2.fq.gz \
-S ${FILENAME}_align.sam

/data/gpfs/home/mholding/software/samtools-1.11/samtools view -bh -@ 8 -f 3 -F 256 -q 40 ${FILENAME}_align.sam | \
/data/gpfs/home/mholding/software/samtools-1.11/samtools sort -@ 8 > ${FILENAME}_align_sorted.bam
/data/gpfs/home/mholding/software/samtools-1.11/samtools index -b ${FILENAME}_align_sorted.bam
#remove large, unnecessary sam file
rm ${FILENAME}_align.sam
```










