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
* [HISAT2 Alignment](#hisat2-alignment)
* [HTSeq Counts of Mapping Reads](#htseq-counts-of-mapping-reads)
* [DESeq2 Analysis of Differential Expression using R](#deseq2-analysis-of-differential-expression-using-r)


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

With the HISAT2 alignments in hand, we can now get counts of reads aligning to each feature (i.e. each gene). 


### HTSeq Counts of Mapping Reads

We will use the Neotoma bryanti contig-level assembly and corresponding annotation file found at: https://osf.io/bsd4h/

The gff3 annotation file required some clean up for proper use with HTSeq

```
#removing duplicate exon lines present in a small number of genes
awk '!a[$9]++' Neotoma_bryanti.Contigs.gff3 > Neotoma_bryanti.Contigs.dedup.gff3

#edited three Sult1a1 genes whose stop_codon line was not fully withiin an exon
#the used a gffutils script to create a shared "gene_id" and "gene_name" value among gene, exon, and mrna lines
python gffutils_edit_exon.py
python gffutils_edit_gene.py
python gffutils_edit_mrna.py
python gffutils_edit_exon_Name.py
rm Neotoma_bryanti.Contigs.dedup.gff3
rm tmp.gff; rm tmp2.gff; rm tmp3.gff
mv tmp4.gff Neotoma_bryanti.Contigs.dedup.geneid.gff3

#make tig names match our genome alignments by removing "Nbry_" prefix on fasta name lines
sed -i "s/Nbry_tig/tig/g" Neotoma_bryanti.Contigs.dedup.geneid.fixup.gff3
```

We need to install the latest HTSeq (v2.0.1)

```
pip install htseq
pip install --upgrade numpy
```

Next, we count the reads aligned at each exon of gene by calling the following slurm array scripts:
```
bash htseq_submit_array.sh

#the above shell scripts calls repeated instances of htseq_array.sh, which submits the following for each sample's alignment:

htseq-count --format bam \
--stranded=reverse \
--type=exon \
--idattr=gene_id \
--add-chromosome-info \
--additional-attr=gene_name \
--order=pos \
--mode=union \
--nonunique=random \
--secondary-alignments=ignore \
--samout=/data/gpfs/assoc/matocqlab/Neotoma_transcriptomics/htseq/${FILENAME}_htseq/${FILENAME}.sam \
-c /data/gpfs/assoc/matocqlab/Neotoma_transcriptomics/htseq/${FILENAME}_htseq/${FILENAME}.counts.tsv \
-n 8 \
/data/gpfs/assoc/matocqlab/Neotoma_transcriptomics/alignments/${FILENAME}_align/${FILENAME}_align_sorted.bam \
/data/gpfs/assoc/matocqlab/Neotoma_transcriptomics/Neotoma_bryanti.Contigs.dedup.geneid.fixup.gff3

```
The output of htseq-count output for each sample is a tab-separated file containing ID and gene_id and gene_name of each gene, and an integer value for the read counts aligning to all exons. 

Prepare a directory with only these counts tsv files and export it to your local machine for DESeq2 analysis:

```
#Pull the htseq-count *tsv output into a single directory for export and read in R
mkdir htseq/Counts_files; cd htseq
for dir in $(ls -d W*); do cp $dir/*tsv Counts_files/; done
zip -r counts.zip Counts_files

#From a terminal open to your local machine run:
cd ~/Desktop
mkdir DESeq2_Neotoma_transcriptomics
scp -r mholding@pronghorn.rc.unr.edu:/data/gpfs/assoc/matocqlab/Neotoma_transcriptomics/htseq/counts.zip ~/Desktop/DESeq2_Neotoma_transcriptomics
```

### DESeq2 Analysis of Differential Expression using R

With the counts files on your personal computer, R quickly processes the differential expression analysis. 

Set up directory structure for analysis in R
```
cd ~/Desktop/DESeq2_Neotoma_transcriptomics/
unzip counts.zip
mkdir DEGs; mkdir Meta_data; mkdir Plots
cd Counts_files
mkdir Cecum; mkdir Foregut; mkdir Liver; mkdir SmallIntestine
cd ../
```

Place the tsv counts files into appropriate subdirectories within Counts_files/ for each tissues type, such as "Cecum", "Liver", etc.
Place meta-data files containing sample information such as Species and Diet into the folder Meta_data/

Now you can start RStudio and run the following:

```
#load R packages for DEG data analysis
library(DESeq2)
library(GenomicFeatures)
library(apeglm)
library(ashr)
library(IHW)
library(PoiClaClu)
library(limma)


#load R packages for data visualization
library(RColorBrewer)
library(pheatmap)
library(factoextra)
library(plotly)
library(stats)
library(MASS)
library(reshape2)
library(viridis)
library(ggdendro)
library(gridExtra)
library(gtable)
library(grid)
library(EnhancedVolcano)
library(expss)
library(ggplot2)



setwd("~/Desktop")
dir.create("DESeq2_Neotoma_transcriptomics")

#set  working directory and subdirectory with counts
setwd("~/Desktop/DESeq2_Neotoma_transcriptomics/")
dir.create("DEGs")
dir.create("Plots")


#Choose which tissue to analyze
directory <- "Counts_files/Cecum/"
#directory <- "Counts_files/Foregut/"
#directory <- "Counts_files/Liver/"
#directory <- "Counts_files/SmallIntestine/"

#get  names of htseq-count output files
files <-  list.files(directory, pattern = ".tsv")

#sample metadata loaded here
#choose proper match to input samples
sampleTable <- read.table("Meta_data/cecum_meta.txt", header = TRUE)

#build group variable to facilitate contrasts (rather than specifying interactions model)
#approach detailed in https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#interactions
sampleTable$group <- factor(paste0(sampleTable$Species,"_",sampleTable$Diet))

#build dds object for DESeq2, design is a linear model formula given variables to test
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = directory,
                                  design= ~ 0 + group)

#Run the differential expression analysis on the dataset
dds <- DESeq(dds)
resultsNames(dds)

#read in table of Name and gene_id for results export
genes <- read.table("gene_names.txt", header = TRUE, sep = "\t",
                    row.names = 1)


#set up contrasts of interest
#first create appropriate design matrix:
design <- model.matrix(~0 + group, data = sampleTable)
colnames(design) <- levels(sampleTable$group)

#compare the average effect of being a bryanti to the average effect of being lepida
con.species<- makeContrasts(lepVSbry = (N_lepida_PRFA + N_lepida_RHCA)/2
                             - (N_bryanti_PRFA + N_bryanti_RHCA)/2,
                             levels=design)

#compare the average effect of eating RHCA to the average effect of eating PRFA
con.diet <- makeContrasts(PRFAvsRHCA = (N_lepida_PRFA + N_bryanti_PRFA)/2
                               - (N_lepida_RHCA + N_bryanti_RHCA)/2,
                               levels=design)

#Is the effect of species different depending on diet consumed
con.interaction <- makeContrasts(Interaction = ((N_bryanti_PRFA - N_lepida_PRFA)
                                                - (N_bryanti_RHCA - N_lepida_RHCA)),
                                 levels=design)



#Get table of average effect of species
#Using FDR cutoff of p < 0.05 and required log2fold change above 2
species_res <- results(dds, contrast=con.species, alpha = 0.05)
summary(species_res)
write.table(merge(as.data.frame(subset(species_res, species_res$padj < 0.05 & abs(species_res$log2FoldChange) > 2)), 
                  genes, by="row.names", all.y=FALSE),
            file = "DEGs/species_res.txt",
            quote = FALSE, sep = "\t", row.names=FALSE)

#plot most significant gene for species:  
ix = which.min(species_res$padj) # most significant
barplot(assay(dds)[ix,],las=2, main=rownames(dds)[ ix  ],col =as.factor(dds$group)  )
plotCounts(dds, gene=rownames(dds[ix,]), las=2, intgroup=c("group"), main=rownames(dds[ix,]), )


#The average effect of diet 
diet_res <- results(dds, contrast = con.diet , alpha = 0.05)
summary(diet_res)
write.table(merge(as.data.frame(subset(diet_res, diet_res$padj < 0.05 & abs(diet_res$log2FoldChange) > 2)), 
                  genes, by="row.names", all.y=FALSE),
            file = "DEGs/diet_res.txt",
            quote = FALSE, sep = "\t", row.names=FALSE)

#plot most significant gene for diet:  
ix = which.min(diet_res$padj) # most significant
barplot(assay(dds)[ix,],las=2, main=rownames(dds)[ ix  ],col =as.factor(dds$group)  )
plotCounts(dds, gene=rownames(dds[ix,]), las=2, intgroup=c("group"), main=rownames(dds[ix,]), )

#The genes that show a significant diet x species interaction
#i.e. "Is the diet effect different across species?"
interaction_result <- results(dds, contrast = con.Interaction, alpha = 0.05)
summary(interaction_result)
write.table(merge(as.data.frame(subset(interaction_result, interaction_result$padj < 0.05 & abs(interaction_result$log2FoldChange) > 2)), 
                  genes, by="row.names", all.y=FALSE),
            file = "DEGs/interaction_result.txt",
            quote = FALSE, sep = "\t", row.names=FALSE)


#plot most significant gene for interaction:  
ix = which.min(interaction_result$padj) # most significant
barplot(assay(dds)[ix,],las=2, main=rownames(dds)[ ix  ],col =as.factor(dds$group)  )
plotCounts(dds, gene=rownames(dds[ix,]), las=2, intgroup=c("group"), main=rownames(dds[ix,]), )


#Plot samples in cluster analysis heatmap
vsd <- vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste( vsd$Species, vsd$Diet, sep="--")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


#Plot samples in PCA space with important genes as loading vectors
#Want to use fviz_pca_biplot() instead of default ggplot from DESeq2
# perform a PCA on the data in assay(x) for the selected genes
CEN = scale(assay(vsd), center = T, scale = T)
pca <- prcomp(t(CEN))
loadings <- as.data.frame(pca$rotation)
# visualize
fviz_pca_biplot(pca, repel = TRUE, axes = c(1,2),
                select.var = list(contrib = 25), #draw top 100 arrows
                #select.var = list(name = c("Q375E", "Q375P")),  #alternative to draw specific substitution loadings
                addEllipses = TRUE,
                habillage = dds$group,
                col.ind = dds$Species,
                ellipse.level=0.95,
                geom=c("point"), pointsize = 3.5,   #change to geom=c("point","text") for sample ID
                ind.shape = dds$Diet,
                ind.fill = dds$Species,
                invisible = c( "quali"), #remove enlarged symbol for group mean
                title = "Woodrat Cecum DEG PCA")

#Make volcano plots for each effect
#example below looks at diet effect in N. lepida
diet_volc <- merge(as.data.frame(diet_res), 
      genes, by="row.names", all.y=FALSE)
EnhancedVolcano(diet_volc,
                lab = diet_volc$gene_id,
                title = "Volcano plot",
                subtitle = bquote(italic("Average Diet Effect")),
                x = 'log2FoldChange',
                y = 'padj',
                #xlim = c(-5.5, 5.5),
                #ylim = c(0, -log10(10e-12)),
                pCutoff = 0.05,
                FCcutoff = 2,
                labSize = 2.0,
                col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                drawConnectors = TRUE,
                widthConnectors = 0.25,
)
```


The lists now in the DEGs/ directory contain the gene_ids (last column) that can be input into the ShinyGO web application for GO analyses of enrichment
http://bioinformatics.sdstate.edu/go/
The file gene_names.txt contains the full gene list from N. bryanti to use as background for the GO enrichment analyses.



