if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("dada2")


#Set the path
path <- getwd() #sets the path to the files
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#Quality assessment 
af<-plotQualityProfile(fnFs, aggregate=TRUE)
aR<-plotQualityProfile(fnRs, aggregate=TRUE)


#Filter and Trim
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

#Learn Errors
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)
#Denoise
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

#Merge the reads and create a sequnence table
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)


seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

#Track the reads
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#Add the taxonomy

wget "https://zenodo.org/record/1172783/files/silva_nr_v132_train_set.fa.gz?download=1"

taxa<-assignTaxonomy(seqtab.nochim, "/content/silva_nr_v132_train_set.fa.gz?download=1", multithread=TRUE)
colnames(taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

#Align the sequences and create a phylogenetic tree
BiocManager::install("phyloseq")
BiocManager::install("DECIPHER")
BiocManager::install("phangorn")

seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)

phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeUPGMA  <- upgma(dm)

#Create a phyloseq object
ps <- phyloseq(tax_table(taxa), sample_data(metadata),
                 otu_table(seqtab.nochim,taxa_are_rows=FALSE),phy_tree(treeUPGMA))


#Read sums table
read_sums<-data.frame(rowSums(otu_table(ps)))
colnames(read_sums)
colnames(read_sums)<-c("reads")
read_sums$sample.names<-rownames(read_sums)

#Normalise the data
ps_rare <- phyloseq::rarefy_even_depth(ps,sample.size=4699, rngseed = 123, replace = FALSE)

#Make a counts table without the sequences (they are replaced with ASV)
table<-otu_table(ps)

ncol(table)

ASV<-rep(c("ASV"),c(999))
ASV<-paste0(ASV,"",seq(1:999))

sequences<-colnames(table)
sequences<-data.frame(sequences, ASV)

colnames(table)<-ASV

counts_table<-table[,-51]

#Rarefaction plots
rare_data<-rarecurve(t(counts_table), step = 20, sample=5000, col = "blue", cex = 0.6, tidy=TRUE)

#Alpha diversity
rich<-estimate_richness(ps_rare, split = TRUE, measures = NULL)
rich

#Beta diversity
beta_dist = phyloseq::distance(ps_rare, method="unifrac", weighted=F)
ordination = ordinate(ps_rare, method="PCoA", distance=beta_dist)
ordination_data<-ordination$vectors
ordination_Eigenvalues<-ordination$values$Eigenvalues

#All sample code
!unzip /content/rarefied_table.qza
all_samples_table<-read.delim("table.from_biom.txt", skip = 1)
rownames(all_samples_table)<-all_samples_table[,1]
all_samples_table<-all_samples_table[,-1]
all_metadata<-read.csv("all_samples_metadata.csv", header=T)
rownames(all_metadata)<-all_metadata[,1]
all_metadata<-all_metadata[,-1]
ps_all <- phyloseq(sample_data(all_metadata),
                 otu_table(all_samples_table, taxa_are_rows = TRUE))

rich_all<-estimate_richness(ps_all, split = TRUE, measures = NULL)




#Beta diversity
beta_dist = phyloseq::distance(ps_all, method="bray")
ordination = ordinate(ps_all, method="PCoA", distance=beta_dist)
ordination_data<-ordination$vectors
ordination_Eigenvalues<-ordination$values$Eigenvalues
