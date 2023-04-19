# this R code has been condensed from the tutorial here: 
	# github.com/CBC-UCONN/RNA-seq-with-reference-genome-and-annotation

# Load the libraries we'll need in the following code:
library("DESeq2")
library("apeglm")
library("pheatmap")
library("tidyverse")
library("ggrepel")
library("ashr")
library("goseq")
library("biomaRt")

######################################################
# Point R to count data, set an output file prefix 
######################################################

# create an object with the directory containing your counts:
	# !!edit this to point to your own count file directory!!
directory <- "../04_counts/counts"

# ensure the count files are where you think they are
list.files(directory)

# create a list of sample files
sampleFiles <- list.files(directory, pattern = ".*counts$")

# sample files don't perfectly match metadata table. need to split names and clean up so we can associate metadata with count files
sampleFiles <- data.frame(sampleFiles=sampleFiles) %>%
  separate(sampleFiles,into=c("athena","year","rna","x1","x2","x3","suffix"),sep="_") %>%
  cbind(.,sampleFiles)
  # replace SUB03 with S03
  sampleFiles$x2 <- str_replace(sampleFiles$x2,"SUB03","S03")

# there is a sample swap. here I'm going to switch the count files to fix it. 
  # These are the count files to be swapped
  # [1] "Athena1_2021_RNA_POST_S15_POSTHA1_S81.counts"
  # [2] "Athena1_2022_RNA_POST_S12_POSTHA1_S75.counts"
sampleFiles[c(30,63),8] <- sampleFiles[c(63,30),8] 


######################################################
# Read the count data into R along with treatment information
######################################################

# load in metadata table
  # read in, split column 1 so we can match up with the sample files
meta <- read.table("../metadata/metadata.txt", header=TRUE, sep="\t") 
  colnames(meta) <- c("tube","sampleData","timePoint","sex","subjectID","group")
  meta <- separate(meta,tube,into=c("athena","year","rna","x1","x2","x3"),sep=" ")

# join meta and sampleFiles tables
metaFiles <- left_join(meta,sampleFiles[,-c(1:3,7)])

# now create a data frame with sample names, file names and treatment information. 
sampleTable <- data.frame(
	sampleName = paste(metaFiles$x1,metaFiles$x2,metaFiles$x3,sep="_"),
	fileName = metaFiles$sampleFiles,
	timepoint = as.factor(metaFiles$timePoint),
  sex = as.factor(metaFiles$sex),
	subject = as.factor(metaFiles$subjectID),
  subject.n = as.factor(c(rep(1,14),rep(2,14),rep(3,14),rep(1,14),rep(2,14),rep(3,14)))
	)

# look at the data frame to ensure it is what you expect:
sampleTable

# create the DESeq data object
  # THE DESIGN FORMULA HERE IS A DUMMY FORMULA. THE DESIGN NEEDS TO BE UPDATED NEXT
    # subject sex is a "linear combination" of subject IDs, leading to a matrix that is not "full rank"
    # per this section of the DESeq2 vignette: https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#group-specific-condition-effects-individuals-nested-within-groups

ddsHTSeq <- DESeqDataSetFromHTSeqCount(
		sampleTable = sampleTable, 
		directory = directory, 
		design = ~ timepoint + sex
		)



# set the design formula properly
# add the dummy subject variable
ddsHTSeq$subject.n <- sampleTable$subject.n
# update the design formula
design(ddsHTSeq) <- ~ sex + sex:subject.n + sex:timepoint


######################################################
# Filter out genes with very low expression
######################################################

# what does expression look like across genes?

# sum counts for each gene across samples
sumcounts <- rowSums(counts(ddsHTSeq))
# take the log
logsumcounts <- log(sumcounts,base=10)
# plot a histogram of the log scaled counts
hist(logsumcounts,breaks=100)

# get genes with summed counts greater than 20
keep <- sumcounts > 252

# keep only the genes for which the vector "keep" is TRUE
ddsHTSeq <- ddsHTSeq[keep,]

######################################################
# Run the statistical analysis
######################################################

dds <- DESeq(ddsHTSeq)


######################################################
# Which results do we want?
######################################################

# in this experiment there are 32 coefficents. you can see them here:

resultsNames(dds)

# we have a huge list of contrasts to do

# first compare males vs females, look for interactions with time point.-----------------------------------------------
  # set up contrast lists for deseq2
  # the first element gives the male/female difference at timepoint 1, the rest give male/female differences that are *significantly different* from timepoint 1
MvFxTime <- list()
    MvFxTime[[1]] <- list("sex_Male_vs_Female")

for(i in 2:14){
    MvFxTime[[i]]<- list(paste("sexMale.timepoint",i,sep=""),paste("sexFemale.timepoint",i,sep=""))
    }

# extract each contrast, put the table in a list
MvFxTime_res <- list()

for(i in 1:14){
    MvFxTime_res[[i]] <- results(dds, contrast=MvFxTime[[i]])
    print(i)
    }

names(MvFxTime_res) <- c("Males-vs-Females-1",paste("Males-vs-Females_",2:14,"vs",1,sep=""))

# now timepoint comparisons within and averaged across sexes, and interactions(?)---------------------------------------

# set up numerator and denominator of timepoints to be contrasted for log-fold changes (i.e. timepoint 2 vs timepoint 1)
num <- c(2,3,6,8,10,12,4,5,5,11,12,10,10,6,8,10,12,13,14)

den <- c(1,1,1,1,1, 1, 3,3,4,10,10,6, 8, 3,3,3, 3, 12,12)


# first, comparisons within males---------------------------------------------------------------------------------------
M_time <- list()

for(i in 1:19){
    if(den[i] == 1){
        M_time[[i]] <- list(paste("sexMale.timepoint",num[i],sep=""))  
        }
    if(den[i] > 1) {
        M_time[[i]] <- list(paste("sexMale.timepoint",num[i],sep=""),paste("sexMale.timepoint",den[i],sep=""))
        }
    }

# extract each contrast, put the table in a list
M_time_res <- list()

for(i in 1:19){
    M_time_res[[i]] <- results(dds, contrast=M_time[[i]])
    print(i)
    }

names(M_time_res) <- paste("Males_", num, "vs", den, sep="")

# next, comparisons within females---------------------------------------------------------------------------------------
F_time <- list()

for(i in 1:19){
    if(den[i] == 1){
        F_time[[i]] <- list(paste("sexFemale.timepoint",num[i],sep=""))  
        }
    if(den[i] > 1) {
        F_time[[i]] <- list(paste("sexFemale.timepoint",num[i],sep=""),paste("sexFemale.timepoint",den[i],sep=""))
        }
    }

# extract each contrast, put the table in a list
F_time_res <- list()

for(i in 1:19){
    F_time_res[[i]] <- results(dds, contrast=F_time[[i]])
    print(i)
    }

names(F_time_res) <- paste("Females_", num, "vs", den, sep="")

# next, averaged across sexes---------------------------------------------------------------------------------------------
MF_time <- list()

for(i in 1:19){
    if(den[i] == 1){
        MF_time[[i]] <- list(c(paste("sexMale.timepoint",num[i],sep=""),paste("sexFemale.timepoint",num[i],sep="")))  
        }
    if(den[i] > 1) {
        MF_time[[i]] <- list(c(paste("sexMale.timepoint",num[i],sep=""),paste("sexFemale.timepoint",num[i],sep="")),c(paste("sexMale.timepoint",den[i],sep=""),paste("sexFemale.timepoint",den[i],sep="")))
        }
    }

# extract each contrast, put the table in a list
MF_time_res <- list()

for(i in 1:19){
    MF_time_res[[i]] <- results(dds, contrast=MF_time[[i]])
    print(i)
    }

names(MF_time_res) <- paste("Males+Females_", num, "vs", den, sep="")



# finally, sex x time interactions (some are redundant with the MvFxTime object ---------------------------------------------------------------------------------------------
MFxtime <- list()

for(i in 1:19){
    if(den[i] == 1){
        MFxtime[[i]] <- list(paste("sexMale.timepoint",num[i],sep=""),paste("sexFemale.timepoint",num[i],sep=""))  
        }
    if(den[i] > 1) {
        MFxtime[[i]] <- list(c(paste("sexMale.timepoint",num[i],sep=""),paste("sexFemale.timepoint",den[i],sep="")),c(paste("sexMale.timepoint",den[i],sep=""),paste("sexFemale.timepoint",num[i],sep="")))
        }
    }

# extract each contrast, put the table in a list
MFxtime_res <- list()

for(i in 1:19){
    MFxtime_res[[i]] <- results(dds, contrast=MFxtime[[i]])
    print(i)
    }

names(MFxtime_res) <- paste("Males-vs-Females_", num, "vs", den, sep="")



# combine results into one gigantic table, long format

# initialize table
bigTable <- c()

# Females + Males
for(i in 1:19){
    temp <- data.frame(contrast=names(MF_time_res)[i],geneID=rownames(MF_time_res[[i]]),as.data.frame(MF_time_res[[i]]))
    bigTable <- rbind(bigTable, temp)
    print(i)
}

# Females
for(i in 1:19){
    temp <- data.frame(contrast=names(F_time_res)[i],geneID=rownames(F_time_res[[i]]),as.data.frame(F_time_res[[i]]))
    bigTable <- rbind(bigTable, temp)
    print(i)
}

# Males
for(i in 1:19){
    temp <- data.frame(contrast=names(M_time_res)[i],geneID=rownames(M_time_res[[i]]),as.data.frame(M_time_res[[i]]))
    bigTable <- rbind(bigTable, temp)
    print(i)
}

# Females vs Males part 1
for(i in 1:14){
    temp <- data.frame(contrast=names(MvFxTime_res)[i],geneID=rownames(MvFxTime_res[[i]]),as.data.frame(MvFxTime_res[[i]]))
    bigTable <- rbind(bigTable, temp)
    print(i)
}


# Females vs Males part 2, include only unique comparisons not already in part 1. that's indexes 7-19
for(i in 7:19){
    temp <- data.frame(contrast=names(MFxtime_res)[i],geneID=rownames(MFxtime_res[[i]]),as.data.frame(MFxtime_res[[i]]))
    bigTable <- rbind(bigTable, temp)
    print(i)
}


# move the table to wide format
bigTableWide <- pivot_wider(
    bigTable,
    id_cols=c("geneID", "baseMean"),
    names_from="contrast",
    values_from=c("log2FoldChange","pvalue","padj")
    )

# reorder the columns. there's got to be a better way

nameorder <- c("geneID","baseMean")
for(i in 1:19){
    nameorder <- c(nameorder, paste(c("log2FoldChange","pvalue","padj"), rep(names(MF_time_res)[i],3),sep="_"))
    nameorder <- c(nameorder, paste(c("log2FoldChange","pvalue","padj"), rep(names(F_time_res)[i],3),sep="_"))
    nameorder <- c(nameorder, paste(c("log2FoldChange","pvalue","padj"), rep(names(M_time_res)[i],3),sep="_"))
    }
for(i in 1:14){
    nameorder <- c(nameorder, paste(c("log2FoldChange","pvalue","padj"), rep(names(MvFxTime_res)[i],3),sep="_"))
    }
for(i in 7:19){
    nameorder <- c(nameorder, paste(c("log2FoldChange","pvalue","padj"), rep(names(MFxtime_res)[i],3),sep="_"))
    }

bigTableWide <- as.data.frame(bigTableWide)[,nameorder]


######################################################
# Quickly summarize results by contrast type
######################################################

bigSummary <- bigTable %>%
    group_by(contrast) %>%
    summarize(significant=sum(padj < 0.1, na.rm=TRUE),Up=sum(log2FoldChange > 0 & padj < 0.1, na.rm=TRUE),Down=sum(log2FoldChange < 0 & padj < 0.1, na.rm=TRUE)) %>%
    data.frame()


######################################################
# Data visualization
######################################################


#############

# PCA plot

# normalized, variance-stabilized transformed counts for visualization
vsd <- vst(dds, blind=FALSE)

plotPCA(vsd, intgroup=c("timepoint","sex"))

# alternatively, using ggplot

dat <- plotPCA(vsd,returnData=TRUE,intgroup=c("timepoint","sex"))

p <- ggplot(dat,aes(x=PC1,y=PC2,col=timepoint,shape=sex))
p <- p + geom_point() + 
	xlab(paste("PC1: ", round(attr(dat,"percentVar")[1],2)*100, "% variation explained", sep="")) + 
	ylab(paste("PC2: ", round(attr(dat,"percentVar")[2],2)*100, "% variation explained", sep="")) +
	geom_label_repel(aes(label=name))
p




##############

# plot counts for individual genes

plotCounts(dds,"ENSG00000153234",intgroup=c("timepoint","sex"),col=factor(dds$sex),pch=as.numeric(dds$subject))

######################################################
# Get gene annotations using biomaRt
######################################################


##############################
# Select a mart and dataset
##############################

# ensembl host:
  # most recent is "https://ensembl.org"
  # to list archived version hosts: listEnsemblArchives()

ensemblhost <- "https://dec2021.archive.ensembl.org"

listMarts(host=ensemblhost)

# create an object for the Ensembl Genes v100 mart
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host=ensemblhost)

# occasionally ensembl will have connectivity issues. we can try an alternative function:
	# select a mirror: 'www', 'uswest', 'useast', 'asia'
	# mart <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", mirror = "useast")

# see a list of datasets within the mart
	# at the time of writing, there were 203
listDatasets(mart)

# figure out which dataset is the killifish
	# be careful using grep like this. verify the match is what you want
searchDatasets(mart,pattern="Human")

# there's only one match, get the name
genedata <- searchDatasets(mart,pattern="Human")[,1]

# create an object for the dataset
gene_mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = ensemblhost, dataset = genedata)


#########################
# Query the mart/dataset
#########################

# filters, attributes and values

# see a list of all "filters" available for the dataset.
	# at the time of writing, over 300
listFilters(gene_mart)

# see a list of all "attributes" available
	# 129 available at the time of writing
listAttributes(mart = gene_mart, page="feature_page")

# we can also search the attributes and filters
searchAttributes(mart = gene_mart, pattern = "ensembl_gene_id")

searchFilters(mart = gene_mart, pattern="ensembl")

# get gene names and transcript lengths when they exist
ann <- getBM(filter="ensembl_gene_id",value=rownames(res1),attributes=c("ensembl_gene_id","description","transcript_length"),mart=gene_mart)

# pick only the longest transcript for each gene ID
ann <- group_by(ann, ensembl_gene_id) %>% 
  summarize(.,description=unique(description),transcript_length=max(transcript_length)) %>%
  as.data.frame()





######################################################
# Write out some results
######################################################

write.csv(bigTableWide,file="")

write.csv(resdata, file = paste0(outputPrefix, "-results-with-normalized.csv"))

# send normalized counts to tab delimited file for GSEA, etc.
write.table(
  as.data.frame(counts(dds),normalized=T), 
  file = paste0(outputPrefix, "_normalized_counts.txt"), 
  sep = '\t'
)

