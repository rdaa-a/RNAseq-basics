## Author: Rodrigo Bedera Garcia
## This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0$

library(ballgown)
library(FactoMineR)
library(factoextra)
library(limma)
library(DESeq2)


## Ballgown needs a file with the folder names of each sample and the group
## they belong to.
pheno.data <- read.csv("pheno_data.csv")

## Load samples
bg.data <- ballgown(dataDir = "samples/", samplePattern = "sample", pData=pheno.data[c(1,10,11,12,2,3,4,5,6,7,8,9),])
bg.data
sampleNames(bg.data)

## Extract FPKM and correct sample order
gene.expression <- gexpr(bg.data)
head(gene.expression)
dim(gene.expression)
gene.expression <- gene.expression[,c(1,5,6,7,8,9,10,11,12,2,3,4)]
gene.names <- rownames(gene.expression)
## Naming columns
colnames(gene.expression) <- c("wtll_1","wtll_2","wtll_3","wthl_1","wthl_2","wthl_3",
                               "vipll_1","vipll_2","vipll_3","viphl_1","viphl_2","viphl_3")
head(gene.expression)

# Maximum orf detection across samples
max(c(sum(gene.expression[,1]>0),sum(gene.expression[,2]>0),sum(gene.expression[,3]>0),sum(gene.expression[,4]>0),
       sum(gene.expression[,5]>0),sum(gene.expression[,6]>0),sum(gene.expression[,7]>0),sum(gene.expression[,8]>0),
       sum(gene.expression[,9]>0),sum(gene.expression[,10]>0),sum(gene.expression[,11]>0),sum(gene.expression[,12]>0)))

## PCA analysis ############
set.seed(123)
colnames(gene.expression) <- c("wt.ll_1","wt.ll_2","wt.ll_3","wt.hl_1","wt.hl_2","wt.hl_3",
                               "vip1-1.ll_1","vip1-1.ll_2","vip1-1.ll_3","vip1-1.hl_1","vip1-1.hl_2","vip1-1.hl_3")

pca.gene.expression <- data.frame(colnames(gene.expression),t(gene.expression))
colnames(pca.gene.expression)[1] <- "condition"


res.pca <- PCA(pca.gene.expression, graph = FALSE,scale.unit = TRUE,quali.sup = 1)

fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 70),main = "")

fviz_pca_ind(res.pca, col.ind = factor(c(rep("WT LL",3),rep("WT HL",3),rep("vip1-1 LL",3),rep("vip1-1 HL",3)),levels=c("WT LL","WT HL","vip1-1 LL","vip1-1 HL")),
             pointsize=2, pointshape=21,fill="black",
             repel = TRUE,
             addEllipses = TRUE,ellipse.type = "confidence",
             legend.title="Conditions", legend="bottom",
             title="Principal Components Analysis",
             show_legend=TRUE,show_guide=TRUE) +
  theme(plot.title=element_text(face="bold",hjust=0.5))





## Normalization: upper quantile
upper.quantiles <- vector(mode="numeric",length=ncol(gene.expression))

for(i in 1:ncol(gene.expression))
{
  upper.quantiles[i] <- quantile(gene.expression[,i],probs=0.75)
}

mean.upper.quantiles <- mean(upper.quantiles)

for(i in 1:ncol(gene.expression))
{
  gene.expression[,i] <- (gene.expression[,i] / upper.quantiles[i]) * mean.upper.quantiles
}

## Log2 transformation
log.gene.expression <- log2(gene.expression+1)

## Alternative normalization
# library(NormalyzerDE)
# design <- data.frame(sample=colnames(gene.expression),
#                      group=c(rep("wtll",3),rep("wthl",3),rep("vipll",3),rep("viphl",3)))
# 
# write.table(x = design,file = "normalyzer_design.tsv",quote = F,row.names = F,
#             sep = "\t")
# gene.expression.1 <- gene.expression +1
# write.table(x = gene.expression.1,file = "vip_gene_expression.tsv",
#             quote = F,row.names = F,
#             sep = "\t")
# 
# normalyzer(jobName = "vip",designPath = "normalyzer_design.tsv",
#            dataPath = "vip_gene_expression.tsv",outputDir = ".")
# 
# 
# normalized.gene.expression <- read.table(file="vip/Quantile-normalized.txt", header=T)
# normalized.gene.expression[is.na(normalized.gene.expression)] <- 0
# 
# head(normalized.gene.expression)
# rownames(normalized.gene.expression) <- gene.names
# nrow(normalized.gene.expression)




## DEG analysis using limma, a separate analysis for each comparison

## vip LL vs wt LL: ##

normalized.wtvsvipll <- log.gene.expression[,c(1,2,3,7,8,9)]

## Experimental design

experimental.design <- model.matrix(~ -1+factor(c(1,1,1,2,2,2)))
colnames(experimental.design) <- c("wtll","vipll")

## Adjustment of each gene expression to a linear model

linear.fit <- lmFit(normalized.wtvsvipll, experimental.design)

## Comparisons to perform

contrast.matrix <- makeContrasts(vipll-wtll,levels=c("wtll","vipll"))

## Fold-change and p-value calculation

contrast.linear.fit <- contrasts.fit(linear.fit, contrast.matrix)
contrast.results <- eBayes(contrast.linear.fit)
nrow(normalized.wtvsvipll)

wtvip.ll <- topTable(contrast.results, number=17741,coef=1,sort.by="logFC")
head(wtvip.ll)

## Give names to objects

log.fold.change <- wtvip.ll$logFC
q.value <- wtvip.ll$adj.P.Val
genes.ids <- rownames(wtvip.ll)
names(log.fold.change) <- genes.ids
names(q.value) <- genes.ids

## Selection of DEGs. Here, fold change of 2 and adjusted p-value 
## of 0.05 were set as threshold

activated.genes <- genes.ids[log.fold.change > log(base = 2,x = 2) & q.value < 0.05]
repressed.genes <- genes.ids[log.fold.change < - log(base = 2,x = 2) & q.value < 0.05]


length(activated.genes)
length(repressed.genes)

write(x=activated.genes,file="activated_genes_ll.txt",sep="\t")
write(x=repressed.genes,file="repressed_genes_ll.txt",sep="\t")

## Volcano plot
log.p.val <- -log10(wtvip.ll$P.Value)
names(log.p.val) <- rownames(wtvip.ll)

plot(log.fold.change,log.p.val,pch=19,col="grey",cex=0.8,
     xlim=c(-3,3),ylim = c(0,8), 
     xlab="log2(Fold-change)",ylab="-log10(p-value)",cex.lab=1.5,
     main="WT vs vip1-1 LL")
lines(x=c(-1,-1),y=c(-10,10),lty=2)
lines(x=c(1,1),y=c(-10,10),lty=2)
lines(x=c(-10,10),y=c(-log10(0.05),-log10(0.05)),lty=2)
points(x = log.fold.change[activated.genes],
       y = log.p.val[activated.genes],col="red",cex=0.8,pch=19)
points(x = log.fold.change[repressed.genes],
       y = log.p.val[repressed.genes],col="blue",cex=0.8,pch=19)
text(x=log.fold.change[activated.genes],2.4,labels = activated.genes,cex=0.7)


## vip HL vs wt HL: same proceedings ##

normalized.wthlvsviphl <- log.gene.expression[,c(4,5,6,10,11,12)]

experimental.design <- model.matrix(~ -1+factor(c(1,1,1,2,2,2)))
colnames(experimental.design) <- c("wthl","viphl")

linear.fit <- lmFit(normalized.wthlvsviphl, experimental.design)

contrast.matrix <- makeContrasts(viphl-wthl,levels=c("wthl","viphl"))

contrast.linear.fit <- contrasts.fit(linear.fit, contrast.matrix)
contrast.results <- eBayes(contrast.linear.fit)

nrow(normalized.wthlvsviphl)

wtvip.hl <- topTable(contrast.results, number=17741,coef=1,sort.by="logFC")
head(wtvip.hl)

log.fold.change <- wtvip.hl$logFC
q.value <- wtvip.hl$adj.P.Val
genes.ids <- rownames(wtvip.hl)
names(log.fold.change) <- genes.ids
names(q.value) <- genes.ids

activated.genes <- genes.ids[log.fold.change > log(x = 2,base = 2) & wtvip.hl$adj.P.Val < 0.05]
repressed.genes <- genes.ids[log.fold.change < - log(x=2,base=2) & wtvip.hl$adj.P.Val < 0.05]
#write(x=art.repressed.genes,file = "repressed_genes.txt",sep = "\t")

length(activated.genes)
length(repressed.genes)

write(x=activated.genes,file="activated_genes_hl.txt",sep="\t")
write(x=repressed.genes,file="repressed_genes.txt",sep="\t")

## Volcano plot
log.p.val <- -log10(wtvip.hl$P.Value)
names(log.p.val) <- rownames(wtvip.hl)

plot(log.fold.change,log.p.val,pch=19,col="grey",cex=0.8,
     xlim=c(-3,3),ylim = c(0,8), 
     xlab="log2(Fold-change)",ylab="-log10(p-value)",cex.lab=1.5,
     main="WT vs vip1-1 HL")
lines(x=c(-1,-1),y=c(-10,10),lty=2)
lines(x=c(1,1),y=c(-10,10),lty=2)
lines(x=c(-10,10),y=c(-log10(0.05),-log10(0.05)),lty=2)
points(x = log.fold.change[activated.genes],
       y = log.p.val[activated.genes],col="red",cex=0.8,pch=19)
points(x = log.fold.change[repressed.genes],
       y = log.p.val[repressed.genes],col="blue",cex=0.8,pch=19)
text(x=log.fold.change[activated.genes],2.4,labels = activated.genes,cex=0.7)


## wt HL vs wt LL: same proceedings ##

normalized.wtllhl <- log.gene.expression[,c(1,2,3,4,5,6)]

experimental.design <- model.matrix(~ -1+factor(c(1,1,1,2,2,2)))
colnames(experimental.design) <- c("wtll","wthl")

linear.fit <- lmFit(normalized.wtllhl, experimental.design)

contrast.matrix <- makeContrasts(wthl-wtll,levels=c("wtll","wthl"))

contrast.linear.fit <- contrasts.fit(linear.fit, contrast.matrix)
contrast.results <- eBayes(contrast.linear.fit)

nrow(normalized.wthlvsviphl)

wt.llhl <- topTable(contrast.results, number=17741,coef=1,sort.by="logFC")
head(wtvip.hl)


## vip HL vs vip LL: same proceedings ##

normalized.vipllhl <- log.gene.expression[,c(7,8,9,10,11,12)]

experimental.design <- model.matrix(~ -1+factor(c(1,1,1,2,2,2)))
colnames(experimental.design) <- c("vipll","viphl")

linear.fit <- lmFit(normalized.vipllhl, experimental.design)

contrast.matrix <- makeContrasts(viphl-vipll,levels=c("vipll","viphl"))

contrast.linear.fit <- contrasts.fit(linear.fit, contrast.matrix)
contrast.results <- eBayes(contrast.linear.fit)

nrow(normalized.wthlvsviphl)

vip.llhl <- topTable(contrast.results, number=17741,coef=1,sort.by="logFC")
head(wtvip.hl)



###########################################################################
## DESeq2 analysis: interaction term ##

## This file is also needed, like in the limma analysis

pheno.data <- read.csv("pheno_data.csv")
pheno.data
pheno.data$genotype <- factor(pheno.data$genotype)
pheno.data$condition <- factor(pheno.data$condition)

## DESeq2 uses counts instead of FPKM
gene.count.matrix <- read.table(file = "gene_count_matrix.csv",header = T,sep = ",",row.names = 1)
head(gene.count.matrix)
gene.count.matrix <- gene.count.matrix[,pheno.data$sample]
gene.ids <- rownames(gene.count.matrix)

## The interaction term analysis is sensitive to very low counts. A gene increasing their 
## expression as a response to high light from 5 to 15 counts will test positive, while
## the same 10 count increase for genes reporting 1000 counts will be insignificant.
## Therefore, it is recommended to filter low counts, as the count detection in these
## genes is proportionally very variable. The mean counts will not be used as DESeq2 input

wtll.cont <- (gene.count.matrix[,1] + gene.count.matrix[,2] + gene.count.matrix[,3]) /3
wthl.cont <- (gene.count.matrix[,4] + gene.count.matrix[,5] + gene.count.matrix[,6]) /3
vipll.cont <- (gene.count.matrix[,7] + gene.count.matrix[,8] + gene.count.matrix[,9]) /3
viphl.cont <- (gene.count.matrix[,10] + gene.count.matrix[,11] + gene.count.matrix[,12]) /3
unified.cont.matrix <- matrix(c(wtll.cont,wthl.cont,vipll.cont,viphl.cont),byrow = F,nrow = length(viphl.cont))

## Filter genes which in no condition have more than 30 mean counts
filtered.gene.count.matrix <- gene.count.matrix[apply(gene.count.matrix,MARGIN = 1,FUN = max) > 30,]
nrow(filtered.gene.count.matrix) #15406
gene.ids <- rownames(filtered.gene.count.matrix)

## Prepare deseq2 object. An interaction term analysis is used, as we are interesting in finding genes 
## with a different response to high light in the mutant, i.e., genes where the expression in mutant HL is not 
## the result of adding the response to HL in wt + effect of genotype (expression change between wt LL and mutant LL).

dds <- DESeqDataSetFromMatrix(countData=filtered.gene.count.matrix, colData=pheno.data, design = ~ genotype + condition + genotype:condition)
dds$condition <- relevel(dds$condition,ref="ll")
dds$genotype <- relevel(dds$genotype,ref="wt")
dds$condition
dds$genotype

## Perform analysis

dds <- DESeq(dds)
resultsNames(dds)

## Explore results

res_lighteffect_to_wt <- results(dds,contrast=c("condition","hl","ll")) 
res_lighteffect_to_vip <- results(dds,list(c("condition_hl_vs_ll","genotypevip.conditionhl")))

## Select desired result

is_diff_lighteffect <- results(dds,name="genotypevip.conditionhl")
is_diff_lighteffect

## Save result in variable

res <- is_diff_lighteffect

log.fold.change <- res$log2FoldChange
q.value <- res$padj
names(log.fold.change) <- gene.ids
names(q.value) <- gene.ids

## Obtain genes with a significantly different response to high light in the mutan strain
## P-value cutoff of 0.01 was selected (this testing is less restricitve than limma testing)
significant.genes.deseq2 <- row.names(res)[q.value < 0.01]
significant.genes.deseq2 <- significant.genes.deseq2[!is.na(significant.genes.deseq2)]

length(significant.genes.deseq2)

write(significant.genes.deseq2,file="diff_sign_lighteffect_filtered_counts_001.txt")



##################################################

#### WORKFLOW: barplots and heatmaps
#### REQUIRES : query_names, gene.expression,
# phytozome query is returned in order, must give sort(query_names) as input 
# Also, if only gene_description is asked, query will remove duplicates. Ask for gene_name1 too

library(biomaRt)
library(httr)
library(jsonlite)
library(xml2)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(plyr)

builder.df.ggplot <- function(expr.matrix,gene,cond.names){
  expr.1 <- unlist(c(expr.matrix[gene, 1:3]))
  expr.2 <- unlist(c(expr.matrix[gene, 4:6]))
  expr.3 <- unlist(c(expr.matrix[gene, 7:9]))
  expr.4 <- unlist(c(expr.matrix[gene, 10:12]))
  
  mean.1 <- mean(expr.1)
  mean.2 <- mean(expr.2)
  mean.3 <- mean(expr.3)
  mean.4 <- mean(expr.4)
  
  sd.1 <- sd(expr.1)
  sd.2 <- sd(expr.2)
  sd.3 <- sd(expr.3)
  sd.4 <- sd(expr.4)
  
  means <- c(mean.1, mean.2,mean.3,mean.4)
  sds <- c(sd.1, sd.2,sd.3,sd.4)
  return(data.frame(cond.names,means,sds))
}
cond.names <- c("WT LL","WT HL","vip1-1 LL","vip1-1 HL")

# Change labs, cond.names' levels
roda.ggplotter <- function(gene.expression,query_names,cond.names,title,titlesize){
  # Start of the function. Variables to give: gene.expression,query_names,cond.names
  myplots <- list()
  for(i in 1:length(query_names)){ 
    #build df for each gene
    ggplot.daata <- builder.df.ggplot(gene.expression,query_names[i],
                                      cond.names =  cond.names)
    
    # Standard deviation of the mean as error bar, parameters explained
    p <- ggplot(ggplot.daata, aes(x=factor(cond.names,levels=c("WT LL","WT HL","vip1-1 LL","vip1-1 HL")),
                                  y=means, fill=cond.names)) + 
      geom_bar(stat="identity", position=position_dodge()) + #creates barplot
      geom_errorbar(aes(ymin=means-sds, ymax=means+sds), width=.2,
                    position=position_dodge(.9)) + # finishes barplot
      xlab("") + ylab("") + ggtitle(as.character(queryresult$gene_description[i])) + # appropiate labs
      guides(fill=guide_legend(title="")) + #legend title
      scale_fill_brewer(palette="Dark2") + # palette, this is from Rcolorbrewer (name is enough)
      theme_minimal() + # minimal theme (remove to have grey background)
      theme(plot.title = element_text(size=titlesize,hjust=0.5,face="bold")) # title font size
    
    myplots[[i]] <- p # adds plot to list
  }
  
  # outside loop, plot in a grid all the plots. Technical difficulties, need to insert in a list all default args.
  # change necessary args
  def.param.list <- list( 
    plotlist = NULL,
    ncol = NULL,
    nrow = NULL,
    labels = NULL,
    label.x = 0,
    label.y = 1,
    hjust = -0.5,
    vjust = 1.5,
    font.label = list(size = 14, color = "black", face = "bold", family = NULL),
    align = c("none", "h", "v", "hv"),
    widths = 1,
    heights = 1,
    legend = "bottom",
    common.legend = T,
    legend.grob = NULL
  )
  
  # Unpack the plots and the arguments and create a figure 
  figure <- splat(ggarrange)(c(myplots,def.param.list))
  annotate_figure(figure,
                  #top = text_grob(as.character(title), color = "black", face = "bold", size = 14),
                  #bottom = text_grob("Conditions", color = "black", face = "bold"),
                  left = text_grob("FPKM", color = "black", rot = 90,face="bold"),
                  #right = "I'm done, thanks :-)!",
                  fig.lab = as.character(title), fig.lab.face = "bold",fig.lab.pos = "top"
  )                 
  # end of function
}

#### USAGE:
query_names <- sort(c("Cre01.g042750", "Cre03.g144807", "Cre03.g149100", "Cre12.g514750"))

# Barplots of selected genes, retrieving gene names from phytozome
mart <- useMart(biomart = 'phytozome_mart_archive', 
                host = "https://phytozome-next.jgi.doe.gov/",
                dataset = 'Creinhardtii_281', 
                port = 443)
queryresult <- getBM(attributes = "gene_description",
                     filters=c("gene_name_filter"),
                     values=query_names,
                     mart=mart,
                     verbose = T)
cond.names <- c("WT LL","WT HL","vip1-1 LL","vip1-1 HL")
roda.ggplotter(gene.expression = gene.expression,query_names = query_names,cond.names = cond.names,"Acetate assimilation",10)

# Heatmap of selected genes
query_norm_expr <- gene.expression[rownames(gene.expression) %in% query_names,]
pheatmap(query_norm_expr,color= rev(brewer.pal(40,"PuOr")),cluster_cols = F,
         scale="row",cluster_rows = T,fontsize_row=8)


