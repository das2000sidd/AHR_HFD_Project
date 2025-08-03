setwd("~/Desktop/PhD_Project_related/PARP7_HFD_Study");


library(tximportData);
library(TxDb.Mmusculus.UCSC.mm39.refGene);
library(BSgenome.Mmusculus.UCSC.mm10);
library(ensembldb);
library(EnsDb.Mmusculus.v79);
library(biomaRt);
library(tibble);
library(tximport);
library(DESeq2);
library(edgeR);
library("PoiClaClu");
library("pheatmap");
library("RColorBrewer");
library(org.Mm.eg.db);


res <- read.table(file="GRCm39_gene_and_transcript_stable_ID_version.txt",
	header = T,
	sep="\t",
	stringsAsFactors = F
	);

res <- res[,c(1,4)];

colnames(res) <- c("GENEID","TXNAME");

res <- res[,c(2,1)];

res_tibble <- as_tibble(res);

files <- list.files(path="/Users/siddhaduio.no/Desktop/PhD_Project_related/PARP7_HFD_Study",
                    full.names = T
                    );

files <- files[c(12:59)];

files <- as.data.frame(files);



files_quant <- paste(files$files,
                     "/quant.sf",
                     sep=""
                     );

samp_tab <- read.csv(file="Comparison_sheet.csv",
                     header = T,
                     stringsAsFactors = F
                     );

samp_tab$Comp_var <- paste(samp_tab$Tissue,samp_tab$Diet,
                           samp_tab$Genotype,
                           sep="_"
                           );

samp_tab <- samp_tab[order(samp_tab$Sequence.name),];


sample_order <- c(paste("BAT_Control_PARP7H532A",1:4,sep=""),
                  paste("BAT_Control_WT",1:4,sep=""),
                  paste("BAT_High fat_PARP7H532A",1:4,sep=""),
                  paste("BAT_High fat_WT",1:4,sep=""),
                  paste("Liver_Control_PARP7H532A",1:4,sep=""),
                  paste("Liver_Control_WT",1:4,sep=""),
                  paste("Liver_High fat_PARP7H532A",1:4,sep=""),
                  paste("Liver_High fat_WT",1:4,sep=""),
                  paste("WAT_Control_PARP7H532A",1:4,sep=""),
                  paste("WAT_Control_WT",1:4,sep=""),
                  paste("WAT_High fat_PARP7H532A",1:4,sep=""),
                  paste("WAT_High fat_WT",1:4,sep="")
                  );

txi <- tximport(files_quant,
                type = "salmon",
                tx2gene = res_tibble
                );

names(txi);
head(txi$counts);

colnames(txi$counts) <- sample_order;

head(txi$counts);



model_mat_geno_treat <- model.matrix( ~ 0 + Comp_var,data=samp_tab);

dds <- DESeqDataSetFromTximport(txi, 
                                samp_tab, 
                                model_mat_geno_treat
                                );

ntd <- normTransform(dds);

library("vsn");

meanSdPlot(assay(ntd));

## Exploratory analysis and visualisation
nrow(dds);

## Generating variance stabilised object
vsd <- vst(dds, blind = FALSE);

head(assay(vsd), 3);

colData(vsd);

library("dplyr");
library("ggplot2");

## Estimating size factor for normalisation
dds <- estimateSizeFactors(dds);

###PCA plot***

plotPCA(vsd, 
	intgroup = c("Genotype")) + 
	ggtitle("PCA by genotype");


## DE analysis using DESeq2
library(Glimma)

dds$group <- factor(paste0(dds$Comp_var));

design(dds) <- ~ group;

dds <- DESeq(dds);

plotMDS(dds);

##Calling results without any arguments will extract the estimated log2 fold changes and p 
##values for the last variable in the design formula

contrast_wt_vs_parp7mut_liver <- c("group", 
                                   "Liver_Control_PARP7H532A", 
                                   "Liver_Control_WT"
                                   );

res_wt_vs_parp7mut_liver <- results(dds, 
                                    contrast=contrast_wt_vs_parp7mut_liver,
                                    pAdjustMethod = "BH",
                                    format = "DataFrame"
                                    );

res_wt_vs_parp7mut_liver_df <- as.data.frame(res_wt_vs_parp7mut_liver);

res_wt_vs_parp7mut_liver_df$Ensembl <- rownames(res_wt_vs_parp7mut_liver_df);

summary(res_wt_vs_parp7mut_liver);

res_wt_vs_parp7mut_liver_df <- res_wt_vs_parp7mut_liver_df[complete.cases(res_wt_vs_parp7mut_liver_df),];## 20069


contrast_wt_vs_parp7mut_highfat_liver <- c("group", 
                                           "Liver_High fat_PARP7H532A", 
                                           "Liver_High fat_WT"
                                           );

res_wt_vs_parp7mut_highfat_liver <- results(dds, 
                                            contrast=contrast_wt_vs_parp7mut_highfat_liver,
                                            pAdjustMethod = "BH",
                                            format = "DataFrame"
                                            );

res_wt_vs_parp7mut_highfat_liver_df <- as.data.frame(res_wt_vs_parp7mut_highfat_liver);

res_wt_vs_parp7mut_highfat_liver_df$Ensembl <- rownames(res_wt_vs_parp7mut_highfat_liver_df);

summary(res_wt_vs_parp7mut_highfat_liver);

res_wt_vs_parp7mut_highfat_liver_df <- res_wt_vs_parp7mut_highfat_liver_df[complete.cases(res_wt_vs_parp7mut_highfat_liver_df),]; ## 16224

contrast_wt_vs_parp7mut_BAT <- c("group", 
                                 "BAT_Control_PARP7H532A", 
                                 "BAT_Control_WT"
                                 );

res_wt_vs_parp7mut_BAT <- results(dds, 
                                  contrast=contrast_wt_vs_parp7mut_BAT,
                                  pAdjustMethod = "BH",
                                  format = "DataFrame"
                                  );

res_wt_vs_parp7mut_BAT_df <- as.data.frame(res_wt_vs_parp7mut_BAT);

res_wt_vs_parp7mut_BAT_df$Ensembl <- rownames(res_wt_vs_parp7mut_BAT_df);

summary(res_wt_vs_parp7mut_BAT);

res_wt_vs_parp7mut_BAT_df <- res_wt_vs_parp7mut_BAT_df[complete.cases(res_wt_vs_parp7mut_BAT_df),];


contrast_wt_vs_parp7mut_highfat_BAT <- c("group", 
                                         "BAT_High fat_PARP7H532A", 
                                         "BAT_High fat_WT"
                                         );

res_wt_vs_parp7mut_highfat_BAT <- results(dds, 
                                          contrast=contrast_wt_vs_parp7mut_highfat_BAT,
                                          pAdjustMethod = "BH",
                                          format = "DataFrame"
                                          );
                                          
res_wt_vs_parp7mut_highfat_BAT_df <- as.data.frame(res_wt_vs_parp7mut_highfat_BAT);

res_wt_vs_parp7mut_highfat_BAT_df$Ensembl <- rownames(res_wt_vs_parp7mut_highfat_BAT_df);

summary(res_wt_vs_parp7mut_highfat_BAT);

res_wt_vs_parp7mut_highfat_BAT_df <- res_wt_vs_parp7mut_highfat_BAT_df[complete.cases(res_wt_vs_parp7mut_highfat_BAT_df),];


contrast_wt_vs_parp7mut_WAT <- c("group", 
                                 "WAT_Control_PARP7H532A", 
                                 "WAT_Control_WT"
                                 );

res_wt_vs_parp7mut_WAT <- results(dds, contrast=contrast_wt_vs_parp7mut_WAT,
                                  pAdjustMethod = "BH",
                                  format = "DataFrame"
                                  );

res_wt_vs_parp7mut_WAT_df <- as.data.frame(res_wt_vs_parp7mut_WAT);

res_wt_vs_parp7mut_WAT_df$Ensembl <- rownames(res_wt_vs_parp7mut_WAT_df);

summary(res_wt_vs_parp7mut_WAT);

res_wt_vs_parp7mut_WAT_df <- res_wt_vs_parp7mut_WAT_df[complete.cases(
                            res_wt_vs_parp7mut_WAT_df
                            ),];


contrast_wt_vs_parp7mut_highfat_WAT <- c("group", 
                                         "WAT_High fat_PARP7H532A", 
                                         "WAT_High fat_WT"
                                         );

res_wt_vs_parp7mut_highfat_WAT <- results(dds, 
                                          contrast=contrast_wt_vs_parp7mut_highfat_WAT,
                                          pAdjustMethod = "BH",format = "DataFrame"
                                          );

res_wt_vs_parp7mut_highfat_WAT_df <- as.data.frame(res_wt_vs_parp7mut_highfat_WAT);

res_wt_vs_parp7mut_highfat_WAT_df$Ensembl <- rownames(res_wt_vs_parp7mut_highfat_WAT_df);

summary(res_wt_vs_parp7mut_highfat_WAT);

res_wt_vs_parp7mut_highfat_WAT_df <- res_wt_vs_parp7mut_highfat_WAT_df[complete.cases(res_wt_vs_parp7mut_highfat_WAT_df),];

## Remaining contrasts asked for
contrast_hfd_vs_chow_liver_PARP7H532A <- c("group", 
                                           "Liver_High fat_PARP7H532A", 
                                           "Liver_Control_PARP7H532A"
                                           );

contrast_hfd_vs_chow_WAT_PARP7H532A <- c("group", 
                                         "WAT_High fat_PARP7H532A", 
                                         "WAT_Control_PARP7H532A"
                                         );

contrast_hfd_vs_chow_BAT_PARP7H532A <- c("group", 
                                         "BAT_High fat_PARP7H532A",
                                         "BAT_Control_PARP7H532A"
                                         );

contrast_hfd_vs_chow_liver_WT <- c("group", 
                                   "Liver_High fat_WT", 
                                   "Liver_Control_WT"
                                   );

contrast_hfd_vs_chow_WAT_WT <- c("group", 
                                 "WAT_High fat_WT", 
                                 "WAT_Control_WT"
                                 );

contrast_hfd_vs_chow_BAT_WT <- c("group", 
                                 "BAT_High fat_WT", 
                                 "BAT_Control_WT"
                                 );


res_hfd_vs_chow_BAT_WT <- results(dds, 
                                  contrast=contrast_hfd_vs_chow_BAT_WT,
                                  pAdjustMethod = "BH",
                                  format = "DataFrame"
                                  );

res_hfd_vs_chow_BAT_WT <- as.data.frame(res_hfd_vs_chow_BAT_WT);

res_hfd_vs_chow_BAT_WT$Ensembl <- rownames(res_hfd_vs_chow_BAT_WT);

summary(res_hfd_vs_chow_BAT_WT);

res_hfd_vs_chow_BAT_WT <- res_hfd_vs_chow_BAT_WT[complete.cases(res_hfd_vs_chow_BAT_WT),];

sig_genes <- subset(res_hfd_vs_chow_BAT_WT,
                    abs(res_hfd_vs_chow_BAT_WT$log2FoldChange) > 1 & res_hfd_vs_chow_BAT_WT$padj < 0.05);

sig_genes$Entrez <- mapIds(org.Mm.eg.db, 
                           sig_genes$Ensembl,
                           keytype="ENSEMBL", 
                           column="ENTREZID",
                           multiVals = "first"
                           );

sig_genes$Symbol <- mapIds(org.Mm.eg.db, 
                           sig_genes$Entrez,
                           keytype="ENTREZID", 
                           column="SYMBOL",
                           multiVals = "first"
                           );

sig_genes$Gene_Name <- mapIds(org.Mm.eg.db, 
                              sig_genes$Entrez,
                              keytype="ENTREZID", 
                              column="GENENAME",
                              multiVals = "first"
                              );

sig_genes$Symbol <- as.character(sig_genes$Symbol);

sig_genes$Gene_Name <- as.character(sig_genes$Gene_Name);

write.table(sig_genes,
            file="Sig_genes_HFD_vs_chow_BAT_WT.txt",
            col.names = T,
            row.names = F,
            sep="\t",
            quote = F
            );


res_wt_vs_parp7mut_liver_df$Entrez <- mapIds(org.Mm.eg.db, 
                                             res_wt_vs_parp7mut_liver_df$Ensembl,
                                             keytype="ENSEMBL", 
                                             column="ENTREZID",
                                             multiVals = "first"
                                             );

res_wt_vs_parp7mut_liver_df$Symbol <- mapIds(org.Mm.eg.db, 
                                             res_wt_vs_parp7mut_liver_df$Entrez,
                                             keytype="ENTREZID", 
                                             column="SYMBOL",
                                             multiVals = "first"
                                             );

res_wt_vs_parp7mut_liver_df$Gene_Name <- mapIds(org.Mm.eg.db, 
                                                res_wt_vs_parp7mut_liver_df$Entrez,
                                                keytype="ENTREZID", 
                                                column="GENENAME",
                                                multiVals = "first"
                                                );

res_wt_vs_parp7mut_liver_df$Symbol <- as.character(res_wt_vs_parp7mut_liver_df$Symbol);

res_wt_vs_parp7mut_liver_df$Gene_Name <- as.character(res_wt_vs_parp7mut_liver_df$Gene_Name);

write.table(res_wt_vs_parp7mut_liver_df,
            file="res_wt_vs_parp7mut_liver_df.txt",
            col.names = T,
            row.names = F,
            sep="\t",
            quote = F
            );

res_wt_vs_parp7mut_highfat_liver_df$Entrez <- mapIds(org.Mm.eg.db, 
                                                     res_wt_vs_parp7mut_highfat_liver_df$Ensembl,
                                                     keytype="ENSEMBL", 
                                                     column="ENTREZID",
                                                     multiVals = "First"
                                                     );

res_wt_vs_parp7mut_highfat_liver_df$Symbol <- mapIds(org.Mm.eg.db, 
                                                     res_wt_vs_parp7mut_highfat_liver_df$Entrez,
                                                     keytype="ENTREZID", 
                                                     column="SYMBOL",
                                                     multiVals = "first"
                                                     );

res_wt_vs_parp7mut_highfat_liver_df$Gene_Name <- mapIds(org.Mm.eg.db, 
                                                        res_wt_vs_parp7mut_highfat_liver_df$Entrez,
                                                        keytype="ENTREZID", 
                                                        column="GENENAME",
                                                        multiVals = "first"
                                                        );

res_wt_vs_parp7mut_highfat_liver_df$Symbol <- as.character(res_wt_vs_parp7mut_highfat_liver_df$Symbol);

res_wt_vs_parp7mut_highfat_liver_df$Gene_Name <- as.character(res_wt_vs_parp7mut_highfat_liver_df$Gene_Name);

write.table(res_wt_vs_parp7mut_highfat_liver_df,
            file="res_wt_vs_parp7mut_highfat_liver_df.txt",
            col.names = T,
            row.names = F,
            sep="\t",
            quote = F
            );

res_wt_vs_parp7mut_BAT_df$Entrez <- mapIds(org.Mm.eg.db, 
                                           res_wt_vs_parp7mut_BAT_df$Ensembl,
                                           keytype="ENSEMBL", 
                                           column="ENTREZID",
                                           multiVals = "first"
                                           );

res_wt_vs_parp7mut_BAT_df$Symbol <- mapIds(org.Mm.eg.db, 
                                           res_wt_vs_parp7mut_BAT_df$Entrez,
                                           keytype="ENTREZID", 
                                           column="SYMBOL",
                                           multiVals = "first"
                                           );

res_wt_vs_parp7mut_BAT_df$Gene_Name <- mapIds(org.Mm.eg.db, 
                                              res_wt_vs_parp7mut_BAT_df$Entrez,
                                              keytype="ENTREZID", 
                                              column="GENENAME",
                                              multiVals = "first"
                                              );


res_wt_vs_parp7mut_BAT_df$Symbol <- as.character(res_wt_vs_parp7mut_BAT_df$Symbol);

res_wt_vs_parp7mut_BAT_df$Gene_Name <- as.character(res_wt_vs_parp7mut_BAT_df$Gene_Name);


write.table(res_wt_vs_parp7mut_BAT_df,file="res_wt_vs_parp7mut_BAT_df.txt",
            col.names = T,
            row.names = F,
            sep="\t",
            quote = F
            );


res_wt_vs_parp7mut_highfat_BAT_df$Entrez <- mapIds(org.Mm.eg.db, 
                                                   res_wt_vs_parp7mut_highfat_BAT_df$Ensembl,
                                                   keytype="ENSEMBL", 
                                                   column="ENTREZID",
                                                   multiVals = "first"
                                                   );

res_wt_vs_parp7mut_highfat_BAT_df$Symbol <- mapIds(org.Mm.eg.db, 
                                                   res_wt_vs_parp7mut_highfat_BAT_df$Entrez,
                                                   keytype="ENTREZID", 
                                                   column="SYMBOL",
                                                   multiVals = "first"
                                                   );

res_wt_vs_parp7mut_highfat_BAT_df$Gene_Name <- mapIds(org.Mm.eg.db, 
                                                      res_wt_vs_parp7mut_highfat_BAT_df$Entrez,
                                                      keytype="ENTREZID", 
                                                      column="GENENAME",
                                                      multiVals = "first"
                                                      );


res_wt_vs_parp7mut_highfat_BAT_df$Symbol <- as.character(res_wt_vs_parp7mut_highfat_BAT_df$Symbol);
res_wt_vs_parp7mut_highfat_BAT_df$Gene_Name <- as.character(res_wt_vs_parp7mut_highfat_BAT_df$Gene_Name);


write.table(res_wt_vs_parp7mut_highfat_BAT_df,
            file = "res_wt_vs_parp7mut_highfat_BAT_df.txt",
            col.names = T,
            row.names = F,
            sep = "\t",
            quote = F
            );


res_wt_vs_parp7mut_WAT_df$Entrez <- mapIds(org.Mm.eg.db, 
                                           res_wt_vs_parp7mut_WAT_df$Ensembl,
                                           keytype = "ENSEMBL", 
                                           column = "ENTREZID",
                                           multiVals = "first"
                                           );

res_wt_vs_parp7mut_WAT_df$Symbol <- mapIds(org.Mm.eg.db, 
                                           res_wt_vs_parp7mut_WAT_df$Entrez,
                                           keytype = "ENTREZID", 
                                           column = "SYMBOL",
                                           multiVals = "first"
                                           );

res_wt_vs_parp7mut_WAT_df$Gene_Name <- mapIds(org.Mm.eg.db, 
                                              res_wt_vs_parp7mut_WAT_df$Entrez,
                                              keytype = "ENTREZID", 
                                              column = "GENENAME",
                                              multiVals = "first"
                                              );


res_wt_vs_parp7mut_WAT_df$Symbol <- as.character(res_wt_vs_parp7mut_WAT_df$Symbol);

res_wt_vs_parp7mut_WAT_df$Gene_Name <- as.character(res_wt_vs_parp7mut_WAT_df$Gene_Name);

write.table(res_wt_vs_parp7mut_WAT_df,
            file = "res_wt_vs_parp7mut_WAT_df.txt",
            col.names = T,
            row.names = F,
            sep = "\t",
            quote = F
            );

res_wt_vs_parp7mut_highfat_WAT_df$Entrez <- mapIds(org.Mm.eg.db, 
                                                   res_wt_vs_parp7mut_highfat_WAT_df$Ensembl,
                                                   keytype = "ENSEMBL", 
                                                   column = "ENTREZID",
                                                   multiVals = "first"
                                                   );

res_wt_vs_parp7mut_highfat_WAT_df$Symbol <- mapIds(org.Mm.eg.db, 
                                                   res_wt_vs_parp7mut_highfat_WAT_df$Entrez,
                                                   keytype = "ENTREZID", 
                                                   column = "SYMBOL",
                                                   multiVals = "first"
                                                   );

res_wt_vs_parp7mut_highfat_WAT_df$Gene_Name <- mapIds(org.Mm.eg.db, 
                                                      res_wt_vs_parp7mut_highfat_WAT_df$Entrez,
                                                      keytype = "ENTREZID", 
                                                      column = "GENENAME",
                                                      multiVals = "first"
                                                      );


res_wt_vs_parp7mut_highfat_WAT_df$Symbol <- as.character(res_wt_vs_parp7mut_highfat_WAT_df$Symbol);

res_wt_vs_parp7mut_highfat_WAT_df$Gene_Name <- as.character(res_wt_vs_parp7mut_highfat_WAT_df$Gene_Name);

write.table(res_wt_vs_parp7mut_highfat_WAT_df,
            file = "res_wt_vs_parp7mut_highfat_WAT_df.txt",
            col.names = T,
            row.names = F,
            sep = "\t",
            quote = F
            );

norm_count <- cpm(dds,
               normalised = TRUE
               );

write.table(rownames(dds),
            file = "Background_genes.txt",
            col.names = F,
            row.names = F,
            sep = "\t",
            quote = F
            );

write.table(norm_count,
            file = "CPM_normalised_counts.txt",
            col.names = T,
            row.names = T,
            sep = "\t",
            quote = F
            );



