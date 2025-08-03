setwd("~/Desktop/PhD_Project_related/AHR_HFD_Study")


library(dplyr);
library(clusterProfiler)
library(msigdbr)
library(org.Mm.eg.db)
library(magrittr)


`%ni%` = Negate(`%in%`);


liver <- read.table(file="res_wt_vs_parp7mut_liver_df.txt",
		header = T,
		sep = "\t",
		stringsAsFactors = F,
		quote = ""
		); ## 8214
					
hf_liver <- read.table(file="res_wt_vs_parp7mut_highfat_liver_df.txt",
		header = T,
		sep = "\t",
		stringsAsFactors = F,
		quote = ""
		); ## 9021

bat <- read.table(file="res_wt_vs_parp7mut_BAT_df.txt",
		header = T,
		sep = "\t",
		stringsAsFactors = F,
		quote = ""); ## 9582
				
hf_bat <- read.table(file="res_wt_vs_parp7mut_highfat_BAT_df.txt",
		header = T,
		sep = "\t",
		stringsAsFactors = F,
		quote = ""
		); ##9286

wat <- read.table(file="res_wt_vs_parp7mut_WAT_df.txt",
		header = T,
		sep = "\t",
		stringsAsFactors = F,
		quote = ""
		);## 9311
				  
hf_wat <- read.table(file="res_wt_vs_parp7mut_highfat_WAT_df.txt",
		header = T,
		sep = "\t",
		stringsAsFactors = F,
		quote = ""
		); ## 9582


liver_up <- subset(hf_liver,hf_liver$log2FoldChange > 1 & hf_liver$padj < 0.05);

bat_up <- subset(hf_bat,hf_bat$log2FoldChange > 1 & hf_bat$padj < 0.05);

wat_up <- subset(hf_wat,hf_wat$log2FoldChange > 1 & hf_wat$padj < 0.05);


liver_up_g <- liver_up$Ensembl;

bat_up_g <- bat_up$Ensembl;

wat_up_g <- wat_up$Ensembl;


length(liver_up_g);

length(bat_up_g);

length(wat_up_g);



liver_dn <- subset(hf_liver,hf_liver$log2FoldChange < -1 & hf_liver$padj < 0.05);

bat_dn <- subset(hf_bat,hf_bat$log2FoldChange < -1 & hf_bat$padj < 0.05);

wat_dn <- subset(hf_wat,hf_wat$log2FoldChange < -1 & hf_wat$padj < 0.05);

liver_dn_g <- liver_dn$Ensembl;

bat_dn_g <- bat_dn$Ensembl;

wat_dn_g <- wat_dn$Ensembl;


length(liver_dn_g);

length(bat_dn_g);

length(wat_dn_g);

## Enrichment analysis
mm_msigdb_df <- msigdbr(species = "Mus musculus");

head(mm_msigdb_df)

#Filter the human data frame to the KEGG pathways that are included in the
# curated gene sets
hs_GO_df <- mm_msigdb_df %>%
  dplyr::filter(
  	gs_cat == "C5", # This is to filter only to the C2 curated gene sets
  	gs_subcat %in% c("GO:BP","GO:CC","GO:MF") # This is because we only want KEGG pathways
  	);

hs_KEGG_df <- mm_msigdb_df %>%
  dplyr::filter(
  	gs_cat == "C2", # This is to filter only to the C2 curated gene sets
  	gs_subcat %in% c("CP:KEGG") # This is because we only want KEGG pathways
  	);


## ONLY RESULTS TO USE
all_background_genes <- c(liver$Ensembl,bat$Ensembl,wat$Ensembl,hf_liver$Ensembl,hf_bat$Ensembl,hf_wat$Ensembl);

all_background_genes < unique(all_background_genes);



## AHRKO
GO_ora_results_wat_hf_up_g <- enricher(
	gene = wat_up_g, # A vector of your genes of interest
	pvalueCutoff = 0.05, # Can choose a FDR cutoff
	pAdjustMethod = "BH",
	universe = all_background_genes,# Method to be used for multiple testing correction
	TERM2GENE = dplyr::select(hs_GO_df,
				gs_name,
				ensembl_gene
				)
	);

View(GO_ora_results_wat_hf_up_g@result);

enrich_plot <- enrichplot::dotplot(GO_ora_results_wat_hf_up_g, showCategory = 15,font.size = 10,title = "GO term enrichment for genes up in WAT tissue PARP7H532A mutant relative to WT high fat diet",orderBy = "p.adjust", decreasing = FALSE);

enrich_plot;


write.table(GO_ora_results_ahrko@result,file = "GO_Enrichment_terms_for_genes_up_in_M2_relative_to_M0_in_AHRKO.txt",col.names = T,row.names = T,sep = "\t",quote = F);




