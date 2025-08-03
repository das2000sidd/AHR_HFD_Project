setwd("~/Desktop/PhD_Project_related/AHR_HFD_Study");

library(VennDiagram);

liver <- read.table(file = "res_wt_vs_parp7mut_liver_df.txt",
					header = T,
					sep = "\t",
					stringsAsFactors = F,
					quote = ""
					); 
					
hf_liver <- read.table(file = "res_wt_vs_parp7mut_highfat_liver_df.txt",
					   header = T,
					   sep="\t",
					   stringsAsFactors = F,
					   quote = ""
					   ); 

bat <- read.table(file = "res_wt_vs_parp7mut_BAT_df.txt",
				  header = T,
				  sep = "\t",
				  stringsAsFactors = F,
				  quote = ""
				  ); 
				  
hf_bat <- read.table(file = "res_wt_vs_parp7mut_highfat_BAT_df.txt",
					 header = T,
					 sep = "\t",
					 stringsAsFactors = F,
					 quote = ""); 

wat <- read.table(file = "res_wt_vs_parp7mut_WAT_df.txt",
				  header = T,
				  sep = "\t",
				  stringsAsFactors = F,
				  quote = ""
				  ); 
				  
hf_wat <- read.table(file = "res_wt_vs_parp7mut_highfat_WAT_df.txt",
					 header = T,
					 sep = "\t",
					 stringsAsFactors = F,
					 quote = ""
					 ); 


liver_up <- subset(liver,
				   liver$log2FoldChange > 1 & liver$padj < 0.05
				   );
					
bat_up <- subset(bat,
				 bat$log2FoldChange > 1 & bat$padj < 0.05
				 );
				 
wat_up <- subset(wat,
				 wat$log2FoldChange > 1 & wat$padj < 0.05
				 );

liver_dn <- subset(liver,
				   liver$log2FoldChange < -1 & liver$padj < 0.05
				   );
				   
bat_dn <- subset(bat,
				bat$log2FoldChange < -1 & bat$padj < 0.05
				);
				
wat_dn <- subset(wat,
				wat$log2FoldChange < -1 & wat$padj < 0.05
				);


hf_liver_up <- subset(hf_liver,
				      hf_liver$log2FoldChange > 1 & hf_liver$padj < 0.05
				      );
				      
hf_bat_up <- subset(hf_bat,
					hf_bat$log2FoldChange > 1 & hf_bat$padj < 0.05
					);
					
hf_wat_up <- subset(hf_wat,
					hf_wat$log2FoldChange > 1 & hf_wat$padj < 0.05
					);

hf_liver_dn <- subset(hf_liver,
					  hf_liver$log2FoldChange < -1 & hf_liver$padj < 0.05
					  );
					  
hf_bat_dn <- subset(hf_bat,
					hf_bat$log2FoldChange < -1 & hf_bat$padj < 0.05
					);
					
hf_wat_dn <- subset(hf_wat,
					hf_wat$log2FoldChange < -1 & hf_wat$padj < 0.05
					);




pdf('Overlap_genes_up_BAT_tissue_2_diets.pdf');

p1 <- venn.diagram(list(
						`BAT, Control` = bat_up$Ensembl,
						`BAT, High fat` = hf_bat_up$Ensembl
						),
				   filename=NULL,
				   cex=2,
				   fill=c("yellow", "green"),
				   main.cex=c(1),
				   cat.cex=c(1.5, 1.5),
				   cat.dist=c(.025, .025),
				   main="Overlap of genes up in BAT tissue for HF and control diet",
				   cat.pos=-c(30, -30)
				   );

grid.draw(p1);

dev.off();


pdf('Overlap_genes_dn_BAT_tissue_2_diets.pdf');

p2 <- venn.diagram(list(
						`BAT, Control` = bat_dn$Ensembl,
						`BAT, High fat` = hf_bat_dn$Ensembl
						),
				   filename=NULL,
				   cex=2,
				   fill=c("yellow", "green"),
				   main.cex=c(1),
				   cat.cex=c(1.5, 1.5),
				   cat.dist=c(.025, .025),
				   main="Overlap of genes up in BAT tissue for HF and control diet",
				   cat.pos=-c(30, -30)
				   );

grid.draw(p2);

dev.off();


pdf('Overlap_genes_up_WAT_tissue_2_diets.pdf');

p3 <- venn.diagram(list(
						`WAT, Control`= wat_up$Ensembl,
						`WAT, High fat`= hf_wat_up$Ensembl
						),
				   filename=NULL,
				   cex=2,
				   fill=c("yellow", "green"),
				   main.cex=c(1),
				   cat.cex=c(1.5, 1.5),
				   cat.dist=c(.025, .025),
				   main="Overlap of genes up in WAT tissue for HF and control diet",
				   cat.pos=-c(30, -30)
				   );

grid.draw(p3);

dev.off();


pdf('Overlap_genes_dn_WAT_tissue_2_diets.pdf');

p4 <- venn.diagram(list(
						`WAT, Control`= wat_dn$Ensembl,
						`WAT, High fat`= hf_wat_dn$Ensembl
						),
				   filename=NULL,
				   cex=2,
				   fill=c("yellow", "green"),
				   main.cex=c(1),
				   cat.cex=c(1.5, 1.5),
				   cat.dist=c(.025, .025),
				   main="Overlap of genes down in WAT tissue for HF and control diet",
				   cat.pos=-c(30, -30)
				   );

grid.draw(p4);

dev.off();


pdf('Overlap_genes_up_Liver_tissue_2_diets.pdf');

p5 <- venn.diagram(list(
						`Liver, Control`= liver_up$Ensembl,
						`Liver, High fat`= hf_liver_up$Ensembl
						),
				   filename=NULL,
				   cex=2,
				   fill=c("yellow", "green"),
				   main.cex=c(1),
				   cat.cex=c(1.5, 1.5),
				   cat.dist=c(.025, .025),
				   main="Overlap of genes up in liver tissue for HF and control diet",
				   cat.pos=-c(30, -30)
				   );

grid.draw(p5);

dev.off();


pdf('Overlap_genes_dn_Liver_tissue_2_diets.pdf');

p6 <- venn.diagram(list(
						`Liver, Control` = liver_dn$Ensembl,
						`Liver, High fat` = hf_liver_dn$Ensembl
						),
				   filename=NULL,
				   cex=2,
				   fill=c("yellow", "green"),
				   main.cex=c(1),
				   cat.cex=c(1.5, 1.5),
				   cat.dist=c(.025, .025),
				   main="Overlap of genes down in liver tissue for HF and control diet",
				   cat.pos=-c(30, -30)
				   );

grid.draw(p6);

dev.off();





