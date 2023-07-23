require(tidyverse)
require(ggpubr)
require(gridExtra)
require(cowplot)
library(readxl)
library(glue)


## Import Brain span data 
rows_metadata <- read_csv("rows_metadata.csv")
columns_metadata <- read_csv("columns_metadata.csv")
expression_matrix <- read_csv("expression_matrix.csv", col_names = FALSE)

# change var name gene_symbol to GENE
rows_metadata <- rows_metadata %>% select(., 1:3, GENE = gene_symbol, 5)

# merge rows_metadata with expression
dat <- cbind(rows_metadata, expression_matrix)

## Check overall sample size 
samplecount <- columns_metadata %>% group_by(donor_id) %>% count() %>% arrange(desc(n))

## Check sample size for each of the available tissue * time point 
structure <- columns_metadata %>% group_by(structure_name) %>% count()
timepoint <- columns_metadata %>% group_by(age) %>% count()

## Recode Age into Weeks
columns_metadata <- columns_metadata %>% mutate(Weeks = ifelse(age =="8 pcw",8, 
                                                            ifelse(age =="9 pcw",9, 
                                                            ifelse(age =="12 pcw",12,
                                                            ifelse(age =="13 pcw",13,
                                                            ifelse(age =="16 pcw",16,
                                                            ifelse(age =="17 pcw",17,
                                                            ifelse(age =="19 pcw",19,
                                                            ifelse(age =="21 pcw",21,
                                                            ifelse(age =="24 pcw",24,
                                                            ifelse(age =="25 pcw",25,
                                                            ifelse(age =="26 pcw",26,
                                                            ifelse(age =="35 pcw",35,
                                                            ifelse(age =="37 pcw",37,
                                                            ifelse(age =="4 mos",53,
                                                            ifelse(age =="10 mos",77,
                                                            ifelse(age =="1 yrs",89,
                                                            ifelse(age =="2 yrs",141,
                                                            ifelse(age =="3 yrs",193,
                                                            ifelse(age =="4 yrs",245,
                                                            ifelse(age =="8 yrs",453,
                                                            ifelse(age =="11 yrs",609,
                                                            ifelse(age =="13 yrs",713,
                                                            ifelse(age =="15 yrs",817,
                                                            ifelse(age =="18 yrs",973,
                                                            ifelse(age =="19 yrs",1025,
                                                            ifelse(age =="21 yrs",1129,
                                                            ifelse(age =="23 yrs",1233,
                                                            ifelse(age =="30 yrs",1597,
                                                            ifelse(age =="36 yrs",1909,
                                                            ifelse(age =="37 yrs",1961,
                                                            ifelse(age =="40 yrs",2117, "no"))))))))))))))))))))))))))))))))

columns_metadata <- columns_metadata %>% mutate(Sex = ifelse(gender == "M", 1, ifelse(gender == "F", 0, "NA")))

## Import candidate driver genes from Supplementary Table 16
drivergenes <- read_excel("Consolidated-genesets.xlsx",sheet = "discordant")

geenesets<-unique(drivergenes$Set)


# Loop over each geneset
for (geeneset in geenesets) {
  
  geeneset_df <- drivergenes[drivergenes$Set == geeneset,]
  geeneset_df <- as.data.frame(unique(unlist(strsplit(geeneset_df$CommonGenes, ","))))
  colnames(geeneset_df) <- c("GENE")
  geeneset_df$C1 <- 1
  
  outname <- glue("GenesetCluster_Discordant_{geeneset}")
  
  ## Merge driver genes to Brain Span genes
  dat2 <- right_join(geeneset_df, dat, by = "GENE")
  
  ## Remove row.number col X1 from dat2 dataset
  dat2 <- dat2 %>% select(-X1)
  
  ## Extract and transpose gene expression columns and merge them
  cog1 <- dat2 %>% group_by(C1) %>% summarise(across(X2:X525, ~ mean(.x, na.rm = TRUE))) %>% select(-1) %>% t() %>% as.data.frame() %>% setNames("C1")
  total <- dat2 %>% summarise(across(X2:X525, ~ mean(.x, na.rm = TRUE))) %>% t() %>% as.data.frame() %>% setNames("TOTALexp")
  
  vlong <- cbind(cog1, total)
  
  dat3 <- cbind(columns_metadata, vlong)
  
  ## Generate Violin plots for visualization
  C1 <- dat3 %>% ggplot(aes(y = C1, x = factor(Weeks, level = c('8','9','12','13','16','17','19','21','24','25','26','35','37','53','77','89','141','193','245','453','609','713','817','973','1025','1129','1233','1597','1909','1961','2117')))) +
    geom_point() + geom_violin() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1)) +
    labs(title = "Common Loci") + theme(axis.text.y = element_text(size = 6)) +
    theme(axis.text.x = element_text(size = 6)) + theme(plot.title = element_text(size = 9)) +
    theme(axis.title.x = element_blank()) + theme(axis.title.y = element_blank())
  
  ## Save database
  write.table(dat3, file = glue("{outname}.txt"), sep = " ", row.names = FALSE, quote = FALSE)
  Common <- dat3 %>% select(donor_id, Weeks, Sex, structure_acronym, Expression = C1) %>% mutate(Trait = "Common")
  
  dat4 <- rbind(Common)
  dat4$Weeks <- as.numeric(dat4$Weeks)
  
  # Create the combined plot
  C1_plot <- C1 + labs(title = "Panel A") + theme(plot.title = element_text(hjust = 0.5))
  
  dat4_plot <- dat4 %>% ggplot(aes(x = Weeks, y = Expression, color = Trait)) +
    geom_point() +
    theme(legend.position = "none") +
    geom_smooth(method = 'lm', formula = y ~ x, linetype = "dashed", size = 0.7) +
    labs(title = "Panel B") + theme(plot.title = element_text(hjust = 0.5))
  
  combined_plot <- plot_grid(C1_plot, dat4_plot, nrow = 2, labels = c("A", "B"), label_size = 14, label_fontface = "bold")
  
  # Save the combined plot as a single TIFF file
  tiff(filename = glue("{outname}_Combined.tiff"), width = 7, height = 8.3, units = "in", pointsize = 10, bg = "white", res = 800)
  print(combined_plot)
  dev.off()
}
