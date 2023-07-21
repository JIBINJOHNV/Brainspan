require(tidyverse)
require(ggpubr)
require(gridExtra)
require(cowplot)
library(readxl)

## Import Brain span data 
## / START

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
drivergenes <- read_excel("Supplementary_Table16_1.xlsx",sheet = "CandidateGenes")

## Merge driver genes to Brain Span genes
dat2 <- right_join(drivergenes, dat, by = "GENE")

## Remove row.number col X1 from dat2 dataset
dat2 <- dat2 %>% select(., -X1)

## Extract and transpose gene expression columns and merge them 
cog1 <- dat2 %>% group_by(C1) %>% summarise(across(X2:X525, ~ mean(.x, na.rm = TRUE))) %>% select(., -1) %>% t() %>% as.data.frame() %>% select(., C1 = V1)
cog3 <- dat2 %>% group_by(C3) %>% summarise(across(X2:X525, ~ mean(.x, na.rm = TRUE))) %>% select(., -1) %>% t() %>% as.data.frame() %>% select(., C3 = V1)
cog4 <- dat2 %>% group_by(C4) %>% summarise(across(X2:X525, ~ mean(.x, na.rm = TRUE))) %>% select(., -1) %>% t() %>% as.data.frame() %>% select(., C4 = V1)
total <- dat2 %>% summarise(across(X2:X525, ~ mean(.x, na.rm = TRUE))) %>% t() %>% as.data.frame() %>% select(., TOTALexp = V1)

vlong <- cbind(cog1, cog3, cog4, total)

dat3 <- cbind(columns_metadata, vlong)

## Generate Violin plots for visualization 
C1 <- dat3 %>% ggplot(., aes(y=C1, x=factor(Weeks, level = c('8','9','12','13','16','17','19','21','24','25','26','35','37','53','77','89','141','193','245','453','609','713','817','973','1025','1129','1233','1597','1909','1961','2117')))) + 
        geom_point() + geom_violin() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + 
        labs(title="Common Loci") + theme(axis.text.y = element_text(size=6)) + 
        theme(axis.text.x = element_text(size=6)) + theme(plot.title = element_text(size=9)) + 
        theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank()) 

C3 <- dat3 %>% ggplot(., aes(y=C3, x=factor(Weeks, level = c('8','9','12','13','16','17','19','21','24','25','26','35','37','53','77','89','141','193','245','453','609','713','817','973','1025','1129','1233','1597','1909','1961','2117')))) + 
        geom_point() + geom_violin() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + 
        labs(title="Discordent Loci") + theme(axis.text.y = element_text(size=6)) + 
        theme(axis.text.x = element_text(size=6)) + theme(plot.title = element_text(size=9)) + 
        theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank())

C4 <- dat3 %>% ggplot(., aes(y=C4, x=factor(Weeks, level = c('8','9','12','13','16','17','19','21','24','25','26','35','37','53','77','89','141','193','245','453','609','713','817','973','1025','1129','1233','1597','1909','1961','2117')))) + 
        geom_point() + geom_violin() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + 
        labs(title="Concordant Loci") + theme(axis.text.y = element_text(size=6)) + 
        theme(axis.text.x = element_text(size=6)) + theme(plot.title = element_text(size=9)) + 
        theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank())

## Save database
write.table(dat3, file = "BrainSpan_CogNCog_MetaLoci_data.txt", sep = " ", row.names=FALSE, quote = FALSE)
Common <- dat3 %>% select(., donor_id, Weeks, Sex, structure_acronym, Expression = C1) %>% mutate(Trait = "Common")
Discordent <- dat3 %>% select(., donor_id, Weeks, Sex, structure_acronym, Expression = C3) %>% mutate(Trait = "Discordent")
Concordant <- dat3 %>% select(., donor_id, Weeks, Sex, structure_acronym, Expression = C4) %>% mutate(Trait = "Concordant")

dat4 <- rbind(Common, Discordent,Concordant)
dat4$Weeks <- as.numeric(dat4$Weeks)

### Extract plots 
tiff(filename = "BrainSpan_PanelBi.tiff", width = 4.25, height = 5, units = "in", pointsize = 10, bg = "white",  res = 800)
    panelb <- ggarrange(C1,C3,C4, ncol=1, nrow=3)
    annotate_figure(panelb, left = text_grob("BrainSpan Normalized Gene Expression", color = "Black", face = "bold", size = 11, rot=90), bottom = text_grob("Weeks", color = "Black", face = "bold", size = 11))
dev.off()

tiff(filename = "BrainSpan_PanelA.tiff", width = 7, height = 3.3, units = "in", pointsize = 10, bg = "white",  res = 800)
    Cog_NonCog <- dat4 %>% ggplot(., aes(x = Weeks, y = Expression, color=Trait)) + geom_point() + facet_grid(rows=vars(Trait)) + theme(legend.position="none") + geom_smooth(method='lm', formula= y~x, linetype="dashed", size = 0.7)
    Cog_NonCog
dev.off()



library("lme4")
library("languageR")
library(lmerTest)

BS0 <- lmer(Expression ~ Trait + (1 | donor_id), data = dat4)
BS1 <- lmer(Expression ~ Trait*Weeks + (1 | donor_id), data = dat4)

BS0 <- lmer(Expression ~ Weeks + Sex + (1 | donor_id), data = dat4)
BS1 <- lmer(Expression ~ Sex + Weeks*Trait + (1 | donor_id), data = dat4)
anova(BS0, BS1, "chisq")

    # mixed model on individual meta-locus
dat3$Weeks <- as.numeric(dat3$Weeks)
GCA1 <- lmer(C1 ~ Weeks + Sex + (1 | donor_id), data = dat3)
GCA3 <- lmer(C3 ~ Weeks + Sex + (1 | donor_id), data = dat3)
GCA4 <- lmer(C4 ~ Weeks + Sex + (1 | donor_id), data = dat3)
### / END
anova(GCA1)


Common <- lmer(C1 ~ Weeks*Sex + (1 | donor_id), data = dat3)
Discordant <- lmer(C3 ~ Weeks*Sex + (1 | donor_id), data = dat3)
Concordant <- lmer(C4 ~ Weeks*Sex + (1 | donor_id), data = dat3)
### / END
anova(GCA1)
summ(GCA1)

anova(Common, Discordant,Concordant, "chisq")




Common <- dat3 %>% select(., donor_id, Weeks, Sex, structure_acronym, Expression = C1) %>% mutate(Trait = "Common")
Discordant <- dat3 %>% select(., donor_id, Weeks, Sex, structure_acronym, Expression = C3) %>% mutate(Trait = "Discordant")
Concordant <- dat3 %>% select(., donor_id, Weeks, Sex, structure_acronym, Expression = C4) %>% mutate(Trait = "Concordant")


ConcordantDiscordant <- dat3 %>% select(., donor_id, Weeks, Sex, structure_acronym, Expression = C1) %>% mutate(Trait = "ConcordantDiscordant")
Discordant <- dat3 %>% select(., donor_id, Weeks, Sex, structure_acronym, Expression = C3) %>% mutate(Trait = "Discordant")
Concordant1 <- dat3 %>% select(., donor_id, Weeks, Sex, structure_acronym, Expression = C4) %>% mutate(Trait = "Concordant1")



dat4 <- rbind(ConcordantDiscordant, Discordant,Concordant1)

dat4 <- rbind(Common, Discordant)

dat4$Weeks <- as.numeric(dat4$Weeks)

library("lme4")
BS0 <- lmer(Expression ~ Trait + (1 | donor_id), data = dat4)
BS1 <- lmer(Expression ~ Trait*Weeks + (1 | donor_id), data = dat4)

BS2 <- lmer(Expression ~ Weeks + Sex + (1 | donor_id), data = dat4)
BS3 <- lmer(Expression ~ Sex + Weeks*Trait + (1 | donor_id), data = dat4)

anova(BS0, BS1, "chisq")

anova(BS0)
summ(BS0)



### Draw Figure 5
### / START
tiff(filename = "BrainSpan_PanelA.tiff", width = 8, height = 5, units = "in", pointsize = 10, bg = "white",  res = 800)
    ggdraw() + draw_plot(Cog_NonCog, x = 0, y = 0, height = 1 ) + draw_plot(C1, x=0.5, y=0.8, height =0.15, width=0.45) + draw_plot(C1, x=0.5, y = 0.08, height = 0.15, width = 0.45) + draw_plot_label(label = "Model Pval= Trait * Weeks = 2.28e-133", x=0, y=0.52, size = 12)
dev.off()

tiff(filename = "BrainSpan_PanelA.tiff", width = 8, height = 5, units = "in", pointsize = 10, bg = "white",  res = 800)
    ggdraw() + draw_plot(Cog_NonCog, x = 0, y = 0, height = 1 ) + draw_plot_label(label = "Model Pval= Trait * Weeks = 2.28e-133")
dev.off()

