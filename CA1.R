library(GEOquery)
library(affy)
library(oligo)
library(limma)
library(tidyverse)



gse126202 <- getGEO('GSE126202')
gse126202 <- gse126202[[1]]
class(gse126202)

#gse126202_celdata <- read.celfiles(list.celfiles(getwd(),                               full.names=T,
#                                                 listGzipped = T))


varLabels(gse126202)

pd <- pData(gse126202)
pd['cel_file'] <- str_split(pd$supplementary_file,"/") %>% map_chr(tail,1)

gse126202_celdata <- read.celfiles(
  paste0(pd$cel_file),         #CEL files already in the work directory (not in the .tar file)#
  phenoData = phenoData(gse126202)
)

class(gse126202_celdata) #HTAFeatureSet

#see how the experiment was designed directly#
a <- pData(gse126202_celdata)[,c("genotype/variation:ch1",
                            "treatment:ch1")]
write.csv(a,
          "C:\\Learning\\BL5631\\Microarray\\pData.csv",
          row.names = F)

library(pd.mta.1.0)
gse126202_eset <- rma(gse126202_celdata)

#not included in the essay#
design <- model.matrix( ~ 0 + gse126202_eset[['treatment:ch1']])
colnames(design) <- levels(as.factor(gse126202_eset[['treatment:ch1']]))
fit <- lmFit(gse126202_eset,design)
fitted.ebayes <- eBayes(fit)
topTable(fitted.ebayes)
summary(decideTests(
  fitted.ebayes[,"Selumetinib, pulsatile"],lfc=1)) 
#strange results returned from the above code, not included in the essay#


###treat1:Selumetinib, continuous;
###treat2:Selumetinib, pulsatile;
###control1:Vehicle for continuous selumetinib group;
###control2:Vehicle for pulsatile selumetinib group


#colnames(design) <- c('Treat1','Treat2','Control1','Control2')
#contrast_matrix <- makeContrasts(Treat1-Control1,
#                                 Treat2-Control2,
#                                 Treat1-Treat2,
#                                 levels=design)
#contrast_matrix_new <- makeContrasts(
#  Treat1,
#  levels = design
#)
#contrast_matrix

#fit <- lmFit(gse126202_eset,design)
#fit2 <- contrasts.fit(fit,contrasts=contrast_matrix)
#fit2 <- eBayes(fit2)
#summary(decideTests(fit2,lfc=1))

library(dplyr) #include in tidyverse package
pd <- pData(gse126202_eset)
pd <- rename(
  pd,cell_type="tissue:ch1",
  treatment="treatment:ch1")
pd$treatment <- as.factor(pd$treatment)
levels(pd$treatment) <- c("Treat1",
                          "Treat2",
                          "Control1",
                          "Control2")
design_1 <- model.matrix(~ 0 + pd$treatment)
colnames(design_1) <- as.factor(levels(pd$treatment))
design_1

#column names too long to process by R
#contrast_matrix_new <- makeContrasts(
  #Selumetinib_continuous=Selumetinib..continuous-Vehicle.for.continuous.selumetinib.group,
  #Selumetinib_pulsatile=Selumetinib..pulsatile-Vehicle.for.pulsatile.selumetinib.group,
  #cell_with_drug=Selumetinib..continuous-Selumetinib..pulsatile,
  #interaction=(Vehicle.for.continuous.selumetinib.group-Vehicle.for.pulsatile.selumetinib.group)-(Selumetinib..continuous-Selumetinib..pulsatile),
  #levels = levels(pd$treatment)
#)
contrast_matrix_1 <- makeContrasts(
  Selumetinib_continuous=Treat1-Control1,
  Selumetinib_pulsatile=Treat2-Control2,
  cell_with_drug=Treat1-Treat2,
  interaction=(Control1-Control2)-(Treat1-Treat2),
  levels = design_1
)



install.packages("knitr")
c <- knitr::kable(contrast_matrix_1)

gse126202_fit <- lmFit(gse126202_eset,design_1)
gse126202_fit2 <- contrasts.fit(gse126202_fit,
                                contrasts=contrast_matrix_1)
gse126202_fit2 <- eBayes(gse126202_fit2)
summary(decideTests(gse126202_fit2,lfc=1))

volcanoplot(gse126202_fit2,
            coef=2,
            main="Volcanoplot of differently expressed genes"
)
interesting_genes <- topTable(gse126202_fit2,
                              number=Inf,
                              p.value = 0.05,
                              lfc=2)


heatmap(
  exprs(gse126202_eset)[rownames(interesting_genes),],
  margins = c(6,10))




BiocManager::install("mta10probeset.db")
BiocManager::install("mta10transcriptcluster.db")
BiocManager::install("org.Mm.eg.db",force = T)
library(org.Mm.eg.db)
library(mta10probeset.db)
library(mta10transcriptcluster.db)
ls('package:mta10probeset.db')
ls('package:mta10transcriptcluster.db')
ls('package:org.Mm.eg.db')
ps <- rownames(topTable(gse126202_fit2))

unlist(mget(ps,mta10probesetSYMBOL))
unlist(mget(ps,mta10transcriptclusterSYMBOL))
unlist(mget(ps,org.Mm.egSYMBOL))

columns(mta10transcriptcluster.db)
columns(org.Mm.eg.db)

keytypes(mta10transcriptcluster.db)
keytypes(org.Mm.eg.db)
head(keys(mta10transcriptcluster.db,
          keytype="PROBEID"))

#no probeID in org.Mm.eg.db package

AnnotationDbi::select(
  mta10transcriptcluster.db,
  ps,c("SYMBOL","ENTREZID","GENENAME"),
  keytype="PROBEID")
