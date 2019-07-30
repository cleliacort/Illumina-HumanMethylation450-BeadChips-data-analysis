##################################################
###########		FIRST STEP		##################
################################################## 
#OUR AIM IS TO LOAD RAW DATA WITH MINFI AND CREATE AN OBJECT CALLE RGSET STORING THE DATA

#clean the workspace and upload the needed packages 
rm(list=ls())
source("https://bioconductor.org/biocLite.R")
biocLite(c("minfi","wateRmelon","shinyMethyl")) 

#we set our workdirectory 
setwd("C:/Users/Celia/Desktop/DRD_2018_LAB/REPORT")

#we upload the the package called minfi which allow us to load the idat files 
#and the experiment samplesheet
library(minfi)

#now we upload the whole directory containing the input file [VERIFICA CHE TIPO
#SONO] e the summurize file calle "SampleSheet" 
#http://127.0.0.1:22629/library/minfi/doc/minfi.html
#su questo sito ci sono le informazioni 
baseDir <- ("C:/Users/Celia/Desktop/DRD_2018_LAB/Report_input")
#quindi nell'oggetto baseDir abbiamo salvato tutti i nomi dei file contenuti nella
#Input_data (verifica cio) che possiamo leggere facendo il seguente comando
list.files(baseDir)

#altro passo è quello di salvare su in oggetto le informazioni del file SampleSheet
#che poi utilizzeremo e lo facciamo cosi
targets <- read.metharray.sheet(baseDir)

#now we want to be able to read inside this information store in the target
RGset <- read.metharray.exp(targets = targets,extend=T)
save(RGset,file="RGset.RData")

##2
#now the aim is to create two different dataframe which respectvely the red and
#green channel
Red <- data.frame(getRed(RGset))
Green <- data.frame(getGreen(RGset))

##3
#retrieve intensity of probe with specific address
Red_address=Red["61760464",]
Red_address
#         X9376538140_R01C02 X9376538140_R03C01 X9376538140_R03C02
#61760464                333                354                302
#         X9376538140_R06C01 X9403904089_R03C01 X9403904089_R04C02
#61760464               5052               8898               3910
#         X9403904089_R05C01 X9403904089_R06C01
#61760464                399                374

Green_address=Green["61760464",]
Green_address
#         X9376538140_R01C02 X9376538140_R03C01 X9376538140_R03C02
#61760464                190                189                213
#         X9376538140_R06C01 X9403904089_R03C01 X9403904089_R04C02
#61760464               7833               5836               8252
#         X9403904089_R05C01 X9403904089_R06C01
#61760464                262                189

##4
#Extract methil and unmethi signal
MSet.raw <- preprocessRaw(RGset)
save(MSet.raw,file="MSet_raw.RData")

##5.1
#quality check
#QC plot 
#Estimate sample-specific quality control (QC) for methylation data.
qc <- getQC(MSet.raw)
plotQC(qc)
#grafico

##5.2
#check negative controll 
#tutt le probes 
#getProbeInfo(RGset, type = "Control")
#check the intensity of negative controls with minfi
controlStripPlot(RGset, controls="NEGATIVE")
#grafico 

#Summarizing methylation data from a RGChannelSet or GenomicRatioSet into 
#a shinyMethylSet needed to launch the interactive interface of shinyMethyl.
#check the intensity of negative controls with shinymethyl
library(shinyMethyl)
summary <- shinySummarize(RGset)
save(summary,file="summary_shinyMethyl.RData")
load('summary_shinyMethyl.RData')
runShinyMethyl(summary)

##5.3
#p-value calculation
detP <- detectionP(RGset)
save(detP,file="detP.RData")
failed <- detP>0.05
summary(failed)
# 9376538140_R01C02 9376538140_R03C01 9376538140_R03C02 9376538140_R06C01
# Mode :logical     Mode :logical     Mode :logical     Mode :logical    
# FALSE:485010      FALSE:485100      FALSE:485119      FALSE:485414     
# TRUE :502         TRUE :412         TRUE :393         TRUE :98         
# 9403904089_R03C01 9403904089_R04C02 9403904089_R05C01 9403904089_R06C01
# Mode :logical     Mode :logical     Mode :logical     Mode :logical    
# FALSE:485391      FALSE:485393      FALSE:485155      FALSE:485071     
# TRUE :121         TRUE :119         TRUE :357         TRUE :441

##6
#b and M value
beta <- getBeta(MSet.raw)
M <- getM(MSet.raw) 

##6.1
#density plot of mean of beta values 
mean_of_beta <- apply(beta,1,mean)
d_mean_of_beta <- density(mean_of_beta,na.rm=T)

jpeg("Density_mean_beta.jpeg", width = 15, height = 15, units = 'cm', res = 300)
plot(d_mean_of_beta,main="Density of Beta Values",col="blue")
dev.off()
#lines(density(beta[,1],na.rm=T),col="red")
#lines(density(beta[,2],na.rm=T),col="green")

#mean of M value graph 
mean_of_M <- apply(M,1,mean)
d_mean_of_M <- density(mean_of_M,na.rm=T)
jpeg("Density_mean_M.jpeg", width = 15, height = 15, units = 'cm', res = 300)
plot(d_mean_of_M,main="Density of M Values",col="blue")
dev.off()

##6.2
pal<- rainbow(8)
#create the boxplot 
jpeg("Boxplot_beta.jpeg", width = 15, height = 15, units = 'cm', res = 300)
par(mar = c(12, 5, 4, 2) + 0.1)
boxplot(beta,show.names=TRUE,main=" Boxplot of Beta values",col=pal,las=2)
dev.off()

jpeg("Boxplot_M.jpeg", width = 15, height = 15, units = 'cm', res = 300)
par(mar = c(12, 5, 4, 2) + 0.1)
boxplot(M,show.names=TRUE,main=" Boxplot of M values",col=pal,las=2)
dev.off()


####NB:ci sono degli INFINITI RICORDI O CERCA DI CAPIRE IL PERCHE!!!!!!!!!
#MSG DI ERRORE
#In bplt(at[i], wid = width[i], stats = z$stats[, i], out = z$out[z$group ==  :
#Dati anomali (Inf, -Inf) nel boxplot 1 non disegnati

##6.3
#smothscatter for the beta value
jpeg("Smooth_beta.jpeg", width = 15, height = 15, units = 'cm', res = 300)
SD_of_beta <- apply(beta,1,sd)
smoothScatter(mean_of_beta,SD_of_beta,xlab="Mean of Beta", ylab="SD of Beta", main="Smoothscatter plot of Beta")
dev.off()

###PROBLEMA
#lines(lowess(mean_of_beta,SD_of_beta), col="red")
#mean_of_beta<-mean_of_beta[!mean_of_beta=="-Inf"]
#SD_of_beta<-SD_of_beta[!is.na(SD_of_beta)]

#smoothscatter for the M values
jpeg("Smooth_M.jpeg", width = 15, height = 15, units = 'cm', res = 300)
SD_of_M <- apply(M,1,sd)
smoothScatter(mean_of_M,SD_of_M,xlab="Mean of M", ylab="SD of M", main="Smoothscatter plot of M")
dev.off()

##6.4
################################################################################
#PCA analysis on beta values
#t return the traspose of the original matrix
MSet.raw <- preprocessRaw(RGset)
beta_original <- getBeta(MSet.raw)

beta_na_omit <- na.omit(beta_original)
summary(beta_na_omit) # per controllare che non ci siano più NA. Se notate i valori del summary sono gli stessi di summary(beta)
pca <- prcomp(t(na.omit(beta_na_omit)),scale=T)
pca$x #it give us a sort of table in which there are the samples and the value for the different componenents 
#we chose to plot the first and the second
jpeg("PCA.jpeg", width = 15, height = 15, units = 'cm', res = 300)
plot(pca$x[,1], pca$x[,2],xlim=c(-1000,1000),ylim=c(-500,500),cex=2,main="PCA for detecting batch effect", xlab="PC1", ylab="PC2")
text(pca$x[,1], pca$x[,2],labels=rownames(pca$x),cex=0.5,pos=1)

points(pca$x[pheno$Group=="A",],col="red",cex=2)
points(pca$x[pheno$Group=="B",],col="blue",cex=2)
dev.off()
#################################################################################

##7
#part of normalized data 
#to can start the normalization process we need to use a sort of code give us
#by the builder of microarray 

#NORMALIZED DATA
#predata divided into I and II
load("C:\\Users\\Celia\\Desktop\\DRD_2018_LAB\\REPORT\\Illumina450Manifest_clean.RData")
dfI <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type=="I",]
dfI <- droplevels(dfI)
dfII <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type=="II",]
dfII <- droplevels(dfII)

preprocessNoob_results <- preprocessNoob(RGset)
beta_preprocessNoob <- getBeta(preprocessNoob_results)
beta_preprocessNoob_I <- beta_preprocessNoob[rownames(beta_preprocessNoob) %in% dfI$IlmnID,]
beta_preprocessNoob_II <- beta_preprocessNoob[rownames(beta_preprocessNoob) %in% dfII$IlmnID,]

mean_of_beta_preprocessNoob_I <- apply(beta_preprocessNoob_I,1,mean)
mean_of_beta_preprocessNoob_II <- apply(beta_preprocessNoob_II,1,mean)
d_mean_of_beta_preprocessNoob_I <- density(mean_of_beta_preprocessNoob_I,na.rm=T)
d_mean_of_beta_preprocessNoob_II <- density(mean_of_beta_preprocessNoob_II,na.rm=T)

sd_of_beta_preprocessNoob_I <- apply(beta_preprocessNoob_I,1,sd)
sd_of_beta_preprocessNoob_II <- apply(beta_preprocessNoob_II,1,sd)
d_sd_of_beta_preprocessNoob_I <- density(sd_of_beta_preprocessNoob_I,na.rm=T)
d_sd_of_beta_preprocessNoob_II <- density(sd_of_beta_preprocessNoob_II,na.rm=T)

#NON NORMALIZED DATA 
beta_I = beta[rownames(beta) %in% dfI$IlmnID,]
beta_II = beta[rownames(beta) %in% dfII$IlmnID,]

mean_beta_I = apply(beta_I,1,mean)
mean_beta_II = apply(beta_II,1,mean)
d_mean_beta_I = density(mean_beta_I, na.rm=T)
d_mean_beta_II = density(mean_beta_II, na.rm=T)
sd_beta_I = apply(beta_I, 1, sd)
sd_beta_II = apply(beta_II, 1, sd)
d_sd_beta_I = density(sd_beta_I, na.rm=T)
d_sd_beta_II = density(sd_beta_II, na.rm=T)

#PLOT TOTALE 
jpeg("PCA.jpeg", width = 15, height = 15, units = 'cm', res = 300)
par(mfrow=c(2,3))
plot(d_mean_beta_I, col="blue", main="Raw data, beta")
lines(d_mean_beta_II,col="red")
plot(d_sd_beta_I, col="blue", main="Raw data, sd")
lines(d_sd_beta_II, col="red")
boxplot(beta,xlab="Samples",show.names=FALSE, main="Raw data, boxplot")

plot(d_mean_of_beta_preprocessNoob_I, col="blue", main="preprocessNoob, beta")
lines(d_mean_of_beta_preprocessNoob_II,col="red")
plot(d_sd_of_beta_preprocessNoob_I, col="blue", main="preprocessNoob, sd")
lines(d_sd_of_beta_preprocessNoob_II, col="red")
boxplot(beta_preprocessNoob,xlab="Samples",show.names=FALSE, main="preprocessNoob, boxplot")
dev.off()

##8
pheno <- read.csv("C:/Users/Celia/Desktop/DRD_2018_LAB/Report_input/SampleSheet_report.csv", sep=",")

#ANOVA TEST
MYanovaFunction <- function(x) {
anova_test <- aov(x~ pheno$Group)
return(summary(anova_test)[[1]][[5]][1])
} 

pValuesAnova <- apply(beta_preprocessNoob,1, MYanovaFunction)

# We can create a data.frame with all the beta values and the pValue column
final_Anova <- data.frame(beta_preprocessNoob, pValuesAnova)

#now we order the data in pValuesAnova
final_Anova <- final_Anova[order(final_Anova$pValuesAnova),]

##9
#now we select only those h
final_Anova_0.05 <- final_Anova[final_Anova$pValuesAnova<=0.05,]
#final_Anova_0.05=57784
#final_Anova=485512
dim(final_Anova_0.05)
dim(final_Anova)

#do some multiple correction 
raw_pValues <- final_Anova[,9]
corrected_pValues_BH <- p.adjust(raw_pValues,"BH")
corrected_pValues_Bonf <- p.adjust(raw_pValues,"bonferroni")
final_Anova_corrected <- data.frame(final_Anova, corrected_pValues_BH, corrected_pValues_Bonf)

final_Anova_0.05_bonfi <- final_Anova_corrected[final_Anova_corrected$corrected_pValues_Bonf<=0.05,]
dim( final_Anova_0.05_bonfi )
#2 11

final_Anova_0.05_BH <- final_Anova_corrected[final_Anova_corrected$corrected_pValues_BH<=0.05,]
dim(final_Anova_0.05_BH)
#7 11 

##10
#MANHATTAN PLOT 
install.packages("gap")
library(gap)

final_Anova_corrected <- data.frame(rownames(final_Anova_corrected),final_Anova_corrected)
head(final_Anova_corrected)

colnames(final_Anova_corrected)[1] <- "IlmnID"

final_Anova_corrected_annotated <- merge(final_Anova_corrected, Illumina450Manifest_clean,by="IlmnID")

save(final_Anova_corrected_annotated,file="final_Anova_corrected_annotated.RData")

db <- data.frame(final_Anova_corrected_annotated$CHR, final_Anova_corrected_annotated$MAPINFO, final_Anova_corrected_annotated$pValuesAnova)
#we have problem because the levels are not in order but we need to have them in order
levels(db$final_Anova_corrected_annotated.CHR)

db$final_Anova_corrected_annotated.CHR <- factor(db$final_Anova_corrected_annotated.CHR,levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"))

levels(db$final_Anova_corrected_annotated.CHR)

jpeg("Manhattan.jpeg", width = 15, height = 15, units = 'cm', res = 300)
palette <- rainbow(24)
mhtplot(db,control=mht.control(colors=palette))
axis(2,cex=0.5)
abline(a=-log10(0.01),b=0)
dev.off()

#VOLCAN PLOT
#la differenza in beta group A e gruppo B sta nel fatto che selezione id che stanno nel gruppo A ma le sonde sono lo stesso numero 
beta <- final_Anova_corrected[,2:9]
beta_groupA <- beta[,pheno$Group=="A"]
mean_beta_groupA <- apply(beta_groupA,1,mean)
beta_groupB <- beta[,pheno$Group=="B"]
mean_beta_groupB <- apply(beta_groupB,1,mean)
delta <- mean_beta_groupB-mean_beta_groupA

toVolcPlot <- data.frame(delta, -log10(final_Anova_corrected$pValuesAnova))
jpeg("Volcano.jpeg", width = 15, height = 15, units = 'cm', res = 300)
plot(toVolcPlot[,1], toVolcPlot[,2],xlab="Delta(B-A)", ylab= "-log10(pValues)", main="Volcano plot",pch=16)#pch=to chose the points
abline(a=-log10(0.01),b=0,col="red")#to put the threshold 

# I want tho highlight the probes (that is, the points), that have a nominal pValue<0.01 and a delta > 0.1)
toHighlight <- toVolcPlot[abs(toVolcPlot[,1])>0.1 & toVolcPlot[,2]>(-log10(0.01)),]
#head(toHighlight)
points(toHighlight[,1], toHighlight[,2],pch=16,cex=0.7,col="orange")
dev.off()

#compare with the corrected pValues BH
toVolcPlot_corrected_BH <- data.frame(delta, -log10(final_Anova_corrected$corrected_pValues_BH))

 plot(toVolcPlot_corrected_BH[,1], toVolcPlot_corrected_BH[,2],pch=16, xlab="Delta of means (B-A)", ylab= "-log10(p-values)", main="Volcano plot, multiple test correction")

abline(a=-log10(0.01),b=0,col="red")

toHighlight_corrected_BH<- toVolcPlot_corrected_BH[abs(toVolcPlot_corrected_BH[,1])>0.1 & toVolcPlot_corrected_BH[,2]>(-log10(0.01)),]
points(toHighlight_corrected_BH[,1], toHighlight_corrected_BH[,2],pch=16,cex=0.7,col="green")

final_Anova_0.05_BH
           corrected_pValues_BH corrected_pValues_Bonf
cg21915962          0.004500757            0.004500757
cg04145681          0.009573323            0.019146646
cg18827503          0.022040639            0.066121918
cg24136690          0.030273672            0.121094687
cg07469647          0.043021380            0.257699658
cg00374839          0.043021380            0.287557467
cg04010471          0.043021380            0.301149659

#we have 7 significant values for BH but when we trasform them into log only 
#two remain significative with pvalue greater then 2 


##11
final_Anova_corrected_annotated_sign <- final_Anova_corrected_annotated[final_Anova_corrected_annotated$pValuesAnova <0.05,]

final_Anova_corrected_annotated_sign <- droplevels(final_Anova_corrected_annotated_sign)

dist.euc=dist(t(final_Anova_corrected_annotated_sign[,2:9])) # euclidean distances between the rows
#distance between sample (traspose matrix)

mds=cmdscale(dist.euc)
pdf("8. MDS.pdf")

jpeg("MDS.jpeg", width = 15, height = 15, units = 'cm', res = 400)
plot(mds, main="MDS of nominal pValues", cex=2, xlab="x coordinates", ylab="y coordinates",xlim=c(-20,18),ylim=c(-4.5,8))
points(mds[pheno$Group=="A",],col="red",cex=2)
points(mds[pheno$Group=="B",],col="blue",cex=2)
legend("topleft",c("A","B"),pch=1,col=c("red","blue"))
coord1<-mds[,1]
coord2<-mds[,2]
text(coord1, coord2, labels=rownames(mds),cex=0.5,pos=1)
dev.off()

