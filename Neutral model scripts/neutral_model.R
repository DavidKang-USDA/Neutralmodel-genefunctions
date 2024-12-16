library(vegan)
library(plyr)
source("Burns.et.al.2016_neutral_model_code.R")


#count<-read.table("Raw data/feature-table-Wolb-Spiro-filt100-neg.txt",header=T)
#tax<-read.table("Raw data/taxonomy_all_unass.txt", header=T)
#meta.all=read.table("Raw data/trained_metadata.txt", header=T)

meta.all=read.table("trained_metadata.txt", header=T)
count<-read.table("table230-Wolb-Spiro-100-neg_NEW.txt", header=T, row.names=1)

tax<-read.table("taxonomy230-silva.txt", header=T)

count<-read.table("community.table.txt",header=T)
tax<-read.table("ASV.table.txt", header=T)

row.names(tax)<-tax$OTU_ID

tcount<-t(count)

sumcount<-apply(tcount, 1, sum)
tcount<-tcount[-c(101, 105),]


### Sp S
#countSpS<-count_neof
#rar<-as.data.frame(rrarefy(countSpS, 10000))

#tcompl.neo.mat<-t(compl.neo.mat)
#rar<-as.data.frame(rrarefy(tcompl.neo.mat, 10000))

#### Rarefy

rar<-as.data.frame(rrarefy(tcount, 10000))
sumcount<-apply(rar, 1, sum)
trar<-as.data.frame(t(rar))


### Count matrix
complRar<-merge(tax,trar,by=0)
complRar.mat<-complRar[, -c(1:9)]
rownames(complRar.mat)<-complRar[,1]  ##row=OTU, col=ind
tcomplRar.mat<-t(complRar.mat)        ##row=ind, col=OTU
#colnames(tcomplRar.mat)<-compl.Rar.tax$OTU_ID

compl.Rar.tax<-complRar[, c(1:9)]
row.names(compl.Rar.tax)<-compl.Rar.tax$OTU_ID
#  levels(compl.Rar.tax$Order)<-c(levels(compl.Rar.tax$Order), "unassigned")
#  compl.Rar.tax$Order[is.na(compl.Rar.tax$Order)]<-"unassigned"
#    levels(compl.Rar.tax$Class)<-c(levels(compl.Rar.tax$Class), "unassigned")
#    compl.Rar.tax$Class[is.na(compl.Rar.tax$Class)]<-"unassigned"
  
  

meta176<-meta.all[-c(101, 105),]


NeutMod<-sncm.fit(tcomplRar.mat, stats=FALSE, taxon=compl.Rar.tax)
NeutModSt<-sncm.fit(tcomplRar.mat, taxon=compl.Rar.tax)

  logp<-log(NeutMod$p)
  colrO<-rainbow(63)
  colrP<-rainbow(17)
  colrC<-rainbow(39)
  colrG<-rainbow(186)
 
 windows(6, 6) 
 par(mgp=c(2.2, 1,0), cex=0.8)
 plot(logp, NeutMod$freq, cex=0.7, pch=19,
      #col =NeutMod$colorPCoABray, 
      col=NeutMod$color,
      xlim=c(-13, -2),
      xlab="log(mean relative abundance)",
      ylab="occurence frequency") 
  points(logp, NeutMod$freq.pred, cex=0.1, col="red", type="l")
  points(logp, NeutMod$pred.upr, cex=0.1, col="red", type="l", lty=2)
  points(logp, NeutMod$pred.lwr, cex=0.1, col="red", type="l", lty=2)
  R2<-expression(R^2~"= 0.54")
  m<-expression(m~"= 2.57e"^-3)
  text(-3, 0.1, R2, cex=1)
  text(-2.8, 0.05, m, cex=1)
  
  legend(-13.4, 1.03, legend= levels_Neutmod_21, pch=19, col=levels_colors, bty="n", cex=0.8)
  #legend(-13.4, 1.03, legend= levels(NeutMod$Order)[1:31], pch=19, col=colrO[1:31], bty="n", cex=0.6)
  #legend(-11.2, 1.03, legend= levels(NeutMod$Order)[32:61], pch=19, col=colrO[32:61], bty="n", cex=0.6)
  #legend(-12, 1, legend= levels(NeutMod$phylum), pch=19, col=colrP, bty="n", cex=0.8)
  #legend(-12, 1, legend= levels(NeutMod$class)[1:20], pch=19, col=colrC[1:20], bty="n", cex=0.6)
  #legend(-10, 1, legend= levels(NeutMod$class)[21:39], pch=19, col=colrC[21:39], bty="n", cex=0.6)

  
levels_Neutmod_21<-c("Actinomycetales","Bacillales","Bacteroidales",
                     "BD7-3", "Burkholderiales","Caulobacterales","Cellvibrionales",
                     "Cytophagales","Enterobacteriales","EW055", "Flavobacteriales",
                     "Lactobacillales", "Neisseriales","Orbales", "Pseudomonadales",
                     "Rhizobiales","Rhodobacterales", "Sphingobacteriales","Sphingomonadales", 
                     "Xanthomonadales", "other")
  
levels_colors<-c("deepskyblue1", "darkorchid3", "darkgrey",
                 "mediumpurple1", "chocolate","burlywood3", "darkorange1",
                 "darkolivegreen2","darkseagreen3","tan","gold", 
                 "darkgreen", "brown3", "darkmagenta", "cornflowerblue", 
                 "darkolivegreen3", "dodgerblue3","deeppink", "darkslategray3",
                 "indianred2", "black")

  
  
##color only for above and below
  NeutMod$colorGraph<-NeutMod$Order 
  levels(NeutMod$colorGraph) <- c(levels(NeutMod$colorGraph),"zzz")
  NeutMod$colorGraph[NeutMod$freq>NeutMod$pred.lwr & NeutMod$freq<NeutMod$pred.upr]<-"zzz"

  
## identify points under or above CI
NeutMod$colorBR<-rep("black", length(NeutMod$p))
NeutMod$colorBR<-ifelse(NeutMod$freq<NeutMod$pred.lwr, "blue",
              ifelse(NeutMod$freq>NeutMod$pred.upr, "red", "black"))
NeutMod$color<-ifelse(NeutMod$Order == "Actinomycetales", "deepskyblue1",
                      ifelse(NeutMod$Order == "Cellvibrionales", "darkorange1",
                             ifelse(NeutMod$Order == "Bacillales", "darkorchid3",
                                    ifelse(NeutMod$Order == "Bacteroidales", "darkgrey",
                                           ifelse(NeutMod$Order == "BD7-3", "mediumpurple1",
                                                  ifelse(NeutMod$Order == "Burkholderiales", "chocolate",
                                                         ifelse(NeutMod$Order == "Caulobacterales", "burlywood3",
                                                                ifelse(NeutMod$Order == "Cytophagales", "darkolivegreen2",
                                                                       ifelse(NeutMod$Order == "Enterobacteriales", "darkseagreen3",
                                                                              ifelse(NeutMod$Order == "EW055", "tan",
                                                                                     ifelse(NeutMod$Order == "Flavobacteriales", "gold",
                                                                                            ifelse(NeutMod$Order == "Lactobacillales", "darkgreen",
                                                                                                   ifelse(NeutMod$Order == "Neisseriales", "brown3",
                                                                                                          ifelse(NeutMod$Order == "Orbales", "darkmagenta",
                                                                                                                 ifelse(NeutMod$Order == "Pseudomonadales", "cornflowerblue",
                                                                                                                        ifelse(NeutMod$Order == "Rhizobiales", "darkolivegreen3",
                                                                                                                               ifelse(NeutMod$Order == "Rhodobacterales", "dodgerblue3",
                                                                                                                                      ifelse(NeutMod$Order == "Sphingobacteriales", "deeppink",
                                                                                                                                             ifelse(NeutMod$Order=="Sphingomonadales", "darkslategray3",
                                                                                                                                                   ifelse(NeutMod$Order == "Xanthomonadales", "indianred2", "black"))))))))))))))))))))


### color only points from PCoA 15 ASVs
NeutMod$colorPCoABray<-NeutMod$Order
NeutMod$colorPCoABray<-ifelse(NeutMod$OTU_ID %in% cap.rel.vecs.bray$OTU_ID, NeutMod$color, "black")

###
NeutMod$colorPCoAJac<-NeutMod$order
NeutMod$colorPCoAJac<-ifelse(NeutMod$OTU_ID %in% cap.rel.vecs.jac.tax$OTU_ID, NeutMod$color, "black")


##### identify above/below neutral model
NeutModabove<-NeutMod[NeutMod$colorBR == "red", c(9:15)]  ##566   => ~51%
NeutModbelow<-NeutMod[NeutMod$colorBR == "blue", c(9:15)]    ##69  => ~6%
NeutModbetween<-NeutMod[NeutMod$colorBR == "black", c(9:15)]  ##472  => ~43%

NeutModaboveOrd<-as.data.frame(NeutModabove$Order)
colnames(NeutModaboveOrd)<-"Order"
NeutModaboveOrd$number_above<-rep(1, length(NeutModaboveOrd$Order))
NeutModaboveOrd.sum<-ddply(NeutModaboveOrd, "Order", numcolwise(sum))
NeutModaboveOrd.sum$average_above<-NeutModaboveOrd.sum$number_above/length(NeutModaboveOrd$number_above)*100

NeutModbelowOrd<-as.data.frame(NeutModbelow$Order)
colnames(NeutModbelowOrd)<-"Order"
NeutModbelowOrd$number_below<-rep(1, length(NeutModbelowOrd$Order))
NeutModbelowOrd.sum<-ddply(NeutModbelowOrd, "Order", numcolwise(sum))
NeutModbelowOrd.sum$average_below<-NeutModbelowOrd.sum$number_below/length(NeutModbelowOrd$number_below)*100

NeutModbetweenOrd<-as.data.frame(NeutModbetween$Order)
colnames(NeutModbetweenOrd)<-"Order"
NeutModbetweenOrd$number_between<-rep(1, length(NeutModbetweenOrd$Order))
NeutModbetweenOrd.sum<-ddply(NeutModbetweenOrd, "Order", numcolwise(sum))
NeutModbetweenOrd.sum$average_between<-NeutModbetweenOrd.sum$number_between/length(NeutModbetweenOrd$number_between)*100

order_neutMod<-merge(NeutModaboveOrd.sum, NeutModbetweenOrd.sum, all=TRUE)
order_neutMod<-merge(order_neutMod, NeutModbelowOrd.sum, all=TRUE)                   
