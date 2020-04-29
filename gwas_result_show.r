rm(list = ls())
setwd("/home_extend/qjy/Project/select_gene/")
if(F){
  BiocManager::install("treeio")
  BiocManager::install("Biostrings")
  BiocManager::install("ggtree")
  install.packages("ggmsa")
  install.packages("seqmagick")
  install.packages("cowplot")
}
library(ggplot2)
library(data.table)
library(dplyr)
library(ggsci)
library(pegas)
library(ggmsa)
library(seqmagick)
library(cowplot)
library(ggtree)
date<-fread("test.result.txt")
pos<-fread("genestr.csv")
a<-dir()
fa<-a[grep("fa",a)]
genename<-unlist(strsplit(fa,"[.]"))[1]
dat<-date
p.adjust<-0.01/mean(table(dat$Trait))
date.sig<-dat%>%dplyr::filter(p<p.adjust)
date.sig<-date.sig[order(date.sig$Pos),]
date.unsig<-dat%>%dplyr::filter(p>=p.adjust)
colnames(date.sig)[27]<-"shape"
colnames(date.unsig)[27]<-"shape"
levels(as.factor(date.sig$Trait))
levels(as.factor(date.sig$Marker))
levels(as.factor(date.sig$structure2))
fwrite(date.sig,"test.sig.txt",row.names = F)
date.sig[which(date.sig$Marker=="snp-1110"),]
res<-date.sig%>%dplyr::select(Trait,Marker,pos,p,MarkerR2,structure2,alleles)
res
fwrite(res,"sig.res.txt",row.names = F)
atg.pos<-as.integer(pos[nrow(pos),4])
res$trait<-matrix(unlist(strsplit(res$Trait,"[_]")),ncol = 3,byrow = T)[,3]

p1<-ggplot()+
    geom_hline(yintercept = -log10(p.adjust),linetype=5,col="red")+
    geom_point(data=date.unsig,aes(x=pos,y=-log10(p)),col="grey60",size=1)+
    geom_point(data=date.sig,aes(x=pos,y=-log10(p),col=Trait,shape=as.factor(shape)),size=1.5)+
    geom_segment(data= pos,aes(x=start,xend=end,y=0,yend=0)) +
    geom_rect(data=pos%>%dplyr::filter(structure1!="intron"),
            aes(xmin=start, xmax=end,ymin=-0.4,ymax=0.4,fill=structure1))+
    scale_fill_aaas()+
    scale_color_d3()+
    ggtitle(paste0(genename," gwas result"))+xlab("")+
    theme_classic()+
    theme(axis.text.x = element_blank())
p1
ggsave(paste0(genename,"_gwas_result.pdf"),width = 8,height = 6)

sigsnp<-date.sig[!duplicated(date.sig$Marker),]
sigsnp<-sigsnp[order(sigsnp$Pos),]
sigsnp
sigsnp.matrix<-t(sigsnp[1:nrow(sigsnp),33:ncol(sigsnp)])
colnames(sigsnp.matrix)<-sigsnp$Marker
sigsnp.matrix
x<-sigsnp.matrix[,1]
x
removeN<-function(x){
max<-names(table(x))[order(as.numeric(table(x)),decreasing = T)[1]]
min<-names(table(x))[order(as.numeric(table(x)),decreasing = T)[2]]
x[x==max]<-2
x[x==min]<-0
x[x!=0&x!=2]<-1
x<-as.numeric(x)
return(x)
}
library(LDheatmap)
require(snpStats)
test<-(apply(sigsnp.matrix, 2,removeN))
colnames(test)<-sigsnp$Marker
row.names(test)<-row.names(sigsnp.matrix)
rgb.palette <- colorRampPalette(rev(c("blue", "orange", "red")), space = "rgb")
input<-test[1:348,]
gdat<-as(input,"SnpMatrix")
LDheatmap(gdat,genetic.distances=sigsnp$Pos,flip=TRUE,color=rgb.palette(18),SNP.name = sigsnp$Marker)
input<-test[349:416,]
gdat<-as(input,"SnpMatrix")
LDheatmap(gdat,genetic.distances=sigsnp$Pos,flip=TRUE,color=rgb.palette(18))
input<-test[417:448,]
gdat<-as(input,"SnpMatrix")
LDheatmap(gdat,genetic.distances=sigsnp$Pos,flip=TRUE,color=rgb.palette(18))


phe<-fread("all_phe_nomissing.txt")
data.sig
colnames(phe)[1]<-"Taxa"
id<-phe$Taxa
colnames(data.sig)
for (i in 1:nrow(data.sig)){
  print(i)
marker<-data.sig[i,"Marker"]
trait<-data.sig[i,"Trait"]
r2<-data.sig[i,"MarkerR2"]*100
alleles<-data.sig[i,"alleles"]
input1<-t(data.sig[i,which(colnames(data.sig)%in%id)])
input1<-data.frame(gen=input1,Taxa=row.names(input1))
colnames(input1)[1]<-"gen"
alleles1<-names(table(input1$gen))[order(table(input1$gen),decreasing = T)[1]]
alleles2<-names(table(input1$gen))[order(table(input1$gen),decreasing = T)[2]]
input2<-phe%>%dplyr::select(Taxa,trait)
input3<-merge.data.frame(input1,input2,by="Taxa")
colnames(input3)[2]<-"gen"
colnames(input3)[3]<-"value"
input4<-input3%>%dplyr::filter(gen==alleles1 | gen==alleles2)
input4$gen<-as.factor(as.character(input4$gen))
p<-t.test(value~gen,input4,alternative =  "greater")$p.value
print(p)
p<-ggplot(input4)+
  geom_boxplot(aes(x=gen,y=value,fill=gen))+
  ggtitle(paste0(trait," ",marker," ",r2,"%"," ",alleles))
ggsave(paste0(trait,"_",marker,".png"))
}
p<-ggplot(input4)+
  geom_boxplot(aes(x=gen,y=value,fill=gen))+
  ggtitle(paste0(trait," ",marker," ",r2,"%"," ",alleles))
ggsave(paste0(trait,"_",marker,".png"))