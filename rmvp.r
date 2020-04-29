rm(list = ls())
library(rMVP)
setwd("Project/select_gene/")
hmp<-"/home_extend/qjy/Project/select_gene/test.hmp.txt"
phe<-"/home_extend/qjy/Project/select_gene/all_phe.rmvp.txt"
kin<-"/home_extend/qjy/Project/select_gene/mvp.kinship.txt"
pc<-"/home_extend/qjy/Project/select_gene/mvp.Q.txt"
if(F){
  MVP.Data(fileHMP =hmp,
           filePhe=phe,
           fileKin=kin,
           filePC=pc,
           out="maize"
  )
}



genotype <- attach.big.matrix("maize.geno.desc")
phenotype <- read.table("maize.phe",head=TRUE)
map <- read.table("maize.geno.map" , head = TRUE)
Kinship<-attach.big.matrix("maize.kin.desc")
Covariates<-attach.big.matrix("maize.pc.desc")

for(i in 2:ncol(phenotype)){
  imMVP <- MVP(
    phe=phenotype[, c(1, i)],
    geno=genotype,
    map=map,
    #K=Kinship,
    #CV.GLM=Covariates,
    #CV.MLM=Covariates,
    #CV.FarmCPU=Covariates,
    nPC.GLM=5,
    nPC.MLM=3,
    priority="speed",
    ncpus=10,
    vc.method="EMMA",
    maxLoop=10,
    method.bin="FaST-LMM",#"FaST-LMM","EMMA", "static"
    permutation.threshold=TRUE,
    permutation.rep=100,
    threshold=0.05,
    method=c("MLM"),file.output = F
  )
  data<-data.frame(imMVP$map,imMVP$mlm.results)
  date<-na.omit(data)
  write.csv(date,paste0("MVP.",colnames(phenotype)[i],".MLM.csv"),row.names = F,quote = F)
}
