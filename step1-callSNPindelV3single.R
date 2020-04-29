rm(list=ls())
getwd()
hmp.file="test2.hmp.txt"
#for (f in hmp.file){print(f)}
for (f in hmp.file){
  #rm(list=ls())
  gene=substring(f,1,regexpr("\\.", f)-1)
  df=read.delim(paste0(f),header=T, stringsAsFactors =F)
  for (i in 1:nrow(df)){
    if (any(df[i,12:ncol(df)]=="-"))
    {df[i,1]=paste("indel",df$pos[i],sep="_");df[i,6]="indel"}
    else{df[i,1]=paste("snp",df$pos[i],sep="_");df[i,6]="snp"}}
  df.new=data.frame()
  i=1
  while(i<nrow(df)){
    cat("i=",i,"\n")
    if (df[i,6]=="snp") {df.new=rbind(df.new,df[i,]);print("1. this is snp");i=i+1} 
    else if (df[i+1,6]=="snp"){
      #df[i,12:ncol(df)][which(df[i,12:ncol(df)]!="-")]=1  #df[i,12:ncol(df)][which(df[i,12:ncol(df)]=="-")]=0
      #df[i,2]=paste(1,0,sep="/")
      df.new=rbind(df.new,df[i,]);print("2. this is single indel");i=i+1} 
    else  {
      print("3. this is muti-indel")
      starti=i
      indel.count=1
      while(df[i,6]=="indel"){
        cond=identical(which(df[i,12:ncol(df)]=="-"),which(df[i+1,12:ncol(df)]=="-"))
        if(cond){
          print("3.1 this is hebing indel")
          i=i+1
          indel.count=indel.count+1
        }
        else {
          #df[starti,12:ncol(df)][which(df[starti,12:ncol(df)]!="-")]=indel.count
          #df[starti,12:ncol(df)][which(df[starti,12:ncol(df)]=="-")]=0
          #df[starti,2]=paste(df[4:5,2],collapse = ";")
          #???????????????????????????
          alle=unlist(strsplit(df[starti:i,2],"/"))
          df[starti,2]= paste(paste(alle[seq(1, length(alle), 2)],collapse = ""),
          paste(alle[seq(2, length(alle), 2)],collapse = ""),sep = "/")
          
          
          df.new=rbind(df.new,df[starti,])
          print("3.2 3.2 3.2");i=i+1;break} 
      }
    }
  }
  cat("combine indel is finished","\n")
  gs=read.table("genestr.csv",sep=",",header=T)
  gs$indel=NA;gs$snp=NA
  for (l in 1:(nrow(gs)-1)){
    a=df.new$assembly.[(df.new$pos<=gs$end[l])&(df.new$pos>=gs$start[l])]
    if (all(is.na(a))){gs[l,6:7]=NA}
    else if(length(a)==0){gs[l,6:7]=0}
    else {gs[l,6:7]=as.numeric(table(a))}
    print(l)
  }
  variation.stat=data.frame(aggregate(gs[,6:7],by=list(gs$structure1),FUN=sum,na.rm=T))
  variation.stat$Group.1=as.character(variation.stat$Group.1)
  variation.stat[nrow(variation.stat)+1,]=c("sum",sum(variation.stat$indel,na.rm = T),sum(variation.stat$snp,na.rm = T))
  df.new$pos=df.new$pos+1
  df.new$center=df.new$pos
  df.new$pos=df.new$pos-gs$start[gs$structure1=="ATG"]
  df.new$rs.=paste0(df.new$assembly.,df.new$pos)
  gs$start[1:(nrow(gs)-1)]=gs$start[1:(nrow(gs)-1)]-gs$start[gs$structure1=="ATG"]
  gs$end[1:(nrow(gs)-1)]=gs$end[1:(nrow(gs)-1)]-gs$start[gs$structure1=="ATG"]
  #variation.stat=aggregate(rs.~assembly.,df.new,FUN=length)
  #names(variation.stat)=c("type","count")
  #df.new$assembly.=NA
  write.table(df.new,paste0(gene,"_snpindel.hmp.txt"),sep="\t",quote=F,row.names = F)
  write.table(variation.stat,paste0(gene,".variation.stat.csv"),sep=",",row.names = F)
  write.table(gs,paste0(gene,".variation.stat2.csv"),sep=",",row.names = F)
  cat(f,"is finished","\n")
}
df.new$pos<-df.new$pos+atg.pos
ddf<-df.new[-grep("././.",df.new$alleles),]

write.table(df.new,"test.hmp.txt",sep="\t",quote=F,row.names = F)
df<-fread("test.hmp.txt")
