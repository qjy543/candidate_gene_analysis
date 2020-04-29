rm(list = ls())
library(rvest)
#BiocManager::install("Biostrings")
library(Biostrings)
### 1. find genestructure online-------
web<-read_html("Phytozome v12.1_ Gene.html") # input web from phytozome
cc<-c()
for(i in 1:11){
  select<-paste0('#isc_76 > div > span:nth-child(',i,')')
  str<-web%>%html_nodes(select)%>%html_text()
  cc<-c(cc,str)
}

set<-DNAStringSet(cc)

names(set)<-c("5UTR",paste0("exon",3:length(cc)-2),"3UTR")

### 2. load gene infomation -------
inbred<-read.table("2_sampleID/Inbreds.txt",stringsAsFactors = F)[,1]
landrace<-read.table("2_sampleID/Landraces.txt",stringsAsFactors = F)[,1]
teosinter<-read.table("2_sampleID/Teosinters.txt",stringsAsFactors = F)[,1]
all<-read.table("2_sampleID/sampleID.txt",stringsAsFactors = F)[,1]

a<-dir()
fa<-a[grep("fa",a)]
fa
for(i in fa){
  names<-unlist(strsplit(i,"[.]"))[1]
  print(names)
  seq<-readDNAMultipleAlignment(i)
  seqq<-seq@unmasked
  seq_inbred<-seqq[which(names(seqq)%in%inbred),]
  seq_landrace<-seqq[which(names(seqq)%in%landrace),]
  seq_teosinter<-seqq[which(names(seqq)%in%teosinter),]
  seq_all<-seqq[which(names(seqq)%in%all),]
  writeXStringSet(seq_inbred,paste0("candidate_gene/",names,"_A.fa"))
  writeXStringSet(seq_landrace,paste0("candidate_gene/",names,"_L.fa"))
  writeXStringSet(seq_teosinter,paste0("candidate_gene/",names,"_T.fa"))
  writeXStringSet(seq_all,paste0("candidate_gene/",names,"_ALL.fa"))
}


### 3. match genestructre with gene information-------

str<-c("5UTR",paste0("exon",3:length(cc)-2),"3UTR")
## 3.1 utr & exon-------
pos<-data.frame()
pre.end<-0
for(i in 1:length(cc)){
  print(i)
  x<-set[[i]]
seq.pos<-seqq[1:3,]
l<-length(x)
start<-subseq(x,start=1,end=5)
end<-subseq(x,start=l-6,end=l)
print(start)
start.df<-data.frame(vmatchPattern(start, seq.pos)$`1.B73GENE`)
if(nrow(start.df)==1){
  start.pos<-as.integer(start.df[1,1])
} else{
  start.c<-as.integer(start.df[,1])
  start.pos<-start.c[start.c>pre.end][1]
}

end.df<-data.frame(vmatchPattern(end,seq.pos)$`1.B73GENE`)
if(nrow(end.df)==1){
  end.pos<-as.integer(end.df[1,2])
  pre.end<-end.pos
} else{
  end.c<-as.integer(end.df[,2])
  end.pos<-max(end.c)
}
pre.end<-end.pos
per.pos<-data.frame(start=start.pos,end=end.pos)
pos<-rbind(pos,per.pos)
}
pos$structure2<-str
pos$structure1<-c("5UTR",rep("exon",length(cc)-2),"3UTR")
## 3.2 intron-------
intron.pos<-data.frame()
for ( i in 2:length(cc)){
a<-pos[i-1,2]
b<-pos[i,1]
c<-b-a
if(c==1){
  print(paste0("no intron between ",pos[i-1,3]," and ",pos[i,3]))
  intron.per.pos<-data.frame()
  } else {
  intron.start.pos<-pos[i-1,2]+1
  intron.end.pos<-pos[i,1]-1
  intron.per.pos<-data.frame(start=intron.start.pos,
                             end=intron.end.pos,
                             structure1="intron",
                             structure2=paste0("intron",i-2))
}
intron.pos<-rbind(intron.pos,intron.per.pos)
}
intron.pos
## 3.3 up & down-------
up.down<-data.frame(start=c(1,max(pos[,2])+1),
                    end = c(min(pos[,1])-1,length(seqq[[1]])),
                    structure1=c("upstream","downstream"),
                    structure2=c("upstream","downstream"))

## 3.4 ATG-----------
atg<-data.frame(start=pos[2,1],
                end = NA,
                structure1="ATG",
                structure2="ATG")
## 3.4 combine and output-------
pos<-rbind(pos,intron.pos)
pos<-rbind(pos,up.down)
pos<-pos[order(pos[,1]),]
pos<-rbind(pos,atg)
pos$gene<-strsplit(fa,"[.]")[[1]][1]
pos<-pos%>%dplyr::select(gene,structure1,structure2,start,end)
write.csv(pos,"genestr.csv",row.names = F,quote = F)


