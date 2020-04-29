library(data.table)

hmp<-fread("test.hmp.txt")
pos<-fread("genestr.csv")
atg.pos<-as.integer(pos[nrow(pos),4])
date<-data.table()
path<-dir("result/")
for (i in path){
  print(i)
  data<-fread("mlm+q/mlm_q_res70.txt")
  data<-data%>%dplyr::filter(Marker!="None")
  data$start<-data$Pos
  data$end<-data$Pos
  pos<-data.table(pos)
  pos<-na.omit(pos)
  data<-data.table(data)
setkey(pos,start,end)
  aa<-foverlaps(data,pos,type="any", which=TRUE)
  data$structure1<-pos[aa$yid,"structure1"]
  data$structure2<-pos[aa$yid,"structure2"]
colnames(hmp)[1]<-"Marker"
  dat<-merge.data.frame(data,hmp,by.x = "Marker")
  date<-rbind(date,dat)
}

fwrite(dat,"test.result.txt",row.names = F)

