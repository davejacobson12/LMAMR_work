args<-commandArgs(TRUE)
data<-read.table(args[1],sep="\t",row.names=1,header=TRUE)

taxa_pos<-grep("Bacteria",names(data))

cat_id<-grep (args[2],names(data))


var_frame<-as.character(levels(as.factor(data[,cat_id])))

k<-length(var_frame)

results_data<-as.data.frame(matrix(NA,length(taxa_pos),length(var_frame)+2))

for(i in 1:length(taxa_pos)){
results_data[i,1:k]<-boxplot(data[,taxa_pos[i]]~data[,cat_id],plot=FALSE)$stats[                                                                                                             3,]
results_data[i,k+1]<-kruskal.test(data[,taxa_pos[i]]~data[,cat_id])$p.value
results_data[i,k+2]<-names(data)[taxa_pos[i]]
}

names(results_data)<-c(var_frame,"pval","Taxon")

results_data$bonf<-p.adjust(results_data$pval,method="bonf")
results_data$fdr<-p.adjust(results_data$pval,method="fdr")

write.table(results_data,file=args[3],sep="\t",quote=FALSE,row.names=FALSE)
