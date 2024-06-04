args <- commandArgs(trailingOnly = TRUE)
# 
# args[1]= '/project/DPDS/Wang_lab/shared/BCR_antigen/code/Cmai/data/example/output/binding_results_with_rank.csv'
# args[2]='/project/DPDS/Wang_lab/shared/BCR_antigen/code/Cmai/data/example/BCR_metrics.csv'
# args[3]=0.03
# args[4]=10
# args[5]=15
# args[6]=20

cmai_score = read.csv(args[1])
bcr_tb = read.csv(args[2])
thre <- as.numeric(args[3])
cutoff_string = args[4]
cutoffs = as.numeric(strsplit(gsub("[()]","",cutoff_string),",")[[1]])
# cutoff1 <- as.numeric(args[4])
# cutoff2 <- as.numeric(args[5])
# cutoff3 <- as.numeric(args[6])

antigens = unique(cmai_score$Antigen)
tmp=merge(cmai_score,bcr_tb,by = 'BCR_id',all.x = T)

binding_score = vector(length = length(antigens))
names(binding_score)=antigens

for (antigen in antigens){
  # Bing: check the irAE section of Yi/Yuqiu's paper. Similar idea
  # scores=sapply(1:length(cmai),function(i) {
    # tmp=cmai[[i]]
    tmp=tmp[tmp$Antigen==antigen,]
    tmp=tmp[order(-tmp$bcr_metric),]
    if (dim(tmp)[1]<20) {binding_score[antigen]=NA}
    # maybe use 0.01, 0.02, 0.04, or 0.05 as rank% cutoff
    # maybe use (5,10,20) as the top BCR clone cutoff
    score=sapply(cutoffs,function(x){
      mean(tmp$Rank[1:x]<thre)
    })
    # a=mean(tmp$Rank[1:cutoff1]<thre)
    # b=mean(tmp$Rank[1:cutoff2]<thre) 
    # c=mean(tmp$Rank[1:cutoff3]<thre)
    binding_score[antigen] = mean(score,na.rm=T)/thre
}

write.csv(as.data.frame(binding_score),args[5])

