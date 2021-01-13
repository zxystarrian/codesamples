set.seed(999)

.libPaths("/dartfs-hpc/rc/home/8/f002s78/R/x86_64-redhat-linux-gnu-library/3.5")
library("EMVC",lib.loc="/dartfs-hpc/rc/home/8/f002s78/R/x86_64-redhat-linux-gnu-library/3.6")
library("GSVA",lib.loc="/dartfs-hpc/rc/home/8/f002s78/R/x86_64-redhat-linux-gnu-library/3.5")
library("stringr",lib.loc="/dartfs-hpc/rc/home/8/f002s78/R/x86_64-redhat-linux-gnu-library/3.5")
library("grpreg",lib.loc="/dartfs-hpc/rc/home/8/f002s78/R/x86_64-redhat-linux-gnu-library/3.5")
library("cvTools",lib.loc="/dartfs-hpc/rc/home/8/f002s78/R/x86_64-redhat-linux-gnu-library/3.5")
library("survival",lib.loc="/dartfs-hpc/rc/home/8/f002s78/R/x86_64-redhat-linux-gnu-library/3.5")
library("survcomp",lib.loc="/dartfs-hpc/rc/home/8/f002s78/R/x86_64-redhat-linux-gnu-library/3.5")
library("glmnet",lib.loc="/dartfs-hpc/rc/home/8/f002s78/R/x86_64-redhat-linux-gnu-library/3.5")
library("GSVA",lib.loc="/dartfs-hpc/rc/home/8/f002s78/R/x86_64-redhat-linux-gnu-library/3.5")

################################################################################
# Parameters input from PBS job scripts
Args <- commandArgs(T)
cohorts<-Args[1]
threshold<-Args[2]
threshold<-as.numeric(threshold)
noise<-as.numeric(Args[3])
assoc<-as.numeric(Args[4])

################################################################################
# set the repetition times
r_num=20
# set the K-fold croos validation
K = 5

################################################################################
# Read the pathway data from MSigDB
p <- c("/dartfs-hpc/rc/lab/F/FrostH/members/xzheng/integrate/all/h.all.v7.2.symbols")
filename <- paste(p,"gmt",sep = ".")
data <- file(filename,open="r")
n <- 1
pathway <- list()

while(TRUE){
  line <- readLines(data,n=1)
  if(length(line) == 0){
    break
  }
  pathway[n] <- line
  n <- n+1
}
close(data)

all_pathway <- list()
names <- c()
for(i in 1:length(pathway)){
  one_pathway <- strsplit(as.character(pathway[i]),"\t")[[1]]
  names[i] <- one_pathway[1]
  one_pathway <- one_pathway[-1]
  one_pathway <- one_pathway[-1]
  all_pathway[[i]] <- one_pathway
}
names(all_pathway) <- names
Pathway<-all_pathway

genes<-as.vector(unlist(Pathway))
genes<-unique(genes)

# Process the gene expression data .
expr<-read.table(file = paste('/dartfs-hpc/rc/lab/F/FrostH/members/xzheng/integrate/all/',cohorts,'_HiSeqV2',sep = ''),header = TRUE,row.names = 1,sep = '\t')
row.names(expr)<-gsub("-","",row.names(expr))

need_genes<-intersect(row.names(expr),genes)
expr<-expr[need_genes,]
patients <- colnames(expr)
n_sample <-length(patients)

# Shrinkage the data to only genes in the pathway
clustergenes <- intersect(need_genes,Pathway[[assoc]])
print(names(all_pathway)[assoc])
print(length(clustergenes))

expr_h <- na.omit(expr[clustergenes,]) # deletion of missing
expr_h <- scale(expr_h) # standardize variables

# Take a random cluster of 20 genes with inter-gene correlation > 0.2
tmp_s<-1
while(tmp_s){
  rg<-sample(clustergenes,20)
  M<-abs(cor(t(expr_h[rg,])))
  mc=sprintf("%.2f", mean(abs(M)))
  if(mc>0.2){
    tmp_s<-0
  }
}
print(mc)

# Shuf the gene rows except the clustered genes
unshufgenes=rg
shufgenes=need_genes[which(!(need_genes %in% unshufgenes))]
shufdata=expr[shufgenes,]
unshufdata=expr[unshufgenes,]
for (i in 1:dim(shufdata)[1]){
  shufdata[i,]<-sample(shufdata[i,])
}
expr<-rbind(shufdata,unshufdata)
row.names(expr)<-c(shufgenes,unshufgenes)
colnames(expr)<-patients

print(length(unshufgenes))
print(length(shufgenes))

# take the normal samples
if(length(grep("01$",names(expr)))!=0){
  expr_n <- as.matrix(expr[,-(grep("01$",names(expr)))])
}
# remove the normal samples
if(length(grep("01$",names(expr)))!=0){
  expr_t <- as.matrix(expr[,grep("01$",names(expr))])
}

# Run EMVC with tumor samples
annotations<-matrix(nrow=length(Pathway),ncol=length(need_genes))
row.names(annotations)<-names
colnames(annotations)<-need_genes
for(i in 1:dim(annotations)[1]){
  annotations[i,]=as.numeric(need_genes %in% Pathway[[i]])
}
data=t(expr_t)
EMVC.results=EMVC(data, annotations, bootstrap.iter=50, k.range=5:15, clust.method="kmeans",kmeans.nstart=1, kmeans.iter.max=10, hclust.method="average", hclust.cor.method="spearman")

filtered.opt.annotations = filterAnnotations(EMVC.results,threshold)
sum(filtered.opt.annotations)

# Get the new pathway definition with EMVA results
newPath<-list()
for(i in 1:dim(annotations)[1]){
  newPath[[i]]<-need_genes[filtered.opt.annotations[i,]!=0]
}
names(newPath)<-row.names(annotations)

num_gene<-length(need_genes)
print(num_gene)
num_pathway<-length(newPath)

######## run GSVA with new pathways and shuffled data
e1 <- cbind.data.frame(expr)
gsva_es <- gsva(as.matrix(e1), newPath, mx.diff=TRUE)
gsva_es <- cbind.data.frame(gsva_es)
names_gsva_es<-row.names(gsva_es)

################################################################################
# Start of the replications
################################################################################
c_all <- c()
c_all2<-list()
cp_all <- list()
cp_all2<-list()
pv_all<-c()
pv_all2<-list()

r=1

while(r<=r_num){
  print(r)

###### Generate the survival data
  T<-matrix(nrow = n_sample,ncol = 3)
  colnames(T)<-c("OS","Censoring","status")
  row.names(T)<-patients

  pick<-unshufgenes
  pick<-unique(pick)
  pick<-intersect(pick,need_genes)
  tt<-t(scale(t(expr),center = TRUE,scale = TRUE))

  if(noise==0){
    for (i in row.names(T)){
      ra=mean(tt[pick,i])+10
      T[i,"OS"]<-round(exp(ra)/10)
    }
  }else{
      for (i in row.names(T)){
        ra=mean(tt[pick,i])+rnorm(1,0,noise)+10
        T[i,"OS"]<-round(exp(ra)/10)
      }
  }

  T[,"Censoring"]<-sample(T[,"OS"],n_sample)
  T[,"status"]<-as.numeric(T[,"OS"]<T[,"Censoring"])
  colnames(T)<-c('OS.time','Censor','OS')
  pheno <- T
  surtime <- as.data.frame(na.omit(pheno[,c('OS.time','OS')]))
  if(length(surtime[-which(surtime$OS.time<=0),1])){
    surtime<-surtime[-which(surtime$OS.time<=0),]
  }

  n_obs <- n_sample
  surtime_c <- surtime

################################################################################
# cross validation
################################################################################

  cvt <- cvFolds(n_obs, K = K, R = 1, type = "random")
  folds <- as.data.frame(cbind(c(rep(rep(1:K),5000))[1:n_obs],cvt$subsets))
  colnames(folds) <- c("fold","index")

  ci = c()
  cpi = list()
  pv = c()

  for (i in 1:cvt$K){
    print(i)

# Split the data
    traning <- folds[folds$fold!=i,][,"index"]
    traning <- traning[order(traning)]
    test <- folds[folds$fold==i,][,"index"]
    test <- test[order(test)]
    td <- as.matrix(t((gsva_es))[traning,])
    colnames(td)<-rownames(gsva_es)
    tt <- as.matrix(t(gsva_es)[test,])
    colnames(tt)<-rownames(gsva_es)
    temp <- surtime_c[traning,]
    colnames(temp) <- c("time","status")
# Fit the univariable Cox model
    left_genes <- names(Pathway)[assoc]
    f4 <- as.formula(paste("Surv(surtime_c$OS.time[traning], surtime_c$OS[traning]) ~ ", paste(left_genes, collapse= "+")))
    try(fit_final2 <- coxph(f4, data=as.data.frame(td)),silent=TRUE)
    coeffs2 <- coef(summary(fit_final2))
    pv[i]=coeffs2[,5]
 # Use the estimated parameters to predict on the test data
    temp <- surtime_c[test,]
    colnames(temp) <- c("time","status")
    test_r<-predict(fit_final2,as.data.frame(tt),type="risk")
    try(con <- concordance.index(x=test_r, surv.time=surtime_c$OS.time[test], surv.event=surtime_c$OS[test]),silent = TRUE)
    try(ci[i] <- con$c.index,silent = TRUE)
    try(cpi[[i]] <- left_genes,silent = TRUE)

  }

# Save the results of each replication
  ci=replace(ci,is.na(ci),0.5)
  try(c_all[r]<-mean(ci,na.rm=TRUE),silent = TRUE)
  try(pv_all[r]<-mean(pv,na.rm=TRUE),silent = TRUE)
  try(c_all2[[r]]<-ci,silent = TRUE)
  try(pv_all2[[r]]<-pv,silent = TRUE)
  try(cp_all[[r]]<-cpi,silent = TRUE)
  if(length(cpi)!=5){
    next
  }
  cp_all2<-c(cp_all2,cpi)

  r=r+1

}
########################################################################################################
# End of replications
########################################################################################################

# Take the mean and sd
try(print(mean(c_all,na.rm=TRUE)),silent = TRUE)
try(print(sd(c_all,na.rm=TRUE)),silent = TRUE)
try(print(mean(pv_all,na.rm=TRUE)),silent = TRUE)
# Take the frequency of used predictors
tp <- as.vector(unlist(cp_all2))
length(tp)/length(cp_all2)
tpt<-table(tp)
length(tpt)

########################################################################################################
# Save the results
w<-paste(cohorts,mean(c_all,na.rm=TRUE),sd(c_all,na.rm=TRUE),length(tp)/length(cp_all2),length(tpt),mean(pv_all,na.rm=TRUE),sd(pv_all,na.rm=TRUE),n_obs,sep="\t")
ff <- file(paste("P_expr_combined_results","_",threshold,".txt",sep=""),open="at")
write.table(w,file=ff,append=T,sep="\n",quote=F,row.names=F,col.names=F)
close(ff)

block=expr_t[newPath[[assoc]],]
save(newPath,unshufdata,block,c_all2,pv_all2,cp_all,cp_all2,w,surtime,gsva_es,c_all,tp,tpt,file=paste(cohorts,"_P_expr","_",threshold,".RData",sep=""))
########################################################################################################
warnings()

