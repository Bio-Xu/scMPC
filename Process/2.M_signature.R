geneSet<-readRDS("GeneSets.rds",refhook = NULL)
source("Functions/2.OptimizeGeneSet.R")

#####GSE62321
load("/WorkSpace/scMetastasis/Data/bulk_data/GSE62321_exp.Rdata")
load("/WorkSpace/scMetastasis/Data/bulk_data/pData62321.Rdata")
system.time({
  cl<- makeCluster(2)      
  registerDoParallel(cl) 
  Result62321 <- foreach(
    i=1:25, 
    .errorhandling = "pass") %dopar% {
      FF<-OptimizeGeneSet(geneSet=unlist(geneSet[i]), erMatrix=GSE62321_exp, pData=pData62321)
      cat(FF)
    }
  stopCluster(cl)
}) 

#####GSE205209
load("/WorkSpace/scMetastasis/Data/bulk_data/GSE205209_exp.Rdata")
load("/WorkSpace/scMetastasis/Data/bulk_data/pData205209.Rdata")
system.time({
  cl<- makeCluster(4)      
  registerDoParallel(cl) 
  Result205209 <- foreach(
    i=1:25,
    .errorhandling = "pass") %dopar% {
      FF<-OptimizeGeneSet(geneSet=unlist(geneSet[i]), erMatrix=GSE205209_exp, pData=pData205209)
      cat(FF)
    }
  stopCluster(cl)
})

#####GSE125989
load("/WorkSpace/scMetastasis/Data/bulk_data/GSE125989_exp.Rdata")
load("/WorkSpace/scMetastasis/Data/bulk_data/pData125989.Rdata")
system.time({
  cl<- makeCluster(4)      
  registerDoParallel(cl)
  Result125989 <- foreach(
    i=1:25,
    .errorhandling = "pass") %dopar% {
      FF<-OptimizeGeneSet(geneSet=unlist(geneSet[i]), erMatrix=GSE125989_exp, pData=pData125989)
      cat(FF)
    }
  stopCluster(cl)
}) 

#####GSE73168
load("/WorkSpace/scMetastasis/Data/bulk_data/GSE73168_exp.Rdata")
load("/WorkSpace/scMetastasis/Data/bulk_data/pData73168.Rdata")
system.time({
  cl<- makeCluster(4)      
  registerDoParallel(cl) 
  Result73168 <- foreach(
    i=1:25, 
    .errorhandling = "pass") %dopar% {
      FF<-OptimizeGeneSet(geneSet=unlist(geneSet[i]), erMatrix=GSE73168_exp, pData=pData73168)
      cat(FF)
    }
  stopCluster(cl)
}) 

#####GSE180186
load("/WorkSpace/scMetastasis/Data/bulk_data/GSE180186_exp.Rdata")
load("/WorkSpace/scMetastasis/Data/bulk_data/pData180186.Rdata")

Result180186 <- list()
for(i in 1:25){
  fit<-try(OptimizeGeneSet(geneSet=unlist(geneSet[i]), erMatrix=GSE180186_exp, pData=pData180186),silent=FALSE)
  if("try-error" %in% class(fit))
  {next}
  else
  {FF<-OptimizeGeneSet(geneSet=unlist(geneSet[i]), erMatrix=GSE180186_exp, pData=pData180186)
  Result180186[[i]]<-FF}
}

#####GSE72718
load("/WorkSpace/scMetastasis/Data/bulk_data/GSE72718_exp.Rdata")
load("/WorkSpace/scMetastasis/Data/bulk_data/pData72718.Rdata")

Result72718 <- list()
for(i in 1:25){
  fit<-try(OptimizeGeneSet(geneSet=unlist(geneSet[i]), erMatrix=GSE72718_exp, pData=pData72718),silent=FALSE)
  if("try-error" %in% class(fit))
  {next}
  else
  {FF<-OptimizeGeneSet(geneSet=unlist(geneSet[i]), erMatrix=GSE72718_exp, pData=pData72718)
  Result72718[[i]]<-FF}
}

#####GSE224235
load("/WorkSpace/scMetastasis/Data/bulk_data/GSE224235_exp.Rdata")
load("/WorkSpace/scMetastasis/Data/bulk_data/pData224235.Rdata")

Result224235<-list()
for(i in 1:25){
  fit<-try(OptimizeGeneSet(geneSet=unlist(geneSet[i]), erMatrix=GSE224235_exp, pData=pData224235),silent=FALSE)
  if("try-error" %in% class(fit))
  {next}
  else
  {FF<-OptimizeGeneSet(geneSet=unlist(geneSet[i]), erMatrix=GSE224235_exp, pData=pData224235)
  Result224235[[i]]<-FF}
}

#####GSE60542
load("/WorkSpace/scMetastasis/Data/bulk_data/GSE60542_exp.Rdata")
load("/WorkSpace/scMetastasis/Data/bulk_data/pData60542.Rdata")
system.time({
  cl<- makeCluster(4)      
  registerDoParallel(cl)
  Result60542 <- foreach(
    i=1:25, 
    .errorhandling = "pass") %dopar% {
      FF<-OptimizeGeneSet(geneSet=unlist(geneSet[i]), erMatrix=GSE60542_exp, pData=pData60542)
    }
  stopCluster(cl)
}) 

#####GSE85258
load("/WorkSpace/scMetastasis/Data/bulk_data/GSE85258_exp.Rdata")
load("/WorkSpace/scMetastasis/Data/bulk_data/pData85258.Rdata")
system.time({
  cl<- makeCluster(4)      
  registerDoParallel(cl)
  Result85258 <- foreach(
    i=1:25,
    .errorhandling = "pass") %dopar% {
      FF<-OptimizeGeneSet(geneSet=unlist(geneSet[i]), erMatrix=GSE85258_exp, pData=pData85258)
      cat(FF)
    }
  stopCluster(cl)
}) 

####### Integrate
source("RankAggregate.RRA.R")
geneList1<-list(row.names(as.data.frame(Result224235[[1]])),
                 row.names(as.data.frame(Result125989[[1]])),
                 row.names(as.data.frame(Result85258[[1]])))

geneList2<-list(row.names(as.data.frame(Result62321[[2]])),row.names(as.data.frame(Result224235[[2]])),
                row.names(as.data.frame(Result125989[[2]])),
                row.names(as.data.frame(Result85258[[2]])))

geneList3<-list(row.names(as.data.frame(Result125989[[3]])),
                row.names(as.data.frame(Result85258[[3]])))

geneList4<-list(row.names(as.data.frame(Result62321[[4]])),row.names(as.data.frame(Result224235[[4]])),
                row.names(as.data.frame(Result125989[[4]])),
                row.names(as.data.frame(Result85258[[4]])))

geneList5<-list(row.names(as.data.frame(Result224235[[5]])),
                row.names(as.data.frame(Result205209[[5]])),
                row.names(as.data.frame(Result85258[[5]])))



geneList6<-list(row.names(as.data.frame(Result85258[[6]])))

geneList7<-list(row.names(as.data.frame(Result224235[[7]])),
                row.names(as.data.frame(Result60542[[7]])),
                row.names(as.data.frame(Result125989[[7]])),
                row.names(as.data.frame(Result85258[[7]])))

geneList8<-NULL

geneList9<-list(row.names(as.data.frame(Result125989[[9]])),row.names(as.data.frame(Result85258[[9]])))

geneList10<-NULL

geneList11<-NULL

geneList12<-list(row.names(as.data.frame(Result85258[[12]])))

geneList13<-list(row.names(as.data.frame(Result85258[[13]])))

geneList14<-list(row.names(as.data.frame(Result85258[[14]])))

geneList15<-NULL

geneList16<-NULL

geneList17<-list(row.names(as.data.frame(Result224235[[17]])),row.names(as.data.frame(Result85258[[17]])))

geneList18<-list(
  row.names(as.data.frame(Result224235[[18]])),
  row.names(as.data.frame(Result180186[[18]])),row.names(as.data.frame(Result125989[[18]])),
  row.names(as.data.frame(Result85258[[18]])))

geneList19<-list(row.names(as.data.frame(Result62321[[19]])),row.names(as.data.frame(Result224235[[19]])),
                 row.names(as.data.frame(Result125989[[19]])),
                 row.names(as.data.frame(Result85258[[19]])) )

geneList20<-list(row.names(as.data.frame(Result62321[[20]])),row.names(as.data.frame(Result224235[[20]])),
                 row.names(as.data.frame(Result180186[[20]])),
                 row.names(as.data.frame(Result85258[[20]])))


geneList21<-list(row.names(as.data.frame(Result62321[[21]])),row.names(as.data.frame(Result224235[[21]])),
                 row.names(as.data.frame(Result125989[[21]])),row.names(as.data.frame(Result85258[[21]])))

geneList22<-list(row.names(as.data.frame(Result125989[[22]])),row.names(as.data.frame(Result85258[[22]])))

geneList23<-NULL

geneList24<-list(row.names(as.data.frame(Result85258[[24]])))

geneList25<-list(row.names(as.data.frame(Result224235[[25]])), row.names(as.data.frame(Result85258[[25]])))

geneSet <- readRDS("GeneSets.rds")
geneSet<-GOterms
set1<-RankAggregate.RRA (geneList1,unlist(geneSet[1]))
set2<-RankAggregate.RRA (geneList2,unlist(geneSet[2]))
set3<-RankAggregate.RRA (geneList3,unlist(geneSet[3]))
set4<-RankAggregate.RRA (geneList4,unlist(geneSet[4]))
set5<-RankAggregate.RRA (geneList5,unlist(geneSet[5]))
set6<-RankAggregate.RRA (geneList6,unlist(geneSet[6]))
set7<-RankAggregate.RRA (geneList7,unlist(geneSet[7]))
set8<-NULL
set9<-RankAggregate.RRA (geneList9,unlist(geneSet[9]))
set10<-NULL
set11<-NULL
set12<-RankAggregate.RRA (geneList12,unlist(geneSet[12]))
set13<-RankAggregate.RRA (geneList13,unlist(geneSet[13]))
set14<-RankAggregate.RRA (geneList14,unlist(geneSet[14]))
set15<-NULL
set16<-NULL
set17<-RankAggregate.RRA (geneList17,unlist(geneSet[17]))
set18<-RankAggregate.RRA (geneList18,unlist(geneSet[18]))
set19<-RankAggregate.RRA (geneList19,unlist(geneSet[19]))
set20<-RankAggregate.RRA (geneList20,unlist(geneSet[20]))
set21<-RankAggregate.RRA (geneList21,unlist(geneSet[21]))
set22<-RankAggregate.RRA (geneList22,unlist(geneSet[22]))
set23<-NULL
set24<-RankAggregate.RRA (geneList24,unlist(geneSet[24]))
set25<-RankAggregate.RRA (geneList25,unlist(geneSet[25]))

source("3.2.EvaluateSignature.R")
load("D:/subject/exp/GSE85258_exp.Rdata")
load("D:/subject/exp/pData/pData85258.Rdata")

sign1<-Classifier(expr=GSE85258_exp, pData=pData85258,genelist=set1$Name)
sign2<-Classifier(expr=GSE85258_exp, pData=pData85258,genelist=set2$Name)
sign3<-Classifier(expr=GSE85258_exp, pData=pData85258,genelist=set3$Name)
sign4<-Classifier(expr=GSE85258_exp, pData=pData85258,genelist=set4$Name)
sign5<-Classifier(expr=GSE85258_exp, pData=pData85258,genelist=set5$Name)
sign6<-Classifier(expr=GSE85258_exp, pData=pData85258,genelist=set6$Name)
sign7<-Classifier(expr=GSE85258_exp, pData=pData85258,genelist=set7$Name)
sign9<-Classifier(expr=GSE85258_exp, pData=pData85258,genelist=set9$Name)
sign12<-Classifier(expr=GSE85258_exp, pData=pData85258,genelist=set12$Name)
sign13<-Classifier(expr=GSE85258_exp, pData=pData85258,genelist=set13$Name)
sign14<-Classifier(expr=GSE85258_exp, pData=pData85258,genelist=set14$Name)
sign17<-Classifier(expr=GSE85258_exp, pData=pData85258,genelist=set17$Name)
sign18<-Classifier(expr=GSE85258_exp, pData=pData85258,genelist=set18$Name)
sign19<-Classifier(expr=GSE85258_exp, pData=pData85258,genelist=set19$Name)
sign20<-Classifier(expr=GSE85258_exp, pData=pData85258,genelist=set20$Name)
sign21<-Classifier(expr=GSE85258_exp, pData=pData85258,genelist=set21$Name)
sign22<-Classifier(expr=GSE85258_exp, pData=pData85258,genelist=set22$Name)
sign24<-Classifier(expr=GSE85258_exp, pData=pData85258,genelist=set24$Name)
sign25<-Classifier(expr=GSE85258_exp, pData=pData85258,genelist=set25$Name)
