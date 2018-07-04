#output_mod<-function(runInfoObj){
dev.off()

runInfoObj=runInfoObj_combine_RE_age10
  nSweeps=runInfoObj$nSweeps

  dissimObj <- calcDissimilarityMatrix_AR(runInfoObj,nSweeps1=nSweeps)
  myheatDissMat(dissimObj)
  clusObj <- calcOptimalClustering_AR(dissimObj,nSweeps1=nSweeps)
  #clusObj$nOutcomes<-3
  riskProfileObj <- calcAvgRiskAndProfile_AR(clusObj,nSweeps1=nSweeps)
  #
  clusterOrderObj <- plotRiskProfile(riskProfileObj, "MVN/Combine_RE_age10_quad/output.png")

dev.off()#}
  dissimObj <- calcDissimilarityMatrix(runInfoObj)
  myheatDissMat(dissimObj)
  clusObj <- calcOptimalClustering(dissimObj)
  clusObj$nOutcomes<-4
  riskProfileObj <- calcAvgRiskAndProfile(clusObj)
  clusterOrderObj <- plotRiskProfile_AR(riskProfileObj, "MVN/Combine_RE_age/output.png",nSweeps1=nSweeps)
  #clusObj$nOutcomes<-4

  file<-"/Users/anais/Documents/Cambridge/code/R/Premium/Premium_on_ADNI/new_PReMiuM/output_200s_1000s_cov_X_ICV/long_z.txt"
  outputZ<-read.table(file,header=F)
  nSubjects <- dim(outputZ)[2]
  nSweeps<-dim(outputZ)[1]



  dissim<- matrix(0,nSubjects,nSubjects)
  for(i in 1:nSubjects){
    if(i<nSubjects){
      for(j in (i+1):nSubjects){
        dissim[i,j]<-length(which(outputZ[,i]==outputZ[,j]))
      }
    }
  }
  dissim <- dissim+t(dissim)
  diag(dissim)<-nSweeps

  recordedNBurn<-0
  nFilter<-0
  onlyLS<-FALSE
  nPredictSubjects<-0
  disSimList<-.Call('calcDisSimMat',file,nSweeps,recordedNBurn,nFilter,nSubjects,
                    nPredictSubjects, onlyLS, PACKAGE = 'PReMiuM')


  col.labels <- c("0", "0.5", "1")
  colours <- colorRampPalette(c("white", "black"))(10)
  dissMat <- vec2mat(file, nrow = nSbj)
  heatmap(1 - dissim, keep.dendro = FALSE, symm = TRUE, Rowv = NULL,
          labRow = FALSE, labCol = FALSE, margins = c(4.5, 4.5))#,
          #col = colours, main = main, xlab = xlab, ylab = ylab)
  plotrix::color.legend(0.95, 0.7, 1, 1, legend = col.labels, colours,
                        gradient = "y", align = "rb")
