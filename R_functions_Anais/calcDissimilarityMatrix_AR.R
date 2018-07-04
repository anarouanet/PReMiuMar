calcDissimilarityMatrix_AR<-function(runInfoObj,onlyLS=FALSE,nSweeps1){


  directoryPath=NULL
  fileStem=NULL
  reportBurnIn=NULL
  nSweeps=NULL
  nFilter=NULL
  nSubjects=NULL
  nPredictSubjects=NULL
  nBurn=NULL
library(PReMiuM)

  for (i in 1:length(runInfoObj)) assign(names(runInfoObj)[i],runInfoObj[[i]])

  assign(names(runInfoObj)[4], min(runInfoObj[[4]],nSweeps1))

  fileName <- file.path(directoryPath,paste(fileStem,'_z.txt',sep=''))

  if (reportBurnIn) {
    recordedNBurn<-nBurn
  } else {
    recordedNBurn<-1
  }

  # Call the C++ to compute the dissimilarity matrix
  disSimList<-.Call('calcDisSimMat',fileName,nSweeps,recordedNBurn,nFilter,nSubjects,
                    nPredictSubjects, onlyLS, PACKAGE = 'PReMiuM')
  if (onlyLS){
    lsOptSweep<-disSimList$lsOptSweep
    disSimMatPred<-NULL
    disSimObj<-list('disSimRunInfoObj'=runInfoObj,'disSimMat'=NA,
                    'disSimMatPred'=NA,'lsOptSweep'=lsOptSweep,'onlyLS'=onlyLS)
  } else {
    disSimMat<-disSimList$disSimMat
    lsOptSweep<-disSimList$lsOptSweep
    disSimMatPred<-NULL
    if(nPredictSubjects>0){
      disSimMatPred<-disSimMat[(1+(nSubjects*(nSubjects-1)/2)):length(disSimMat)]
      disSimMat<-disSimMat[1:(nSubjects*(nSubjects-1)/2)]
    }
    disSimObj<-list('disSimRunInfoObj'=runInfoObj,'disSimMat'=disSimMat,
                    'disSimMatPred'=disSimMatPred,'lsOptSweep'=lsOptSweep,'onlyLS'=onlyLS)
  }

  return(disSimObj)
}
