library(PReMiuMar)


set.seed(1)
print('Equispaced measurements')
clusterSummary_quad <- clusSummaryLongitudinalDiscrete_Quadratic_kernel()
#clusterSummary <- clusSummaryLongitudinalDiscrete()
inputs_quad <- generateSampleDataFile(clusterSummary_quad)
head(inputs$inputLongData)
hist(inputs_quad$inputLongData$outcome)

par(mfrow=c(1,1))
plot(inputs_quad$inputLongData$time[inputs_quad$inputLongData$ID==1],
     inputs_quad$inputLongData$outcome[inputs_quad$inputLongData$ID==1],
     type='l',xlab='Age',ylab='outcome',ylim=c(min(inputs_quad$inputLongData$outcome),max(inputs_quad$inputLongData$outcome)))
for(i in 2:max(inputs_quad$inputLongData$ID))
  lines(inputs_quad$inputLongData$time[inputs_quad$inputLongData$ID==i],
        inputs_quad$inputLongData$outcome[inputs_quad$inputLongData$ID==i])

hyper <- list()
hyper$muLSignal<-0
hyper$sigmaLSignal<-1
hyper$muLLengthscale<--3
hyper$sigmaLLengthscale<-1
hyper$muLNoise<--4
hyper$sigmaLNoise<-1
hyper$muLTime<-2
hyper$sigmaLTime<-1

runInfoObj_ROB_long <- PReMiuMar::profRegr(yModel = inputs$yModel, xModel = inputs$xModel,
                                           nSweeps = 1000, nBurn = 0,#nSweeps = 5000, nBurn = 500,
                                           data = inputs$inputData,
                                           output = "output/kernel_Quad",
                                           covNames = inputs$covNames[1],
                                           outcomeT = NA, nClusInit = 20,
                                           run = TRUE, longData=inputs$inputLongData,
                                           hyper=hyper,
                                           kernel = "Quadratic",
                                           seed = 1)

runInfoObj=runInfoObj_ROB_long
nSweeps=runInfoObj$nSweeps

dissimObj <- calcDissimilarityMatrix(runInfoObj)
clusObj <- calcOptimalClustering(dissimObj)
myheatDissMat(dissimObj)
#clusObj$nOutcomes<-2
riskProfileObj <- PReMiuMar::calcAvgRiskAndProfile(clusObj)

clusterOrderObj_50 <- plotRiskProfile(riskProfileObj, "output/kernel_Quad.png")


myheatDissMat <- function (dissimObj, main = NULL, xlab = NULL, ylab = NULL)
{
  nSbj <- dissimObj$disSimRunInfoObj$nSubjects
  col.labels <- c("0", "0.5", "1")
  colours <- colorRampPalette(c("white", "black"))(10)
  dissMat <- vec2mat(dissimObj$disSimMat, nrow = nSbj)
  heatmap(1 - dissMat, keep.dendro = FALSE, symm = TRUE, Rowv = NULL,
          labRow = FALSE, labCol = FALSE, margins = c(0.5,0.5),
          col = colours, main = main, xlab = xlab, ylab = ylab)
  plotrix::color.legend(0.89, 0.7, 0.94, 1, legend = col.labels, colours,
                        gradient = "y", align = "rb",cex=1.3)
}

