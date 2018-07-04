plotRiskProfile_AR_chapter<-function(riskProfObj,outFile,showRelativeRisk=F,orderBy=NULL,whichClusters=NULL,whichCovariates=NULL,
                                     useProfileStar=F,riskLim=NULL,nSweeps1){

  riskProfClusObj=NULL
  clusObjRunInfoObj=NULL
  includeResponse=NULL
  yModel=NULL
  profileStar=NULL
  xModel=NULL
  whicCov=NULL
  nCategoriesY=NULL
  cluster=NULL
  prob=NULL
  meanProb=NULL
  fillColor=NULL
  lowerProb=NULL
  upperProb=NULL
  meanRisk=NULL
  lowerRisk=NULL
  upperRisk=NULL
  clusterSize=NULL
  mu=NULL
  meanMu=NULL
  lowerMu=NULL
  upperMu=NULL
  sigma=NULL
  meanSigma=NULL
  lowerSigma=NULL
  upperSigma=NULL
  weibullFixedShape=NULL
  nu=NULL
  meanNu=NULL
  lowerNu=NULL
  upperNu=NULL
  library("PReMiuM")
  library("grid")
  library("ggplot2")
  for (i in 1:length(riskProfObj)) assign(names(riskProfObj)[i],riskProfObj[[i]])
  for (i in 1:length(riskProfClusObj)) assign(names(riskProfClusObj)[i],riskProfClusObj[[i]])
  for (i in 1:length(clusObjRunInfoObj)) assign(names(clusObjRunInfoObj)[i],clusObjRunInfoObj[[i]])

  assign(names(clusObjRunInfoObj)[4],min(nSweeps1,clusObjRunInfoObj[[4]]))



  if (nClusters==1) stop("Cannot produce plots because only one cluster has been found.")

  plotRiskFlag <- 2
  chapter <- 1
  if(chapter == 1 && yModel=="Longitudinal"){
    includeResponse<-FALSE
    plotRiskFlag <- 0
  }else if(chapter == 1 && yModel=="MVN"){
    includeResponse<-TRUE
    plotRiskFlag <- 0
  }

  if(includeResponse){
    if(yModel=="Normal"){
      showRelativeRisk<-F
    }
    if(yModel=="Longitudinal"){
      plotRiskFlag <- 1
    }
  }
  if(useProfileStar){
    profile<-profileStar
  }
  if(!is.null(whichCovariates)){
    if (!is.numeric(whichCovariates)){
      whichCovariatesTmp<-vector()
      for (k in 1:length(whichCovariates)){
        whichCovariatesTmp[k]<-which(riskProfClusObj$clusObjRunInfoObj$covNames==whichCovariates[k])
      }
      whichCovariates<-whichCovariatesTmp
    }
    if(xModel=='Discrete'){
      profile<-profile[,,whichCovariates,]
      nCategories<-nCategories[whichCovariates]
      covNames<-covNames[whichCovariates]
      nCovariates<-length(whichCovariates)
    }else if(xModel=='Normal'){
      profile<-profile[,,whichCovariates]
      profileStdDev<-profileStdDev[,,whichCovariates,whichCovariates]
      covNames<-covNames[whichCovariates]
      nCovariates<-length(whichCovariates)
    }else if(xModel=='Mixed'){
      nDiscreteCovsAll <- nDiscreteCovs
      nContinuousCovsAll <- nContinuousCovs
      whichDiscreteCovs <- whichCovariates[which(whichCovariates<=nDiscreteCovs)]
      whichContinuousCovs <- whichCovariates[which(whichCovariates>nDiscreteCovs)]
      discreteCovs <- discreteCovs[whichDiscreteCovs]
      nDiscreteCovs <- length(discreteCovs)
      continuousCovs <- continuousCovs[whichContinuousCovs-nDiscreteCovsAll]
      nContinuousCovs <- length(continuousCovs)
      profilePhi<-profilePhi[,,whichDiscreteCovs,]
      nCategories<-nCategories[whichDiscreteCovs]
      profileMu<-profileMu[,,whichContinuousCovs-nDiscreteCovsAll]
      profileStdDev<-profileStdDev[,,whichContinuousCovs-nDiscreteCovsAll,whichContinuousCovs-nDiscreteCovsAll]
      covNames<-c(discreteCovs,continuousCovs)
      nCovariates<-length(covNames)

    }
  }

  png(outFile,width=1350,height=450)
  #pdf(outFile, width=20,height=5)
  orderProvided<-F

  if(!is.null(orderBy)){
    if(!includeResponse){
      if(orderBy!='Empirical'&&orderBy!='ClusterSize'&&!orderBy%in%covNames){
        if(is.numeric(orderBy)){
          if(length(orderBy)==nClusters){
            orderProvided<-T
            meanSortIndex<-orderBy
          }else{
            cat("Order vector provided not of same length as number of clusters. Reverting to default ordering.\n")
            orderBy<-NULL
          }
          orderBy<-NULL
        }
        #orderBy<-NULL
      }
    }else{
      if(orderBy!='Risk'&&orderBy!='Empirical'&&orderBy!='ClusterSize'&&!orderBy%in%covNames){
        if(is.numeric(orderBy)){
          if(length(orderBy)==nClusters){
            orderProvided<-T
            meanSortIndex<-orderBy
          }else{
            cat("Order vector provided not of same length as number of clusters. Reverting to default ordering.\n")
            orderBy<-NULL
          }
          orderBy<-NULL
        }
        #orderBy<-NULL
      }
    }

  }

  # Set up the layout for the plot
  ##//RJ remove 'Risk' from Longitudinal summary

  if(chapter==0){
    plotLayout<-grid.layout(ncol = nCovariates+plotRiskFlag, nrow = 6)
  }else{
    plotLayout<-grid.layout(ncol = nCovariates, nrow = 3)
  }

  grid.newpage()
  pushViewport(viewport(layout = plotLayout))
  if(!orderProvided){
    if(!is.null(risk)){
      if(is.null(orderBy)){
        # Default is to order by posterior theta risk
        # Compute the means
        orderStat<-apply(risk,2,median)
      }else{
        if(orderBy=='Risk'){
          orderStat<-apply(risk,2,median)
        }else if(orderBy=='Empirical'){
          orderStat<-empiricals
        }else if(orderBy=='ClusterSize'){
          orderStat<-clusterSizes
        }else{
          whichCov<-match(orderBy,covNames)
          if(xModel=='Normal'){
            orderStat<-apply(profile[,,whichCov],2,median)
          }else{
            # This assumes that there is some order to the categories
            # and then uses an expected value
            tmpMat<-profile[,,whichCov,1]
            if(nCategories[whichCov]>1){
              for(k in 2:nCategories[whichCov]){
                tmpMat<-tmpMat+k*profile[,,whichCov,k]
              }
            }
            orderStat<-apply(tmpMat,2,median)
          }
        }
      }
    }else{
      if(is.null(orderBy)){
        # Default is to order by empirical risk
        orderStat<-empiricals
      }else{
        if(orderBy=='Empirical'){
          orderStat<-empiricals
        }else if(orderBy=='ClusterSize'){
          orderStat<-clusterSizes
        }else{
          whichCov<-match(orderBy,covNames)
          if(xModel=='Normal'){
            orderStat<-apply(profile[,,whichCov],2,median)
          }else{
            # This assumes that there is some order to the categories
            # and then uses an expected value
            tmpMat<-profile[,,whichCov,1]
            if(nCategories[whichCov]>1){
              for(k in 2:(nCategories[whichCov])){
                tmpMat<-tmpMat+k*profile[,,whichCov,k]
              }
            }
            orderStat<-apply(tmpMat,2,median)
          }
        }
      }
    }
    # Sort into ascending mean size
    meanSortIndex<-order(orderStat,decreasing=F)
    meanSortIndex<-1:nClusters
  }
  if(includeResponse){
    # Reorder the risk matrix
    riskDim<-dim(risk)
    risk<-array(risk[,meanSortIndex,],dim=riskDim)
    if(showRelativeRisk){
      for(c in nClusters:1){
        risk[,c,]<-risk[,c,]/risk[,1,]
      }
    }
    # reorder the nu matrix
    if (yModel=="Survival"&&!weibullFixedShape){
      nuDim<-dim(nuArray)
      nuArray<-array(nuArray[,meanSortIndex],dim=nuDim)
    }
    ##//RJ reorder the L matrix
    if (yModel=="Longitudinal"){
      LDim<-dim(LArray)
      LArray<-array(LArray[,meanSortIndex,],dim=LDim)
    }
    if (yModel=="MVN"){
      MVNmuDim<-dim(MVNmuArray)
      MVNSigmaDim<-dim(MVNSigmaArray)
      MVNmuArray<-array(MVNmuArray[,meanSortIndex,],dim=MVNmuDim)
      MVNSigmaArray<-array(MVNSigmaArray[,meanSortIndex,],dim=MVNSigmaDim)
    }
  }

  # Reorder the cluster sizes
  clusterSizes<-clusterSizes[meanSortIndex]
  # Reorder the empiricals
  empiricals<-empiricals[meanSortIndex]
  meanEmpirical<-sum(empiricals*clusterSizes)/sum(clusterSizes)
  if(includeResponse){
    # Recompute the means and now also credible intervals
    riskMeans<-apply(risk,2,mean,trim=0.005)
    riskMean<-sum(riskMeans*clusterSizes)/sum(clusterSizes)
    riskLower<-apply(risk,2,quantile,0.05)
    riskUpper<-apply(risk,2,quantile,0.95)
    # The next line is to avoid outliers spoiling plot scales
    plotMax<-max(riskUpper)

    # Get the plot colors
    riskColor<-ifelse(riskLower>rep(riskMean,nClusters),"high",
                      ifelse(riskUpper<rep(riskMean,nClusters),"low","avg"))
    if (yModel=="Categorical"){
      riskDF<-data.frame("risk"=c(),"category"=c(),"cluster"=c(),"meanRisk"=c(),
                         "lowerRisk"=c(),"upperRisk"=c(),"fillColor"=c())
    } else {
      riskDF<-data.frame("risk"=c(),"cluster"=c(),"meanRisk"=c(),
                         "lowerRisk"=c(),"upperRisk"=c(),"fillColor"=c())
    }
    if (yModel=="Survival"&&!weibullFixedShape){
      # Recompute the means and now also credible intervals
      nuMeans<-apply(nuArray,2,mean,trim=0.005)
      nuMean<-sum(nuMeans*clusterSizes)/sum(clusterSizes)
      nuLower<-apply(nuArray,2,quantile,0.05)
      nuUpper<-apply(nuArray,2,quantile,0.95)

      nuDF<-data.frame("nu"=c(),"cluster"=c(),"meanNu"=c(),
                       "lowerNu"=c(),"upperNu"=c(),"fillColor"=c())
    }
    if (yModel=="Longitudinal"){
      # Recompute the means and now also credible intervals
      LMeans <- matrix(0,ncol=3,nrow=dim(LArray)[2])
      LMean <- c()
      for(i in 1:3){
        LMeans[,i]<-apply(LArray[,,i],2,mean,trim=0.005)
        LMean[i]<-sum(LMeans*clusterSizes)/sum(clusterSizes)
      }
    }
    if (yModel=="MVN"){
      # Recompute the means and now also credible intervals
      MVNmuMeans<-apply(MVNmuArray,2,mean,trim=0.005)
      MVNmuMean<-sum(MVNmuMeans*clusterSizes)/sum(clusterSizes)
      MVNmuLower<-apply(MVNmuArray,2,quantile,0.05)
      MVNmuUpper<-apply(MVNmuArray,2,quantile,0.95)
      MVNmuDF<-data.frame("MVNmu"=c(),"cluster"=c(),"meanMVNmu"=c(),
                          "lowerMVNmu"=c(),"upperMVNmu"=c(),"fillColor"=c())
    }

  }else{
    riskColor<-ifelse(empiricals>rep(meanEmpirical,length(empiricals)),"high",
                      ifelse(empiricals<rep(meanEmpirical,nClusters),"low","avg"))
  }

  if(is.null(whichClusters)){
    whichClusters<-1:nClusters
  }
  nClusters<-length(whichClusters)

  empiricalDF<-data.frame("empiricals"=c(),"meanEmpirical"=c(),"cluster"=c(),"fillColor"=c())
  sizeDF<-data.frame("clusterSize"=c(),"cluster"=c(),"fillColor"=c())
  # Restructure the data for plotting
  for(c in whichClusters){
    if(includeResponse){
      if (yModel=="Categorical"){
        plotRisk<-risk[,c,]
        nPoints<-dim(plotRisk)[1]
        for (k in 1:nCategoriesY){
          riskDF<-rbind(riskDF,data.frame("risk"=plotRisk[,k],
                                          "category"=rep(k,nPoints),
                                          "cluster"=rep(c,nPoints),
                                          "meanRisk"=rep(riskMean,nPoints),
                                          "lowerRisk"=rep(riskLower[c],nPoints),
                                          "upperRisk"=rep(riskUpper[c],nPoints),
                                          "fillColor"=rep(riskColor[c],nPoints)))
        }
      } else if (yModel=="Longitudinal"||yModel=="MVN"){##//RJ same as empirical
        plotRisk<-mean(risk[,c,1])
        nPoints<-length(plotRisk)
        riskDF<-rbind(riskDF,data.frame("risk"=plotRisk,"cluster"=rep(c,nPoints),
                                        "meanRisk"=rep(riskMean,nPoints),
                                        "lowerRisk"=rep(riskLower[c],nPoints),
                                        "upperRisk"=rep(riskUpper[c],nPoints),
                                        "fillColor"=rep(riskColor[c],nPoints)))
      } else {
        plotRisk<-risk[,c,]
        plotRisk<-plotRisk[plotRisk<plotMax]
        nPoints<-length(plotRisk)
        riskDF<-rbind(riskDF,data.frame("risk"=plotRisk,"cluster"=rep(c,nPoints),
                                        "meanRisk"=rep(riskMean,nPoints),
                                        "lowerRisk"=rep(riskLower[c],nPoints),
                                        "upperRisk"=rep(riskUpper[c],nPoints),
                                        "fillColor"=rep(riskColor[c],nPoints)))
      }
      if (yModel=="Survival"&&!weibullFixedShape){
        plotNu<-nuArray[,c]
        nPoints<-length(plotNu)
        nuDF<-rbind(nuDF,data.frame("nu"=plotNu,"cluster"=rep(c,nPoints),
                                    "meanNu"=rep(nuMean,nPoints),
                                    "lowerNu"=rep(nuLower[c],nPoints),
                                    "upperNu"=rep(nuUpper[c],nPoints),
                                    "fillColor"=rep(riskColor[c],nPoints)))
      }
    }
    empiricalDF<-rbind(empiricalDF,
                       data.frame("empiricals"=empiricals[c],
                                  "meanEmpirical"=meanEmpirical,"cluster"=c,"fillColor"=riskColor[c]))
    sizeDF<-rbind(sizeDF,
                  data.frame("clusterSize"=clusterSizes[c],"cluster"=c,"fillColor"=riskColor[c]))
  }

  if(includeResponse){
    if(yModel=='Categorical'){
      riskDF<-data.frame("prob"=c(),"cluster"=c(),"category"=c(),"meanProb"=c(),
                         "lowerProb"=c(),"upperProb"=c(),"fillColor"=c())
      for(k in 1:nCategoriesY){
        probMat<-risk[,,k]
        nPoints<-nrow(probMat)
        probMeans<-apply(probMat,2,mean)
        probMean<-sum(probMeans*clusterSizes)/sum(clusterSizes)
        probLower<-apply(probMat,2,quantile,0.05)
        probUpper<-apply(probMat,2,quantile,0.95)

        # Get the plot colors
        probColor<-ifelse(probLower>rep(probMean,length(probLower)),"high",
                          ifelse(probUpper<rep(probMean,length(probUpper)),"low","avg"))

        for(c in whichClusters){
          riskDF<-rbind(riskDF,data.frame("prob"=probMat[,c],"cluster"=rep(c,nPoints),
                                          "category"=rep(k-1,nPoints),
                                          "meanProb"=rep(probMean,nPoints),
                                          "lowerProb"=rep(probLower[c],nPoints),
                                          "upperProb"=rep(probUpper[c],nPoints),
                                          "fillColor"=rep(probColor[c],nPoints)))
          rownames(riskDF)<-seq(1,nrow(riskDF),1)

        }
      }

      plotObj<-ggplot(riskDF)
      plotObj<-plotObj+facet_wrap(~category,ncol=1,as.table=F,scales="free_y")
      plotObj<-plotObj+geom_hline(aes(x=as.factor(cluster),y=prob,yintercept=meanProb))
      plotObj<-plotObj+geom_boxplot(aes(x=as.factor(cluster),y=prob,fill=as.factor(fillColor)),outlier.size=0.5)
      plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=lowerProb,colour=as.factor(fillColor)),size=1.5)
      plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=upperProb,colour=as.factor(fillColor)),size=1.5)
      plotObj<-plotObj+
        scale_fill_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
        scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
        theme(legend.position="none")+labs(x="Cluster")+theme(axis.title.x=element_text(size=10))
      plotObj<-plotObj+labs(y="Probability")+theme(axis.title.y=element_text(size=10,angle=90))
      plotObj<-plotObj+labs(title=ifelse(showRelativeRisk,'Relative Risk','Risk'),plot.title=element_text(size=10))
      # Margin order is (top,right,bottom,left)
      plotObj<-plotObj+theme(plot.margin=unit(c(0.5,0.15,0.5,0.15),'lines'))+
        theme(plot.margin=unit(c(0,0,0,0),'lines'))
      print(plotObj,vp=viewport(layout.pos.row=1:6,layout.pos.col=2))
      #		}else if (yModel=="Survival"&&!weibullFixedShape){
      #			rownames(riskDF)<-seq(1,nrow(riskDF),by=1)
      #
      #			# Create the risk plot
      #			plotObj<-ggplot(riskDF)
      #			plotObj<-plotObj+geom_hline(aes(x=as.factor(cluster),y=risk,yintercept=meanRisk))
      #			plotObj<-plotObj+geom_boxplot(aes(x=as.factor(cluster),y=risk,fill=as.factor(fillColor)),outlier.size=0.5)
      #			plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=lowerRisk,colour=as.factor(fillColor)),size=1.5)
      #			plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=upperRisk,colour=as.factor(fillColor)),size=1.5)
      #			plotObj<-plotObj+scale_fill_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
      #				scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
      #				theme(legend.position="none")+
      #				labs(x="Cluster",y=ifelse(showRelativeRisk,'RR',
      #				ifelse(yModel=="Categorical"||yModel=="Bernoulli"||yModel=="Binomial","Probability","E[Y]")))
      #			plotObj<-plotObj+theme(axis.title.y=element_text(size=10,angle=90),axis.title.x=element_text(size=10))
      #			plotObj<-plotObj+labs(title=ifelse(showRelativeRisk,'Relative Risk','Risk'),plot.title=element_text(size=10))
      #			# Margin order is (top,right,bottom,left)
      #			plotObj<-plotObj+theme(plot.margin=unit(c(0,0,0,0),'lines'))+theme(plot.margin=unit(c(0.5,0.15,0.5,0.15),'lines'))
      #			print(plotObj,vp=viewport(layout.pos.row=1:3,layout.pos.col=2))
      #
      #			rownames(nuDF)<-seq(1,nrow(nuDF),by=1)
      #			# Create the nu plot
      #			plotObj<-ggplot(nuDF)
      #			plotObj<-plotObj+geom_hline(aes(x=as.factor(cluster),y=nu,yintercept=meanNu))
      #			plotObj<-plotObj+geom_boxplot(aes(x=as.factor(cluster),y=nu,fill=as.factor(fillColor)),outlier.size=0.5)
      #			plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=lowerNu,colour=as.factor(fillColor)),size=1.5)
      #			plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=upperNu,colour=as.factor(fillColor)),size=1.5)
      #			plotObj<-plotObj+scale_fill_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
      #				scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
      #				theme(legend.position="none")+
      #				labs(x="Cluster",y="Shape Parameter")
      #			plotObj<-plotObj+theme(axis.title.y=element_text(size=10,angle=90),axis.title.x=element_text(size=10))
      #			plotObj<-plotObj+labs(title="",plot.title=element_text(size=10))
      #			# Margin order is (top,right,bottom,left)
      #			plotObj<-plotObj+theme(plot.margin=unit(c(0,0,0,0),'lines'))+theme(plot.margin=unit(c(0.5,0.15,0.5,0.15),'lines'))
      #			print(plotObj,vp=viewport(layout.pos.row=4:6,layout.pos.col=2))
    } else if (yModel!='Longitudinal' && yModel!='MVN'){ ##//RJ
      rownames(riskDF)<-seq(1,nrow(riskDF),by=1)
      # Create the risk plot
      plotObj<-ggplot(riskDF)
      plotObj<-plotObj+geom_hline(aes(x=as.factor(cluster),y=risk,yintercept=meanRisk))
      plotObj<-plotObj+geom_boxplot(aes(x=as.factor(cluster),y=risk,fill=as.factor(fillColor)),outlier.size=0.5)
      if (!is.null(riskLim)) plotObj<-plotObj+coord_cartesian(ylim = riskLim)
      plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=lowerRisk,colour=as.factor(fillColor)),size=1.5)
      plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=upperRisk,colour=as.factor(fillColor)),size=1.5)
      plotObj<-plotObj+scale_fill_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
        scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
        theme(legend.position="none")+
        labs(x="Cluster",y=ifelse(showRelativeRisk,'RR',
                                  ifelse(yModel=="Categorical"||yModel=="Bernoulli"||yModel=="Binomial","Probability","E[Y]")))
      plotObj<-plotObj+theme(axis.title.y=element_text(size=10,angle=90),axis.title.x=element_text(size=10))
      plotObj<-plotObj+labs(title=ifelse(showRelativeRisk,'Relative Risk','Risk'),plot.title=element_text(size=10))
      # Margin order is (top,right,bottom,left)
      plotObj<-plotObj+theme(plot.margin=unit(c(0,0,0,0),'lines'))+theme(plot.margin=unit(c(0.5,0.15,0.5,0.15),'lines'))
      print(plotObj,vp=viewport(layout.pos.row=1:6,layout.pos.col=2))
    }
  }

  # Create a bar chart of cluster empiricals
  if((!is.null(yModel))){
    if(yModel!="Categorical"){
      plotObj<-ggplot(empiricalDF)
      plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=empiricals,colour=as.factor(fillColor)),size=3)
      plotObj<-plotObj+geom_hline(aes(x=as.factor(cluster),y=empiricals,yintercept=meanEmpirical))
      plotObj<-plotObj+scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
        theme(legend.position="none")
      plotObj<-plotObj+labs(title='Empirical Data',plot.title=element_text(size=10))
      plotObj<-plotObj+theme(axis.title.x=element_text(size=10),axis.title.y=element_text(size=10,angle=90))
      plotObj<-plotObj+
        labs(y=ifelse(yModel=="Bernoulli","Proportion of cases",
                      ifelse(yModel=="Binomial","Avg Proportion of occurrence",
                             ifelse(yModel=="Poisson","Avg Count",
                                    ifelse(yModel=="Survival","Avg Survival Time",
                                           ifelse(yModel=="Categorical","Avg Proportion of occurrence","Avg Y"))))),x="Cluster")
      plotObj<-plotObj+theme(plot.margin=unit(c(0,0,0,0),'lines'))+theme(plot.margin=unit(c(0.15,0.5,0.5,1),'lines'))
      print(plotObj,vp=viewport(layout.pos.row=1:3,layout.pos.col=1))
    }
  }

  if(chapter==0){
    # Create a bar chart of cluster sizes
    plotObj<-ggplot(sizeDF)
    plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=clusterSize,colour=as.factor(fillColor)),size=3)
    plotObj<-plotObj+scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+theme(legend.position="none")
    plotObj<-plotObj+labs(title="Size",plot.title=element_text(size=10))
    plotObj<-plotObj+theme(axis.title.x=element_text(size=10),axis.title.y=element_text(size=10,angle=90))
    plotObj<-plotObj+labs(y="No. of Subjects",x="Cluster")
    plotObj<-plotObj+theme(plot.margin=unit(c(0,0,0,0),'lines'))+theme(plot.margin=unit(c(0.15,0.5,0.5,1),'lines'))
    print(plotObj,vp=viewport(layout.pos.row=4:6,layout.pos.col=1))
  }

  # Loop over the covariates
  for(j in 1:nCovariates){
    if(xModel=='Discrete'){
      profileDF<-data.frame("prob"=c(),"cluster"=c(),"category"=c(),"meanProb"=c(),
                            "lowerProb"=c(),"upperProb"=c(),"fillColor"=c())
      for(k in 1:nCategories[j]){
        probMat<-profile[,meanSortIndex,j,k]
        nPoints<-nrow(probMat)
        probMeans<-apply(probMat,2,mean)
        probMean<-sum(probMeans*clusterSizes)/sum(clusterSizes)
        probLower<-apply(probMat,2,quantile,0.05)
        probUpper<-apply(probMat,2,quantile,0.95)
        # Get the plot colors
        probColor<-ifelse(probLower>rep(probMean,length(probLower)),"high",
                          ifelse(probUpper<rep(probMean,length(probUpper)),"low","avg"))


        for(c in whichClusters){
          profileDF<-rbind(profileDF,data.frame("prob"=probMat[,c],"cluster"=rep(c,nPoints),
                                                "category"=rep(k-1,nPoints),
                                                "meanProb"=rep(probMean,nPoints),
                                                "lowerProb"=rep(probLower[c],nPoints),
                                                "upperProb"=rep(probUpper[c],nPoints),
                                                "fillColor"=rep(probColor[c],nPoints)))
          rownames(profileDF)<-seq(1,nrow(profileDF),1)

        }
      }
      plotObj<-ggplot(profileDF)
      plotObj<-plotObj+facet_wrap(~category,ncol=1,as.table=F,scales="free_y")
      plotObj<-plotObj+geom_hline(aes(x=as.factor(cluster),y=prob,yintercept=meanProb))
      plotObj<-plotObj+geom_boxplot(aes(x=as.factor(cluster),y=prob,fill=as.factor(fillColor)),outlier.size=0.5)
      plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=lowerProb,colour=as.factor(fillColor)),size=1.5)
      plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=upperProb,colour=as.factor(fillColor)),size=1.5)
      plotObj<-plotObj+
        scale_fill_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
        scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
        theme(legend.position="none")+labs(x="Cluster")+theme(axis.title.x=element_text(size=10))
      if(j==1){
        plotObj<-plotObj+labs(y="Probability")+theme(axis.title.y=element_text(size=10,angle=90))
      }else{
        plotObj<-plotObj+theme(axis.title.y=element_blank())
      }
      plotObj<-plotObj+labs(title=covNames[j],plot.title=element_text(size=10))
      plotObj<-plotObj+theme(plot.margin=unit(c(0.5,ifelse(j==nCovariates,1,0),0.5,ifelse(j==1,0.5,0)),'lines'))+
        theme(plot.margin=unit(c(0,0,0,0),'lines'))
      print(plotObj,vp=viewport(layout.pos.row=1:6,layout.pos.col=j+plotRiskFlag))

    }else if(xModel=='Normal'){
      # Plot the means
      profileDF<-data.frame("mu"=c(),"cluster"=c(),"muMean"=c(),
                            "lowerMu"=c(),"upperMu"=c(),"fillColor"=c())
      muMat<-profile[,meanSortIndex,j]
      muMeans<-apply(muMat,2,mean)
      muMean<-sum(muMeans*clusterSizes)/sum(clusterSizes)
      muLower<-apply(muMat,2,quantile,0.05)
      muUpper<-apply(muMat,2,quantile,0.95)
      # The next line is to avoid outliers spoiling plot scales
      plotMax<-max(muUpper)
      plotMin<-min(muLower)

      # Get the plot colors
      muColor<-ifelse(muLower>rep(muMean,length(muLower)),"high",
                      ifelse(muUpper<rep(muMean,length(muUpper)),"low","avg"))
      muColor<-c("high","avg","low")
      for(c in whichClusters){
        plotMu<-muMat[,c]
        plotMu<-plotMu[plotMu<plotMax&plotMu>plotMin]
        nPoints<-length(plotMu)
        profileDF<-rbind(profileDF,data.frame("mu"=plotMu,"cluster"=rep(c,nPoints),
                                              "meanMu"=rep(muMean,nPoints),
                                              "lowerMu"=rep(muLower[c],nPoints),
                                              "upperMu"=rep(muUpper[c],nPoints),
                                              "fillColor"=rep(muColor[c],nPoints)))
      }
      rownames(profileDF)<-seq(1,nrow(profileDF),1)
      plotObj<-ggplot(profileDF)
      plotObj<-plotObj+geom_hline(aes(x=as.factor(cluster),y=mu,yintercept=meanMu))
      plotObj<-plotObj+geom_boxplot(aes(x=as.factor(cluster),y=mu,fill=as.factor(fillColor)),outlier.size=0.5)
      plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=lowerMu,colour=as.factor(fillColor)),size=1.5)
      plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=upperMu,colour=as.factor(fillColor)),size=1.5)
      plotObj<-plotObj+
        scale_fill_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
        scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
        theme(legend.position="none")+labs(x="Cluster")+theme(axis.title.x=element_text(size=10))
      if(j==1){
        plotObj<-plotObj+labs(y="Mean")+theme(axis.title.y=element_text(size=10,angle=90))
      }else{
        plotObj<-plotObj+theme(axis.title.y=element_blank())
      }
      plotObj<-plotObj+labs(title=covNames[j],plot.title=element_text(size=10))
      plotObj<-plotObj+
        theme(plot.margin=unit(c(0.5,ifelse(j==nCovariates,1,0),0.5,ifelse(j==1,0.5,0)),'lines'))+
        theme(plot.margin=unit(c(0,0,0,0),'lines'))

      print(plotObj,vp=viewport(layout.pos.row=1:3,layout.pos.col=j+plotRiskFlag))
      # Plot the variances
      profileDF<-data.frame("sigma"=c(),"cluster"=c(),"sigmaMean"=c(),
                            "lowerSigma"=c(),"upperSigma"=c(),"fillColor"=c())
      sigmaMat<-profileStdDev[,meanSortIndex,j,j]
      sigmaMeans<-apply(sigmaMat,2,mean)
      sigmaMean<-sum(sigmaMeans*clusterSizes)/sum(clusterSizes)
      sigmaLower<-apply(sigmaMat,2,quantile,0.05)
      sigmaUpper<-apply(sigmaMat,2,quantile,0.95)
      # The next line is to avoid outliers spoiling plot scales
      plotMax<-max(sigmaUpper)

      # Get the plot colors
      sigmaColor<-ifelse(sigmaLower>rep(sigmaMean,length(sigmaLower)),"high",
                         ifelse(sigmaUpper<rep(sigmaMean,length(sigmaUpper)),"low","avg"))
      for(c in whichClusters){
        plotSigma<-sigmaMat[,c]
        plotSigma<-plotSigma[plotSigma<plotMax]
        nPoints<-length(plotSigma)
        profileDF<-rbind(profileDF,data.frame("sigma"=plotSigma,"cluster"=rep(c,nPoints),
                                              "meanSigma"=rep(sigmaMean,nPoints),
                                              "lowerSigma"=rep(sigmaLower[c],nPoints),
                                              "upperSigma"=rep(sigmaUpper[c],nPoints),
                                              "fillColor"=rep(sigmaColor[c],nPoints)))
      }
      rownames(profileDF)<-seq(1,nrow(profileDF),1)

      plotObj<-ggplot(profileDF)
      plotObj<-plotObj+geom_hline(aes(x=as.factor(cluster),y=sigma,yintercept=meanSigma))
      plotObj<-plotObj+geom_boxplot(aes(x=as.factor(cluster),y=sigma,fill=as.factor(fillColor)),outlier.size=0.5)
      plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=lowerSigma,colour=as.factor(fillColor)),size=1.5)
      plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=upperSigma,colour=as.factor(fillColor)),size=1.5)
      plotObj<-plotObj+
        scale_fill_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
        scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
        theme(legend.position="none")+labs(x="Cluster")+theme(axis.title.x=element_text(size=10))
      if(j==1){
        plotObj<-plotObj+labs(y="Std Dev")+theme(axis.title.y=element_text(size=10,angle=90))
      }else{
        plotObj<-plotObj+theme(axis.title.y=element_blank())
      }
      plotObj<-plotObj+
        theme(plot.margin=unit(c(0.5,ifelse(j==nCovariates,1,0),0.5,ifelse(j==1,0.5,0)),'lines'))+
        theme(plot.margin=unit(c(0,0,0,0),'lines'))
      print(plotObj,vp=viewport(layout.pos.row=4:6,layout.pos.col=j+plotRiskFlag))

    }else if(xModel=='Mixed'){ # ici
      if (j<=nDiscreteCovs){
        profileDF<-data.frame("prob"=c(),"cluster"=c(),"category"=c(),"meanProb"=c(),
                              "lowerProb"=c(),"upperProb"=c(),"fillColor"=c())
        for(k in 1:nCategories[j]){
          if (nDiscreteCovs==1) {
            probMat<-profilePhi[,meanSortIndex,1,k]
          } else {
            probMat<-profilePhi[,meanSortIndex,j,k]
          }
          nPoints<-nrow(probMat)
          probMeans<-apply(probMat,2,mean)
          probMean<-sum(probMeans*clusterSizes)/sum(clusterSizes)
          probLower<-apply(probMat,2,quantile,0.05)
          probUpper<-apply(probMat,2,quantile,0.95)

          # Get the plot colors
          probColor<-ifelse(probLower>rep(probMean,length(probLower)),"high",
                            ifelse(probUpper<rep(probMean,length(probUpper)),"low","avg"))


          for(c in whichClusters){
            profileDF<-rbind(profileDF,data.frame("prob"=probMat[,c],"cluster"=rep(c,nPoints),
                                                  "category"=rep(k-1,nPoints),
                                                  "meanProb"=rep(probMean,nPoints),
                                                  "lowerProb"=rep(probLower[c],nPoints),
                                                  "upperProb"=rep(probUpper[c],nPoints),
                                                  "fillColor"=rep(probColor[c],nPoints)))
            rownames(profileDF)<-seq(1,nrow(profileDF),1)

          }
        }
        profileDF2<-profileDF[which(profileDF$category==1),]

        covNames<-c("Gender","Education","APOE4")
        plotObj<-ggplot(profileDF2)
        #plotObj<-plotObj+facet_wrap(~category,ncol=1,as.table=F,scales="free_y")
        plotObj<-plotObj+geom_hline(aes(x=as.factor(cluster),y=prob,yintercept=meanProb))
        plotObj<-plotObj+geom_boxplot(aes(x=as.factor(cluster),y=prob,fill=as.factor(cluster)),outlier.size=0.5)
        plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=lowerProb,colour=as.factor(fillColor)),size=1.5)
        plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=upperProb,colour=as.factor(fillColor)),size=1.5)
        plotObj<-plotObj+
          #scale_fill_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
          #scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
          theme(legend.position="none")+labs(x="Clusters")+theme(axis.title.x=element_text(size=19))
        if(j==1){
          plotObj<-plotObj+labs(y="Probability")+theme(axis.title.y=element_text(size=20,angle=90))
        }else{
          plotObj<-plotObj+theme(axis.title.y=element_blank())
        }
        plotObj<-plotObj+labs(title=covNames[j],plot.title=element_text(size=10))+
          theme(plot.title = element_text(size=19))
        plotObj<-plotObj+theme(plot.margin=unit(c(0.5,ifelse(j==nCovariates,0.5,0),0.5,ifelse(j==1,0.5,0)),'lines'))

        #+theme(plot.margin=unit(c(0,0,1,1),'lines'))
        print(plotObj,vp=viewport(layout.pos.row=1:3,layout.pos.col=j+plotRiskFlag))
      } else {
        # Plot the means
        profileDF<-data.frame("mu"=c(),"cluster"=c(),"muMean"=c(),
                              "lowerMu"=c(),"upperMu"=c(),"fillColor"=c())
        if (nContinuousCovs==1){
          muMat<-profileMu[,meanSortIndex,1]
        } else {
          muMat<-profileMu[,meanSortIndex,(j-nDiscreteCovs)]
        }
        muMeans<-apply(muMat,2,mean)
        muMean<-sum(muMeans*clusterSizes)/sum(clusterSizes)
        muLower<-apply(muMat,2,quantile,0.05)
        muUpper<-apply(muMat,2,quantile,0.95)
        # The next line is to avoid outliers spoiling plot scales
        plotMax<-max(muUpper)
        plotMin<-min(muLower)

        # Get the plot colors
        muColor<-ifelse(muLower>rep(muMean,length(muLower)),"high",
                        ifelse(muUpper<rep(muMean,length(muUpper)),"low","avg"))
        for(c in whichClusters){
          plotMu<-muMat[,c]
          plotMu<-plotMu[plotMu<plotMax&plotMu>plotMin]
          nPoints<-length(plotMu)
          profileDF<-rbind(profileDF,data.frame("mu"=plotMu,"cluster"=rep(c,nPoints),
                                                "meanMu"=rep(muMean,nPoints),
                                                "lowerMu"=rep(muLower[c],nPoints),
                                                "upperMu"=rep(muUpper[c],nPoints),
                                                "fillColor"=rep(muColor[c],nPoints)))
        }
        rownames(profileDF)<-seq(1,nrow(profileDF),1)

        covNames[4:9]<-c("s.Ventricles","s.Hippocampus",
                         "s.Entorhinal","s.Fusiform","s.MidTemp","s.WholeBrain")

        plotObj<-ggplot(profileDF)#+theme_bw()
        plotObj<-plotObj+geom_hline(aes(x=as.factor(cluster),y=mu,yintercept=meanMu))
        plotObj<-plotObj+geom_boxplot(aes(x=as.factor(cluster),y=mu,fill=as.factor(cluster)),outlier.size=0.5)
        plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=lowerMu,colour=as.factor(fillColor)),size=1.5)
        plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=upperMu,colour=as.factor(fillColor)),size=1.5)
        plotObj<-plotObj+
          #scale_fill_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
          #scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
          theme(legend.position="none")+labs(x="Clusters")+theme(axis.title.x=element_text(size=19))
        if(j==4){
          plotObj<-plotObj+labs(y="Mean")+theme(axis.title.y=element_text(size=20,angle=90))
        }else{
          plotObj<-plotObj+theme(axis.title.y=element_blank())
        }
        plotObj<-plotObj+labs(title=covNames[j],plot.title=element_text(size=10))+
          theme(plot.title = element_text(size=19))
        plotObj<-plotObj+
          theme(plot.margin=unit(c(0.5,ifelse(j==nCovariates,0.5,0),0.5,ifelse(j==4,0.5,0)),'lines'))
        #+ theme(plot.margin=unit(c(0,1,1,1),'lines'))

        print(plotObj,vp=viewport(layout.pos.row=1:3,layout.pos.col=j+plotRiskFlag))

        # Plot the variances
        if(chapter==0){

          profileDF<-data.frame("sigma"=c(),"cluster"=c(),"sigmaMean"=c(),
                                "lowerSigma"=c(),"upperSigma"=c(),"fillColor"=c())
          if (nContinuousCovs==1){
            sigmaMat<-profileStdDev[,meanSortIndex,1,1]
          } else {
            sigmaMat<-profileStdDev[,meanSortIndex,(j-nDiscreteCovs),(j-nDiscreteCovs)]
          }
          sigmaMeans<-apply(sigmaMat,2,mean)
          sigmaMean<-sum(sigmaMeans*clusterSizes)/sum(clusterSizes)
          sigmaLower<-apply(sigmaMat,2,quantile,0.05)
          sigmaUpper<-apply(sigmaMat,2,quantile,0.95)
          # The next line is to avoid outliers spoiling plot scales
          plotMax<-max(sigmaUpper)

          # Get the plot colors
          sigmaColor<-ifelse(sigmaLower>rep(sigmaMean,length(sigmaLower)),"high",
                             ifelse(sigmaUpper<rep(sigmaMean,length(sigmaUpper)),"low","avg"))
          for(c in whichClusters){
            plotSigma<-sigmaMat[,c]
            plotSigma<-plotSigma[plotSigma<plotMax]
            nPoints<-length(plotSigma)
            profileDF<-rbind(profileDF,data.frame("sigma"=plotSigma,"cluster"=rep(c,nPoints),
                                                  "meanSigma"=rep(sigmaMean,nPoints),
                                                  "lowerSigma"=rep(sigmaLower[c],nPoints),
                                                  "upperSigma"=rep(sigmaUpper[c],nPoints),
                                                  "fillColor"=rep(sigmaColor[c],nPoints)))
          }
          rownames(profileDF)<-seq(1,nrow(profileDF),1)
          plotObj<-ggplot(profileDF)
          plotObj<-plotObj+geom_hline(aes(x=as.factor(cluster),y=sigma,yintercept=meanSigma))
          plotObj<-plotObj+geom_boxplot(aes(x=as.factor(cluster),y=sigma,fill=as.factor(fillColor)),outlier.size=0.5)
          plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=lowerSigma,colour=as.factor(fillColor)),size=1.5)
          plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=upperSigma,colour=as.factor(fillColor)),size=1.5)
          plotObj<-plotObj+
            scale_fill_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
            scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
            theme(legend.position="none")+labs(x="Cluster")+theme(axis.title.x=element_text(size=10))
          if(j==1){
            plotObj<-plotObj+labs(y="Std Dev")+theme(axis.title.y=element_text(size=10,angle=90))
          }else{
            plotObj<-plotObj+theme(axis.title.y=element_blank())
          }
          plotObj<-plotObj+
            theme(plot.margin=unit(c(0.5,ifelse(j==nCovariates,1,0),0.5,ifelse(j==1,0.5,0)),'lines'))+
            theme(plot.margin=unit(c(0,0,0,0),'lines'))
          print(plotObj,vp=viewport(layout.pos.row=4:6,layout.pos.col=j+plotRiskFlag))
        }
      }

    }
  }
  dev.off()

  if(yModel=='MVN'){
    png(paste(strsplit(outFile,"\\.")[[1]][1],'-MVN.png',sep=""),width=800,height=400)
    plotLayout<-grid.layout(ncol = nOutcomes, nrow = 3)
    grid.newpage()
    pushViewport(viewport(layout = plotLayout))

    for(j in 1:nOutcomes){
      # Plot the means
      profileDF<-data.frame("mu"=c(),"cluster"=c(),"muMean"=c(),
                            "lowerMu"=c(),"upperMu"=c(),"fillColor"=c())
      muMat<-MVNmuArray[,,j]
      muMeans<-apply(muMat,2,mean)
      muMean<-sum(muMeans*clusterSizes)/sum(clusterSizes)
      muLower<-apply(muMat,2,quantile,0.05)
      muUpper<-apply(muMat,2,quantile,0.95)
      # The next line is to avoid outliers spoiling plot scales
      plotMax<-max(muUpper)
      plotMin<-min(muLower)

      # Get the plot colors
      muColor<-ifelse(muLower>rep(muMean,length(muLower)),"high",
                      ifelse(muUpper<rep(muMean,length(muUpper)),"low","avg"))
      for(c in whichClusters){
        plotMu<-muMat[,c]
        plotMu<-plotMu[plotMu<plotMax&plotMu>plotMin]
        nPoints<-length(plotMu)
        profileDF<-rbind(profileDF,data.frame("mu"=plotMu,"cluster"=rep(c,nPoints),
                                              "meanMu"=rep(muMean,nPoints),
                                              "lowerMu"=rep(muLower[c],nPoints),
                                              "upperMu"=rep(muUpper[c],nPoints),
                                              "fillColor"=rep(muColor[c],nPoints)))
      }


      if(chapter==1){
        print("transformation!!!")
        profileDF[c("mu"        ,"meanMu"    ,"lowerMu"    ,"upperMu")]<-
          profileDF[c("mu"        ,"meanMu"    ,"lowerMu"    ,"upperMu")]* 10.880659
      }


      if(chapter==1)
        outcome<-c("Random intercept", "Random slope")
      rownames(profileDF)<-seq(1,nrow(profileDF),1)

      plotObj<-ggplot(profileDF)
      plotObj<-plotObj+geom_hline(aes(x=as.factor(cluster),y=mu,yintercept=meanMu))
      plotObj<-plotObj+geom_boxplot(aes(x=as.factor(cluster),y=mu,fill=as.factor(cluster)),outlier.size=0.5)
      plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=lowerMu,colour=as.factor(fillColor)),size=1.5)
      plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=upperMu,colour=as.factor(fillColor)),size=1.5)
      plotObj<-plotObj+
        #scale_fill_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
        #scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
        theme(legend.position="none")+labs(x="Clusters")+theme(axis.title.x=element_text(size=20))
      if(j==1){
        plotObj<-plotObj+labs(y="Mean")+theme(axis.title.y=element_text(size=20,angle=90))
      }else{
        plotObj<-plotObj+theme(axis.title.y=element_blank())
      }
      plotObj<-plotObj+labs(title=outcome[j],plot.title=element_text(size=20))+
        theme(plot.title = element_text(size=20))
      plotObj<-plotObj+
        theme(plot.margin=unit(c(0.5,ifelse(j==nOutcomes,1,0),0.5,ifelse(j==1,0.5,0)),'lines'))
      #+ theme(plot.margin=unit(c(0,0,0,0),'lines'))
      print(plotObj,vp=viewport(layout.pos.row=1:3,layout.pos.col=j))


      # Plot the variances
      profileDF<-data.frame("sigma"=c(),"cluster"=c(),"sigmaMean"=c(),
                            "lowerSigma"=c(),"upperSigma"=c(),"fillColor"=c())
      # MVNSigmaArray is a vector of the lower triangular covariance matrix, ordered by row.
      # The index of the variance is the j-th triangular number.
      sigmaMat<-sqrt(MVNSigmaArray[,,j*(j+1)/2])
      sigmaLower<-apply(sigmaMat,2,quantile,0.05)
      sigmaUpper<-apply(sigmaMat,2,quantile,0.95)
      sigmaMeans<-apply(sigmaMat,2,mean)
      sigmaMean<-sum(sigmaMeans*clusterSizes)/sum(clusterSizes)
      # The next line is to avoid outliers spoiling plot scales
      plotMax<-max(sigmaUpper)

      if(chapter==0){

        # Get the plot colors
        sigmaColor<-ifelse(sigmaLower>rep(sigmaMean,length(sigmaLower)),"high",
                           ifelse(sigmaUpper<rep(sigmaMean,length(sigmaUpper)),"low","avg"))
        for(c in whichClusters){
          plotSigma<-sigmaMat[,c]
          plotSigma<-plotSigma[plotSigma<plotMax]
          nPoints<-length(plotSigma)
          profileDF<-rbind(profileDF,data.frame("sigma"=plotSigma,"cluster"=rep(c,nPoints),
                                                "meanSigma"=rep(sigmaMean,nPoints),
                                                "lowerSigma"=rep(sigmaLower[c],nPoints),
                                                "upperSigma"=rep(sigmaUpper[c],nPoints),
                                                "fillColor"=rep(sigmaColor[c],nPoints)))
        }
        rownames(profileDF)<-seq(1,nrow(profileDF),1)
        plotObj<-ggplot(profileDF)
        plotObj<-plotObj+geom_hline(aes(x=as.factor(cluster),y=sigma,yintercept=meanSigma))
        plotObj<-plotObj+geom_boxplot(aes(x=as.factor(cluster),y=sigma,fill=as.factor(fillColor)),outlier.size=0.5)
        plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=lowerSigma,colour=as.factor(fillColor)),size=1.5)
        plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=upperSigma,colour=as.factor(fillColor)),size=1.5)
        plotObj<-plotObj+
          scale_fill_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
          scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
          theme(legend.position="none")+labs(x="Cluster")+theme(axis.title.x=element_text(size=10))
        if(j==1){
          plotObj<-plotObj+labs(y="Std Dev")+theme(axis.title.y=element_text(size=10,angle=90))
        }else{
          plotObj<-plotObj+theme(axis.title.y=element_blank())
        }
        plotObj<-plotObj+
          theme(plot.margin=unit(c(0.5,ifelse(j==nOutcomes,1,0),0.5,ifelse(j==1,0.5,0)),'lines'))+
          theme(plot.margin=unit(c(0,0,0,0),'lines'))
        print(plotObj,vp=viewport(layout.pos.row=4:6,layout.pos.col=j))
      }
    }
    dev.off()
  }
  return(meanSortIndex)
}
