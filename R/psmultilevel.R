#' Propensity Score Weighting Estimator for Clustered Data
#'
#' \code{psmultilevel} is a general function, making calls to
#' \code{multi.ate} and \code{multi.att} depending on the estimand of interest.
#' It implements a variety of algorithms for calculating propensity score
#' weighting estimators and their standard errors for treatment effects in a
#' clustered data setting. Here, the treatment is assigned on an individual
#' level and not cluster level. This function has an optional threshold for the minimum
#' number of units in each treatment group of each cluster - clusters that don't
#' meet this requirement are removed before the propensity scores are
#' estimated.
#' \cr
#' \cr
#' \code{psmultilevel} offers four weighting
#' estimators: marginal weighting, clustered weighting, stratification, and
#' doubly-robust. Standard errors include closed-form formula when available
#' (Lunceford and Davidian 2004), as well as values obtained by bootstrapping.
#'
#' @param data Data containing covariates except for the outcome, the treatment,
#'   and the cluster identification.
#' @param Y A vector containing the outcome of interest. This can be either
#'   binary or continuous, but will need to be specified for methods involving
#'   outcome models.
#' @param Z A vector indicating the observations which are in the treated and
#'   control group. This is binary-coded: 1 = Treated, 0 = Untreated.
#' @param cluster A vector of factored cluster identification.
#' @param y.formula Outcome formula, required for Doubly-Robust. The treatment variable should be called "Z" and the cluster variable should be named "cluster".
#' @param pscore.formula Propensity score formula. The treatment variable should be called "Z" and the cluster variable should be named "cluster".
#' @param estimand A character string for the estimand, either "ATE", the sample
#'   average treatment effect, or "ATT", the sample average treatment effect for
#'   the treated.
#' @param method A character string for the method of interest. Options: "marwt"
#'   for marginal weighting, "clwt" for clustered weighting, "stratification"
#'   for stratification/ subclassification, and "DR" for doubly-robust.
#' @param modeltype A character string identifying the model of the cluster-specific
#' effects for both the propensity score and the outcome. Options: "fixed" for a fixed effects model, and
#'   "random" for a random effects model. Default: fixed effects model. Note that for random
#'   effects models, function \code{glmer} of package \code{lme4} is utilized, and for large models it can be very slow to fit. To speed up, we add two options to the \code{glmer} defaults: \code{nAGQ=0} and \code{control = glmerControl(optimizer = "nloptwrap")}. See the documentation of \code{glmer} for the details about these options.
#' @param y.type A character string classifying the outcome type. Options: "binary"/
#'   "continuous".
#' @param se.report A logical flag for whether the standard error for the
#'   estimator should be reported.
#' @param nboot The number of bootstrap samples to be drawn. Default value is
#'   75.
#' @param nsub A scalar denoting the number of subclasses to be used for
#'   stratification. Default value is 6.
#' @param ps.cut A vector of length 2 containing the lower and upper limit for the
#'   propensity scores cut-offs. Observations beyond these limits are removed. Default
#'   values are (0.05,0.95).
#' @param cl.cutoff A scalar denoting the minimum number of observations required in both
#'   the treatment and the control groups. Clusters with fewer observations than
#'   the threshold in either group will be removed before the estimation
#'   process. Default value is 0.
#' @return A list with the elements: \item{estimate}{Value of the estimator of
#'   interest.} \item{se}{The standard error of the estimator. Closed-form formula
#'   is available for ATE Doubly-Robust (Lunceford and Davidian 2004). Bootstrapping is implemented
#'   to calculate SEs for the rest of the estimators.}
#'
#' @seealso \code{\link{multi.ate}}, \code{\link{multi.att}}
#' @references Lunceford, J. K. and Davidian, M. (2004), Stratification and
#'   weighting via the propensity score in estimation of causal treatment
#'   effects: a comparative study. Statist. Med., 23: 2937–2960.
#'   doi:10.1002/sim.1903
#'
#' @examples
#' y.formula.re<-function(){
#' return (Y~Z+X1+X2+(1|cluster))
#' }
#'
#' ps.formula.re<-function(){
#' return(Z~ X1+X2+ (1|cluster))
#' }
#' data<-simdata[,c(1:2)]
#'
#' psmultilevel(data=data,Y=simdata$Y,Z=simdata$Z,cluster=simdata$cluster,y.formula=NULL,
#' pscore.formula=ps.formula.re(),estimand='ATE',method='marwt',modeltype='random',y.type=NULL,
#' se.report=FALSE,ps.cut=c(0.05,0.95),cl.cutoff=0)
#'
#' \dontrun{
#' psmultilevel(data=data,Y=simdata$Y,Z=simdata$Z,cluster=simdata$cluster,
#' y.formula = y.formula.re(),pscore.formula=ps.formula.re(),estimand='ATE',
#' method='DR',modeltype='random',y.type='continuous',se.report=TRUE,
#' ps.cut=c(0.05,0.95),cl.cutoff=3)
#' }
#' @export

psmultilevel<-function(data,Y,Z,cluster,y.formula=NULL,pscore.formula,
                       estimand,method,modeltype='fixed',y.type,se.report=FALSE,nboot=75,
                       nsub=6,ps.cut=c(0.05,0.95),cl.cutoff=0){

  results<-NULL

  if (estimand=='ATT'){

    results<-multi.att(data=data,Y=Y,Z=Z,cluster=cluster,y.formula=y.formula,pscore.formula=pscore.formula,
                       method=method,modeltype=modeltype,y.type=y.type,se.report=se.report,
                       nboot=nboot,nsub=nsub,ps.cut=ps.cut,cl.cutoff=cl.cutoff)

  } #End ATT section
  else if(estimand=='ATE'){

    results<-multi.ate(data=data,Y=Y,Z=Z,cluster=cluster,y.formula=y.formula,pscore.formula=pscore.formula,
                       method=method,modeltype=modeltype,y.type=y.type,se.report=se.report,
                       nboot=nboot,nsub=nsub,ps.cut=ps.cut,cl.cutoff=cl.cutoff)

  }
  return(results)
}



#' ATE Propensity Score Weighting Estimator for Clustered Data
#'
#' \code{multi.ate} implements a variety of algorithms for calculating
#' propensity score weighting estimators and their standard errors for the
#' sample average treatment effect (ATE) in a clustered data setting. Here, the
#' treatment is assigned on an individual level and not cluster level.
#' \code{multi.ate} has an optional threshold for the minimum number of units in
#' each treatment group of each cluster - clusters that don't meet this
#' requirement are removed before the propensity scores are estimated.
#'\cr
#'\cr
#' \code{multi.ate} offers four weighting estimators: marginal weighting,
#' clustered weighting, stratification, and doubly-robust. Standard errors
#' include closed-form formula when available (Lunceford and Davidian 2004), as
#' well as values obtained by bootstrapping. \code{psmultilevel} makes calls to
#' \code{multi.ate} as needed.
#'
#' @param data Data containing covariates except for the outcome, the treatment,
#'   and the cluster identification.
#' @param Y A vector containing the outcome of interest. This can be either
#'   binary or continuous, but will need to be specified for methods involving
#'   outcome models.
#' @param Z A vector indicating the observations which are in the treated and
#'   control group. This is binary-coded: 1 = Treated, 0 = Untreated.
#' @param cluster A vector of factored cluster identification.
#' @param y.formula Outcome formula, required for Doubly-Robust. The treatment variable should be called "Z" and the cluster variable should be called "cluster".
#' @param pscore.formula Propensity score formula. The treatment variable should be called "Z" and the cluster variable should be called "cluster".
#' @param method A character string for the method of interest. Options: "marwt"
#'   for marginal weighting, "clwt" for clustered weighting, "stratification"
#'   for stratification/ subclassification, and "DR" for doubly-robust.
#' @param modeltype A character string identifying the model of the cluster-specific
#' effects for both the propensity score and the outcome. Options: "fixed" for a fixed effects model, and
#'   "random" for a random effects model. Default: fixed effects model. Note that for random
#'   effects models, function \code{glmer} of package \code{lme4} is utilized, and for large models it can be very slow to fit. To speed up, we add two options to the \code{glmer} defaults: \code{nAGQ=0} and \code{control = glmerControl(optimizer = "nloptwrap")}. See the documentation of \code{glmer} for the details about these options.
#' @param y.type A character string classifying the outcome type. Options: "binary"/
#'   "continuous".
#' @param se.report A logical flag for whether the standard error for the
#'   estimator should be reported.
#' @param nboot The number of bootstrap samples to be drawn. Default value is
#'   75.
#' @param nsub A scalar denoting the number of subclasses to be used for
#'   stratification. Default value is 6.
#' @param ps.cut A vector of length 2 containing the lower and upper limit for the
#'   propensity scores cut-offs. Observations beyond these limits are removed. Default
#'   values are (0.05,0.95).
#' @param cl.cutoff A scalar denoting the minimum number of observations required in both
#'   the treatment and the control groups. Clusters with fewer observations than
#'   the threshold in either group will be removed before the estimation
#'   process. Default value is 0.
#' @return A list with the elements: \item{estimate}{Value of the estimator of
#'   interest.} \item{se}{The standard error of the estimator. Closed-form formula
#'   is available for ATE DR (Lunceford and Davidian 2004). Bootstrapping is implemented
#'   to calculate SE for the rest of the estimators.}
#'
#' @references Lunceford, J. K. and Davidian, M. (2004), Stratification and
#'   weighting via the propensity score in estimation of causal treatment
#'   effects: a comparative study. Statist. Med., 23: 2937–2960.
#'   doi:10.1002/sim.1903
#'
#' @seealso \code{\link{multi.att}}
#'
#' @examples
#' y.formula.re<-function(){
#' return (Y~Z+X1+X2+(1|cluster))
#' }
#'
#' ps.formula.re<-function(){
#' return(Z~ X1+X2+ (1|cluster))
#' }
#' data<-simdata[,c(1:2)]
#'
#' multi.ate(data=data,Y=simdata$Y,Z=simdata$Z,cluster=simdata$cluster,y.formula=NULL,
#' pscore.formula=ps.formula.re(),method='marwt',modeltype='random',y.type=NULL,
#' se.report=FALSE,ps.cut=c(0.05,0.95),cl.cutoff=0)
#
#' \dontrun{
#' #Method = Cluster weighting, obtain bootstrapped SE
#' multi.ate(data=data,Y=simdata$Y,Z=simdata$Z,cluster=simdata$cluster,
#' y.formula = y.formula.re(),pscore.formula=ps.formula(),method='clwt',modeltype='random',
#' se.report=TRUE,nboot=75,ps.cut=c(0.05,0.95),cl.cutoff=0)
#' }
#' @export
multi.ate<-function(data,Y,Z,cluster,y.formula=NULL,pscore.formula,
                    method,modeltype='fixed',y.type,se.report=FALSE,nboot=75,
                    nsub=6,ps.cut=c(0.05,0.95),cl.cutoff=0){
  data<-cbind(Y,Z,data,cluster=cluster)

  #Remove clusters with too few observations
  clus.table<-as.data.frame(table(data$cluster))
  colnames(clus.table)<-c('Cluster','Size')
  clus.table<-clus.table[clus.table$Size!=0,]

  temp<-aggregate(data$Z,by=list(data$cluster),FUN=sum)
  colnames(temp)<-c('Cluster','N1')

  clus.table<-merge(clus.table,temp,by='Cluster')
  remove(temp)
  clus.table$N0<-clus.table$Size - clus.table$N1
  clus.table<-clus.table[clus.table$N1>cl.cutoff & clus.table$N0>cl.cutoff,]
  clus.table$Cluster<-factor(clus.table$Cluster)

  data<-data[data$cluster %in% clus.table$Cluster,] #
  data$cluster<-factor(data$cluster)

  #Estimate propensity scores
  if(modeltype=='fixed'){
    ps.reg<-glm(pscore.formula,data=data,family=binomial(link='logit'))
    data$pscore<-predict(ps.reg,type='response')

  }else if(modeltype=='random'){
    ps.reg<-glmer(pscore.formula,data=data,family=binomial(link='logit'),
                  nAGQ=0,
                  control=glmerControl(optimizer = "nloptwrap"))
    data$pscore<-fitted(ps.reg)
  }

  estimate<-NULL
  se<-NULL
  if(method=='marwt'){
    datawt<-data[data$pscore>ps.cut[1] & data$pscore<ps.cut[2],]

    #
    #Calculate ATE weights
    datawt$wts<-with(datawt,((1/pscore)^Z)*((1/(1-pscore))^(1-Z)))

    #Calculate ATE estimator
    estimate<-with(datawt,sum(Y*wts*Z)/sum(wts*Z)- sum(Y*wts*(1-Z))/sum(wts*(1-Z)))

    if(se.report==TRUE){ #ATE marginal weighting: bootstrap for SE
       estimate.boot<-rep(0,nboot)

        #Resample
        numclus<-nlevels(data$cluster)
        for (n in 1:nboot){
          #Sample the clusters with replacement
          clus.sam<-sample(clus.table$Cluster,replace=TRUE)
          data.boot<-NULL
          for(m in 1:numclus){
            data.cur<-data[data$cluster==clus.sam[m],]
            data.cur$cluster<-as.factor(m)
            data.boot<-rbind(data.boot,data.cur)
          }

          if(modeltype=='fixed'){
            reg.ps<-glm(pscore.formula,data=data.boot,family=binomial(link='logit'))
            data.boot$pscore<-predict(reg.ps,type='response')

          }else if(modeltype=='random'){
            reg.ps<-glmer(pscore.formula,data=data.boot,family=binomial(link='logit'),
                          nAGQ=0,
                          control=glmerControl(optimizer = "nloptwrap"))
            data.boot$pscore<-predict(reg.ps,type='response')
          }
          #cutoff
          datawt<-data.boot[data.boot$pscore>ps.cut[1] & data.boot$pscore<ps.cut[2],]

          #Calculate ATE weights
          datawt$wts<-with(datawt,((1/pscore)^Z)*((1/(1-pscore))^(1-Z)))

          #Calculate ATE
          estimate.boot[n]<-with(datawt,sum(Y*wts*Z)/sum(wts*Z))-
            with(datawt,sum(Y*wts*(1-Z))/sum(wts*(1-Z)))

        }
        se<-sd(estimate.boot)

    }
  } else if(method=='clwt'){
    datawt<-data[data$pscore>ps.cut[1] & data$pscore<ps.cut[2],]

    #Calculate ATE weights
    datawt$wts<-with(datawt,((1/pscore)^Z)*((1/(1-pscore))^(1-Z)))

    #Calculate total weights in each cluster
    totwt<-aggregate(datawt$wts,by=list(datawt$cluster),FUN=sum)$x

    #Fixed effects ps
    cluster.list<-unique(datawt$cluster) #Get list of unique clusters
    wt.clus<-rep(0,length(cluster.list))
    for (i in 1:length(cluster.list)){
      datawt.temp<-datawt[datawt$cluster==cluster.list[i],]

      #Estimate ATE for each cluster
      wt.clus[i]<-with(datawt.temp,sum(Y*wts*Z)/sum(wts*Z) -sum(Y*wts*(1-Z))/sum(wts*(1-Z)))
    }
    tempmat<-cbind(wt.clus,totwt)
    tempmat<-tempmat[!rowSums(!is.finite(tempmat)),]      #Exclude clusters with only treatment or control
    #Calculate ATE
    estimate<-sum(tempmat[,1]*tempmat[,2])/sum(tempmat[,2])

    if(se.report==TRUE){
        estimate.boot<-rep(0,nboot)

        #Resample
        numclus<-nlevels(data$cluster)
        for (n in 1:nboot){
          #Sample the clusters with replacement
          clus.sam<-sample(clus.table$Cluster,replace=TRUE)
          data.boot<-NULL
          for(m in 1:numclus){
            data.cur<-data[data$cluster==clus.sam[m],]
            data.cur$cluster<-as.factor(m)
            data.boot<-rbind(data.boot,data.cur)
          }

          if(modeltype=='fixed'){
            reg.ps<-glm(pscore.formula,data=data.boot,family=binomial(link='logit'))
            data.boot$pscore<-predict(reg.ps,type='response')

          }else if(modeltype=='random'){
            reg.ps<-glmer(pscore.formula,data=data.boot,family=binomial(link='logit'),
                          nAGQ=0,
                          control=glmerControl(optimizer = "nloptwrap"))
            data.boot$pscore<-predict(reg.ps,type='response')
          }
          #cutoff
          datawt<-data.boot[data.boot$pscore>ps.cut[1] & data.boot$pscore<ps.cut[2],]


          #Calculate ATE weights
          datawt$wts<-with(datawt,((1/pscore)^Z)*((1/(1-pscore))^(1-Z)))

          totwt<-aggregate(datawt$wts,by=list(datawt$cluster),FUN=sum)$x

          #Fixed effects ps
          cluster.list<-unique(datawt$cluster) #Get list of unique clusters
          wt.clus<-rep(0,length(cluster.list))
          for (i in 1:length(cluster.list)){
            idx.temp<-which(datawt$cluster==cluster.list[i])

            #Estimate ATE for each cluster
            wt.clus[i]<-with(datawt,sum(Y[idx.temp]*wts[idx.temp]*Z[idx.temp])/sum(wts[idx.temp]*Z[idx.temp]) -
                               sum(Y[idx.temp]*wts[idx.temp]*(1-Z[idx.temp]))/sum(wts[idx.temp]*(1-Z[idx.temp])))
          }
          tempmat<-cbind(wt.clus,totwt)
          tempmat<-tempmat[!rowSums(!is.finite(tempmat)),]      #Exclude clusters with only treatment or control
          #Calculate ATE
          estimate.boot[n]<-sum(tempmat[,1]*tempmat[,2])/sum(tempmat[,2])

        }
        se<-sd(estimate.boot)

    }   #End se.report==TRUE

  } else if (method=='stratification'){
    qtile<-quantile(data$pscore,probs=seq(0,1,length=nsub+1))
    qtile[1]<-qtile[1]-0.0001
    qtile[nsub+1]<-qtile[nsub+1]+0.0001


    #Classify individuals based on their estimated pscore
    subclass<-rep(length(data$Y),0)
    for (l in 1:nsub){
      idx.temp<- which(data$pscore>qtile[l] & data$pscore<=qtile[l+1])
      subclass[idx.temp]<-l
    }

    #Calculate the within-subclass estimate of ATE
    subtable<-table(subclass)
    sub.est<-rep(0,length(subtable))
    for (m in 1:length(subtable)){
      data.temp<-data[subclass==rownames(subtable)[m],]
      sub.est[m]<-with(data.temp,mean(Y[Z==1])-mean(Y[Z==0]))
    }

    tempmat<-cbind(sub.est,subtable)
    tempmat<-tempmat[!rowSums(!is.finite(tempmat)),]    #Exclude subclasses with only treatment or controls

    #Take their average weighted by the number of units in each subclass
    estimate<-sum(tempmat[,1]*tempmat[,2])/sum(tempmat[,2])

    if(se.report==TRUE){
        estimate.boot<-rep(0,nboot)

        #Resample
        numclus<-nlevels(data$cluster)
        for (n in 1:nboot){
          #Sample the clusters with replacement
          clus.sam<-sample(clus.table$Cluster,replace=TRUE)
          data.boot<-NULL
          for(m in 1:numclus){
            data.cur<-data[data$cluster==clus.sam[m],]
            data.cur$cluster<-as.factor(m)
            data.boot<-rbind(data.boot,data.cur)
          }

          if(modeltype=='fixed'){
            reg.ps<-glm(pscore.formula,data=data.boot,family=binomial(link='logit'))
            data.boot$pscore<-predict(reg.ps,type='response')

          }else if(modeltype=='random'){
            reg.ps<-glmer(pscore.formula,data=data.boot,family=binomial(link='logit'),
                          nAGQ=0,
                          control=glmerControl(optimizer = "nloptwrap"))
            data.boot$pscore<-predict(reg.ps,type='response')
          }

          qtile.boot<-quantile(data.boot$pscore,probs=seq(0,1,length=nsub+1))
          qtile.boot[1]<-qtile.boot[1]-0.0001
          qtile.boot[nsub+1]<-qtile.boot[nsub+1]+0.0001


          #Classify individuals based on their estimated pscore
          subclass<-rep(length(data.boot$Y),0)
          for (l in 1:nsub){
            idx.temp<- which(data.boot$pscore>qtile.boot[l] & data.boot$pscore<=qtile.boot[l+1])
            subclass[idx.temp]<-l
          }

          #Calculate the within-subclass estimate of ATE
          subtable<-table(subclass)
          sub.est<-rep(0,length(subtable))
          for (m in 1:length(subtable)){
            data.temp<-data.boot[subclass==rownames(subtable)[m],]
            sub.est[m]<-with(data.temp,mean(Y[Z==1])-mean(Y[Z==0]))
          }

          tempmat<-cbind(sub.est,subtable)
          tempmat<-tempmat[!rowSums(!is.finite(tempmat)),]    #Exclude subclasses with only treatment or controls

          #Take their average weighted by the number of  units in each subclass
          estimate.boot[n]<-sum(tempmat[,1]*tempmat[,2])/sum(tempmat[,2])

        }
        se<-sd(estimate.boot)

    } #End se.report==TRUE
  } else if(method=='DR'){
    datawt<-data[data$pscore>ps.cut[1] & data$pscore<ps.cut[2],]
    data.temp<-datawt
    data.temp$Z<-0

    if(modeltype=='fixed'){
      if(y.type=='binary'){
        ####Fixed Effects
        regDR.y<-glm(y.formula,data=datawt,family=binomial(link='logit'))
        y.0<-predict(regDR.y,newdata=data.temp,type='response')
        data.temp$Z<-1
        y.1<-predict(regDR.y,newdata=data.temp,type='response')

      }else if (y.type=='continuous'){
        regDR.y<-lm(y.formula,data=datawt)
        y.0<-predict(regDR.y,newdata=data.temp)
        data.temp$Z<-1
        y.1<-predict(regDR.y,newdata=data.temp)
      }
    }else if(modeltype=='random'){
      if(y.type=='binary'){
        regDR.y<-glmer(y.formula,
                       data=datawt,family=binomial(link='logit'),
                       nAGQ=0,
                       control=glmerControl(optimizer = "nloptwrap"))
        y.0<-predict(regDR.y,newdata=data.temp,type='response')
        data.temp$Z<-1
        y.1<-predict(regDR.y,newdata=data.temp,type='response')
      } else if(y.type=='continuous'){
        regDR.y<-lmer(y.formula,data=datawt)
        y.0<-predict(regDR.y,newdata=data.temp)
        data.temp$Z<-1
        y.1<-predict(regDR.y,newdata=data.temp)
      }
    }

    #Impute potential outcomes
    datawt$estDR.y0<-y.0
    datawt$estDR.y1<-y.1
    #ATE Estimator
    #DR estimators

    datawt$drest<-with(datawt,((Z*Y/pscore) - (Z-pscore)*estDR.y1/pscore)-
                         ((1-Z)*Y/(1-pscore)+(Z-pscore)*estDR.y0/(1-pscore)))
    estimate<-mean(datawt$drest)

    if (se.report==TRUE){
      se<-sqrt(sum((datawt$drest - estimate)^2)/(dim(datawt)[1]^2))
    }
  } #End method = DR
  results<-list('estimate' = estimate,'se'=se)
  return(results)
}

#' ATT Propensity Score Weighting Estimator for Clustered Data
#'
#' \code{multi.att} implements a variety of algorithms for calculating
#' propensity score weighting estimators and their standard errors for the
#' sample average treatment effect for the treated (ATT) in a clustered data
#' setting. Here, the treatment is assigned on an individual level and not
#' cluster level. \code{multi.att} has an optional threshold for the minimum
#' number of units in each treatment group of each cluster - clusters that don't
#' meet this requirement are removed before the propensity scores are
#' estimated.
#' \cr
#' \cr
#' \code{multi.att} offers four weighting estimators: marginal weighting,
#' clustered weighting, stratification, and doubly-robust. Standard
#' errors are obtained by bootstrapping. \code{psmultilevel} makes
#' calls to \code{multi.att} as needed.
#'
#' @param data Data containing covariates except for the outcome, the treatment,
#'   and the cluster identification.
#' @param Y A vector containing the outcome of interest. This can be either
#'   binary or continuous, but will need to be specified for methods involving
#'   outcome models.
#' @param Z A vector indicating the observations which are in the treated and
#'   control group. This is binary-coded: 1 = Treated, 0 = Untreated.
#' @param cluster A vector of factored cluster identification.
#' @param y.formula Outcome formula, required for Doubly-Robust. The treatment variable should be called "Z", the cluster variable should be called "cluster", and the outcome
#' @param pscore.formula Propensity score formula. The treatment variable should be called "Z" and the cluster variable should be called "cluster".
#' @param method A character string for the method of interest. Options: "marwt"
#'   for marginal weighting, "clwt" for clustered weighting, "stratification"
#'   for stratification/ subclassification, and "DR" for doubly-robust.
#' @param modeltype A character string identifying the model of the cluster-specific
#' effects for both the propensity score and the outcome. Options: "fixed" for a fixed effects model, and
#'   "random" for a random effects model. Default: fixed effects model. Note that for random
#'   effects models, function \code{glmer} of package \code{lme4} is utilized, and for large models it can be very slow to fit. To speed up, we add two options to the \code{glmer} defaults: \code{nAGQ=0} and \code{control = glmerControl(optimizer = "nloptwrap")}. See the documentation of \code{glmer} for the details about these options.
#' @param y.type A character string classifying the outcome type. Options: "binary"/
#'   "continuous".
#' @param se.report A logical flag for whether the standard error for the
#'   estimator should be reported.
#' @param nboot The number of bootstrap samples to be drawn. Default value is
#'   75.
#' @param nsub A scalar denoting the number of subclasses to be used for
#'   stratification. Default value is 6.
#' @param ps.cut A vector of length 2 containing the lower and upper limit for the
#'   propensity scores cut-offs. Observations beyond these limits are removed. Default
#'   values are (0.05,0.95).
#' @param cl.cutoff A scalar denoting the minimum number of observations required in both
#'   the treatment and the control groups. Clusters with fewer observations than
#'   the threshold in either group will be removed before the estimation
#'   process. Default value is 0.
#
#' @return A list with the elements: \item{estimate}{Value of the estimator of
#'   interest.} \item{se}{The standard error of the estimator, obtained by bootstrapping.}
#'
#' @seealso \code{\link{multi.ate}}
#'
#' @examples
#' y.formula.re<-function(){
#' return (Y~Z+X1+X2+(1|cluster))
#' }
#'
#' ps.formula.re<-function(){
#' return(Z~ X1+X2+ (1|cluster))
#' }
#' data<-simdata[,c(1:2)]
#'
#' multi.att(data=data,Y=simdata$Y,Z=simdata$Z,cluster=simdata$cluster,y.formula=NULL,
#' pscore.formula=ps.formula.re(),method='marwt',modeltype='random',y.type=NULL,
#' se.report=FALSE,ps.cut=c(0.05,0.95),cl.cutoff=0)
#
#' \dontrun{
#' #Method = Cluster weighting, obtain bootstrapped SE
#' multi.att(data=data,Y=simdata$Y,Z=simdata$Z,cluster=simdata$cluster,
#' y.formula = y.formula.re(),pscore.formula=ps.formula(),method='clwt',modeltype='random',
#' se.report=TRUE,nboot=75,ps.cut=c(0.05,0.95),cl.cutoff=0)
#' }
#' @export

multi.att<-function(data,Y,Z,cluster,y.formula=NULL,pscore.formula,
                            method,modeltype='fixed',y.type,se.report=FALSE,nboot=75,
                            nsub=6,ps.cut=c(0.05,0.95),cl.cutoff=0){
  data<-cbind(Y,Z,data,cluster=cluster)

  #Remove clusters with too few observations
  clus.table<-as.data.frame(table(data$cluster))
  colnames(clus.table)<-c('Cluster','Size')
  clus.table<-clus.table[clus.table$Size!=0,]

  temp<-aggregate(data$Z,by=list(data$cluster),FUN=sum)
  colnames(temp)<-c('Cluster','N1')

  clus.table<-merge(clus.table,temp,by='Cluster')
  remove(temp)
  clus.table$N0<-clus.table$Size - clus.table$N1
  clus.table<-clus.table[clus.table$N1>cl.cutoff & clus.table$N0>cl.cutoff,]
  clus.table$Cluster<-factor(clus.table$Cluster)

  data<-data[data$cluster %in% clus.table$Cluster,] #
  data$cluster<-factor(data$cluster)


  #Estimate propensity scores
  if(modeltype=='fixed'){
    ps.reg<-glm(pscore.formula,data=data,family=binomial(link='logit'))
    data$pscore<-predict(ps.reg,type='response')

  }else if(modeltype=='random'){
    ps.reg<-glmer(pscore.formula,data=data,family=binomial(link='logit'),
                  nAGQ=0,
                  control=glmerControl(optimizer = "nloptwrap"))
    data$pscore<-predict(ps.reg,type='response')
  }
  estimate<-NULL
  se<-NULL

  if (method=='marwt'){
    datawt<-data[data$pscore>ps.cut[1] & data$pscore<ps.cut[2],]

    #
    #Calculate ATT weights
    datawt$wts<-with(datawt,(1^Z)*((pscore/(1-pscore))^(1-Z)))

    #Calculate ATT
    estimate<-with(datawt,sum(Y*wts*Z)/sum(wts*Z)- sum(Y*wts*(1-Z))/sum(wts*(1-Z)))

    if(se.report==TRUE){
        estimate.boot<-rep(0,nboot)

        #Resample
        numclus<-nlevels(data$cluster)
        for (n in 1:nboot){
          #Sample the clusters with replacement
          clus.sam<-sample(clus.table$Cluster,replace=TRUE)
          data.boot<-NULL
          for(m in 1:numclus){
            data.cur<-data[data$cluster==clus.sam[m],]
            data.cur$cluster<-as.factor(m)
            data.boot<-rbind(data.boot,data.cur)
          }
          if(modeltype=='fixed'){
            reg.ps<-glm(pscore.formula,data=data.boot,family=binomial(link='logit'))
            data.boot$pscore<-predict(reg.ps,type='response')

          }else if(modeltype=='random'){
            reg.ps<-glmer(pscore.formula,data=data.boot,family=binomial(link='logit'),
                          nAGQ=0,
                          control=glmerControl(optimizer = "nloptwrap"))
            data.boot$pscore<-predict(reg.ps,type='response')
          }
          #cutoff
          datawt<-data.boot[data.boot$pscore>ps.cut[1] & data.boot$pscore<ps.cut[2],]

          #Calculate ATT weights
          datawt$wts<-with(datawt,(1^Z)*((pscore/(1-pscore))^(1-Z)))

          #Calculate ATT
          estimate.boot[n]<-with(datawt,sum(Y*wts*Z)/sum(wts*Z)- sum(Y*wts*(1-Z))/sum(wts*(1-Z)))
        }
        se<-sd(estimate.boot)

    }

  } else if(method=='clwt'){
    datawt<-data[data$pscore>ps.cut[1] & data$pscore<ps.cut[2],]

    #Calculate ATT weights
    datawt$wts<-with(datawt,(1^Z)*((pscore/(1-pscore))^(1-Z)))

    #Calculate total weights in each cluster
    totwt<-aggregate(datawt$wts,by=list(datawt$cluster),FUN=sum)$x

    #Fixed effects ps
    cluster.list<-unique(datawt$cluster) #Get list of unique clusters
    wt.clus<-rep(0,length(cluster.list))
    for (i in 1:length(cluster.list)){
      datawt.temp<-datawt[datawt$cluster==cluster.list[i],]

      #Estimate ATT for each cluster
      wt.clus[i]<-with(datawt.temp,sum(Y*wts*Z)/sum(wts*Z) -sum(Y*wts*(1-Z))/sum(wts*(1-Z)))
    }
    tempmat<-cbind(wt.clus,totwt)
    tempmat<-tempmat[!rowSums(!is.finite(tempmat)),]      #Exclude clusters with only treatment or control
    #Calculate ATT
    estimate<-sum(tempmat[,1]*tempmat[,2])/sum(tempmat[,2])

    if(se.report==TRUE){
        estimate.boot<-rep(0,nboot)

        #Resample
        numclus<-nlevels(data$cluster)
        for (n in 1:nboot){
          #Sample the clusters with replacement
          clus.sam<-sample(clus.table$Cluster,replace=TRUE)
          data.boot<-NULL
          for(m in 1:numclus){
            data.cur<-data[data$cluster==clus.sam[m],]
            data.cur$cluster<-as.factor(m)
            data.boot<-rbind(data.boot,data.cur)
          }

          if(modeltype=='fixed'){
            reg.ps<-glm(pscore.formula,data=data.boot,family=binomial(link='logit'))
            data.boot$pscore<-predict(reg.ps,type='response')

          }else if(modeltype=='random'){
            reg.ps<-glmer(pscore.formula,data=data.boot,family=binomial(link='logit'),
                          nAGQ=0,
                          control=glmerControl(optimizer = "nloptwrap"))
            data.boot$pscore<-predict(reg.ps,type='response')
          }
          #cutoff
          datawt<-data.boot[data.boot$pscore>ps.cut[1] & data.boot$pscore<ps.cut[2],]

          #Calculate ATT weights
          datawt$wts<-with(datawt,(1^Z)*((pscore/(1-pscore))^(1-Z)))

          #Calculate total weights in each cluster
          totwt<-aggregate(datawt$wts,by=list(datawt$cluster),FUN=sum)$x

          #Fixed effects ps
          cluster.list<-unique(datawt$cluster) #Get list of unique clusters
          wt.clus<-rep(0,length(cluster.list))
          for (i in 1:length(cluster.list)){
            idx.temp<-which(datawt$cluster==cluster.list[i])

            #Estimate ATT for each cluster
            wt.clus[i]<-with(datawt,sum(Y[idx.temp]*wts[idx.temp]*Z[idx.temp])/sum(wts[idx.temp]*Z[idx.temp]) -
                               sum(Y[idx.temp]*wts[idx.temp]*(1-Z[idx.temp]))/sum(wts[idx.temp]*(1-Z[idx.temp])))
          }
          tempmat<-cbind(wt.clus,totwt)
          tempmat<-tempmat[!rowSums(!is.finite(tempmat)),]      #Exclude clusters with only treatment or control
          #Calculate ATT
          estimate.boot[n]<-sum(tempmat[,1]*tempmat[,2])/sum(tempmat[,2])
        }
        se<-sd(estimate.boot)
    }
  }else if(method=='stratification'){
    qtile<-quantile(data$pscore[data$Z==1],probs=seq(0,1,length=nsub+1))
    qtile[1]<-qtile[1]-0.01
    qtile[nsub+1]<-qtile[nsub+1]+0.01


    #Classify individuals based on their estimated pscore
    subclass<-rep(0,nrow(data))
    for (l in 1:nsub){
      idx.temp<- which(data$pscore>qtile[l] & data$pscore<=qtile[l+1])
      subclass[idx.temp]<-l
    }
    #Remove units that are not classified
    idxdrop<-which(is.na(subclass))

    if(length(idxdrop)>0){
      subclass<-subclass[-idxdrop]
      data1<-data[-idxdrop,]
    }else{
      data1<-data
    }


    #Calculate the within-subclass estimate of ATT
    subtable<-table(subclass,data1$Z)
    subtable<-subtable[!(apply(subtable, 1, function(y) any(y == 0))),] #Remove subclasses with only 1 group

    sub.est<-rep(0,nrow(subtable))
    for (m in 1:nrow(subtable)){
      data.temp<-data1[subclass==as.numeric(rownames(subtable)[m]),]
      sub.est[m]<-with(data.temp,mean(Y[Z==1])-mean(Y[Z==0]))
    }

    #Take their average weighted by the number of treated units in each subclass
    estimate<-sum(sub.est*subtable[,2])/sum(subtable[,2])

    if(se.report==TRUE){
        estimate.boot<-rep(0,nboot)

        #Resample
        numclus<-nlevels(data$cluster)
        for (n in 1:nboot){
          #Sample the clusters with replacement
          clus.sam<-sample(clus.table$Cluster,replace=TRUE)
          data.boot<-NULL
          for(m in 1:numclus){
            data.cur<-data[data$cluster==clus.sam[m],]
            data.cur$cluster<-as.factor(m)
            data.boot<-rbind(data.boot,data.cur)
          }

          if(modeltype=='fixed'){
            reg.ps<-glm(pscore.formula,data=data.boot,family=binomial(link='logit'))
            data.boot$pscore<-predict(reg.ps,type='response')

          }else if(modeltype=='random'){
            reg.ps<-glmer(pscore.formula,data=data.boot,family=binomial(link='logit'),
                          nAGQ=0,
                          control=glmerControl(optimizer = "nloptwrap"))
            data.boot$pscore<-predict(reg.ps,type='response')
          }

          qtile<-quantile(data.boot$pscore[data.boot$Z==1],probs=seq(0,1,length=nsub+1))
          qtile[1]<-qtile[1]-0.01
          qtile[nsub+1]<-qtile[nsub+1]+0.01


          #Classify individuals based on their estimated pscore
          subclass<-rep(dim(data.boot)[1],0)
          for (l in 1:nsub){
            idx.temp<- which(data.boot$pscore>qtile[l] & data.boot$pscore<=qtile[l+1])
            subclass[idx.temp]<-l
          }

          idxdrop<-which(is.na(subclass))
          if(length(idxdrop)>0){
            subclass<-subclass[-idxdrop]
            data.boot1<-data.boot[-idxdrop,]
          }else{
            data.boot1<-data.boot
          }

          #Calculate the within-subclass estimate of ATT
          subtable<-table(subclass,data.boot1$Z)
          subtable<-subtable[!(apply(subtable, 1, function(y) any(y == 0))),] #Remove subclasses with only 1 group

          sub.est<-rep(0,length(subtable[,2]))
          for (m in 1:length(subtable[,2])){
            data.temp<-data.boot1[subclass==rownames(subtable)[m],]
            sub.est[m]<-with(data.temp,mean(Y[Z==1])-mean(Y[Z==0]))
          }

          #Take their average weighted by the number of treated units in each subclass
          estimate.boot[n]<-sum(sub.est*subtable[,2])/sum(subtable[,2])
        }
        se<-sd(estimate.boot)

    }

  }else if(method=='DR'){
    datawt<-data[data$pscore>ps.cut[1] & data$pscore<ps.cut[2],]
    data.temp<-datawt
    data.temp$Z<-0

    if(modeltype=='fixed'){
      if(y.type=='binary'){
        ####Fixed Effects
        regDR.y<-glm(y.formula,data=datawt,family=binomial(link='logit'))
        y.0<-predict(regDR.y,newdata=data.temp,type='response')

      }else if (y.type=='continuous'){
        regDR.y<-lm(y.formula,data=datawt)
        y.0<-predict(regDR.y,newdata=data.temp)
      }
    }else if(modeltype=='random'){
      if(y.type=='binary'){
        regDR.y<-glmer(y.formula,
                       data=datawt,family=binomial(link='logit'),
                       nAGQ=0,
                       control=glmerControl(optimizer = "nloptwrap"))
        y.0<-predict(regDR.y,newdata=data.temp,type='response')
      } else if(y.type=='continuous'){
        regDR.y<-lmer(y.formula,data=datawt)
        y.0<-predict(regDR.y,newdata=data.temp)
      }
    }


    #Impute potential outcomes
    datawt$estDR.y0<-y.0

    #ATT Estimator
    estimate<-with(datawt,sum((Y - estDR.y0)*(Z-pscore)/(1-pscore))/sum(Z))


    if(se.report==TRUE){
        estimate.boot<-rep(0,nboot)

        #Resample
        numclus<-nlevels(data$cluster)
        for (n in 1:nboot){
          #Sample the clusters with replacement
          clus.sam<-sample(clus.table$Cluster,replace=TRUE)
          data.boot<-NULL
          for(m in 1:numclus){
            data.cur<-data[data$cluster==clus.sam[m],]
            data.cur$cluster<-as.factor(m)
            data.boot<-rbind(data.boot,data.cur)
          }
          if(modeltype=='fixed'){
            reg.ps<-glm(pscore.formula,data=data.boot,family=binomial(link='logit'))
            data.boot$pscore<-predict(reg.ps,type='response')

          }else if(modeltype=='random'){
            reg.ps<-glmer(pscore.formula,data=data.boot,family=binomial(link='logit'),
                          nAGQ=0,
                          control=glmerControl(optimizer = "nloptwrap"))
            data.boot$pscore<-predict(reg.ps,type='response')
          }
          #cutoff
          datawt<-data.boot[data.boot$pscore>ps.cut[1] & data.boot$pscore<ps.cut[2],]

          data.temp<-datawt
          data.temp$Z<-0

          if(modeltype=='fixed'){
            #y.formula.boot<-update(y.formula,~. -cluster+clusterid)
            if(y.type=='binary'){
              ####Fixed Effects
              regDR.y<-glm(y.formula,data=datawt,family=binomial(link='logit'))
              y.0<-predict(regDR.y,newdata=data.temp,type='response')

            }else if (y.type=='continuous'){
              regDR.y<-lm(y.formula,data=datawt)
              y.0<-predict(regDR.y,newdata=data.temp)
            }
          }else if(modeltype=='random'){
            #y.formula.boot<-update(y.formula,~. -(1|cluster)+(1|clusterid))
            if(y.type=='binary'){
              regDR.y<-glmer(y.formula,
                             data=datawt,family=binomial(link='logit'),
                             nAGQ=0,
                             control=glmerControl(optimizer = "nloptwrap"))
              y.0<-predict(regDR.y,newdata=data.temp,type='response')
            } else if(y.type=='continuous'){
              regDR.y<-lmer(y.formula,data=datawt)
              y.0<-predict(regDR.y,newdata=data.temp)
            }
          }
          #Impute potential outcomes
          datawt$estDR.y0<-y.0

          #ATT Estimator
          estimate.boot[n]<-with(datawt,sum((Y - estDR.y0)*(Z-pscore)/(1-pscore))/sum(Z))

        }
        se<-sd(estimate.boot)

    } #End if se is required

  } #End if method = DR

  results<-list('estimate' = estimate,'se'=se)
  return(results)
  }



