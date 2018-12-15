
#' Semi-Supervised Clustering Algorithm
#'
#' Following the standard K-Means algorithm, we start by assigning each sample a random cluster value. Howewver, the known control values are assigned their own dosage. Then this algorithm performs consistent iterations of computing the mean of each cluster and then relabeling each sample to it's closest cluster. 
#' @param get.Prob TRUE/FALSE value indicating if probabilities are desired alongside dosage calls. Defaults to TRUE, but use FALSE for more efficiency if probability is not needed.f
#' @param outlier.threshold Used to determine the probability values. Samples outside with standard deviation larger than this value will have very low probability and will be thrown out in determining the remaining samples' probability.
#' @param converge.threshold If the percent change between two iterations is within this value then the algorithm will stop. Makes this more effecient.
#' @param max.iter Max number of iterations allowed to run for this algorithm. If set too low a warning will occur.
#' @param Labels The initial labels of all the data. Unknown samples are marked NA, while the controls are dosage +1 (since 0 is invalid). Leave all as NA for standard K-Means.
#' @param nClust Expected total number of clusters (i.e the possible dosage values). Default is 5 for tetraploids.
#' @param data Dataset used for clustering. Must be an \emph{n} x \emph{p} numeric matrix; where \emph{n} is the number of samples and \emph{p} is the number of features.
#' @keywords SemiSupervised, Clustering, K-Means, Cluster, Supervised, Semi
#' @export
MyClustering <- function(data, nClust = 5, Labels = rep(NA,nrow(data)), max.iter = 1000, converge.threshold = 0.05, outlier.threshold = 2, get.Prob = TRUE){
  
  if(class(data) != "data.frame"){
    warning(paste0('Warning data is of the class ',class(data),', it should be a data.frame'))
    return(NULL)
  }
  # K-Means does not work with missing data, so these rows are removed
  X <- data[which(apply(data,1,function(x){sum(is.na(x),na.rm = TRUE) == 0})),]
  badSamples <- NULL
  for (each.col in seq(ncol(X))){
    badSamples <- union(badSamples,which(is.na(X[,each.col])))
  }
  
  # Normalizing the data first improves convergence
  x <- apply(X,2,function(x){x/mean(x,na.rm = TRUE)})
  labels <- as.double(Labels)[which(apply(data,1,function(x){sum(is.na(x),na.rm = TRUE) == 0}))] + 1;
  d <- ncol(x) # of Features
  n <- nrow(x) # of Samples
  
  # Start by computing the center of each cluster
  dat <- data.frame(x, temp.labels = rep(0,nrow(x)))
  ClustMeans <- data.frame(matrix(rep(as.double(NA),nClust*ncol(x)),ncol = ncol(x)));
  names(ClustMeans) <- names(dat)[1:d]
  rownames(ClustMeans) <- as.integer(1:nClust)
  Clusters <- list()
  
  for (i in 1:nClust){          # Initialize the List of each cluster
    Clusters[[length(Clusters)+1]] <- matrix(as.double(NA),ncol = d,nrow = n)
  }
  
  if (n != length(labels)){    # Did not provide enough labels for the data
    warning("Labels Not the Same Size as non-NA Data")
  }
  
  # Set any known values to their assigned classes
  for (i in 1:length(labels)){ # Find Starting Labels
    if (!is.na(labels[i])){
      for (j in 1:d){
        Clusters[[labels[i]]][i,j] <- x[i,j]
      }
      dat$temp.labels[i] <- labels[i];
    }else{
      choice <- sample(seq(nClust),1)
      for (j in 1:d){
        Clusters[[choice]][i,j] <- x[i,j]
      }
      dat$temp.labels[i] <- choice
    }  
  }
  
  Num.Changes <- 99999; iters <- 1;
  while ((Num.Changes/n) > converge.threshold && iters < max.iter){ # converge on iteration (bad) or threshold (good)
    Num.Changes <- 0
    iters <- iters + 1;
    for (each.row in seq(nrow(ClustMeans))){ # Re-compute means
      ClustMeans[each.row,] <- apply(dat[which(dat$temp.labels == each.row),seq(d)],2,mean) 
    }
    
    # Change non-control labels
    labels.to.change <- dat[which(is.na(labels)),]
    changed.labels <- apply(labels.to.change,1,function(x){ChooseClosestCluster(ClustMeans, x[1:d], nClust)})
    Num.Changes <- sum(dat$temp.labels[which(is.na(labels))] != changed.labels)
    dat$temp.labels[which(is.na(labels))] <- changed.labels
    #print(Num.Changes)
  }
  if (iters >= max.iter){
    warning(paste("Clustering stopped prematurely with",Num.Changes,"changes of clusters on the last iteration.\n   Consider increasing max.iter"))
  }
  
  # Get likelihood of calls as an F distribution
  dat$temp.labels <- dat$temp.labels - 1;
  if (get.Prob){
    probs <- GetConf(fitAllClust = list(df = dat, Means = ClustMeans),limit = outlier.threshold, index = 1:d, lbl.rng = min(dat$temp.labels,na.rm = TRUE):max(dat$temp.labels,na.rm = TRUE));
  }else{
    probs <- NULL
  }
  data$Labels <- NA;
  data$Labels[which(apply(data[,1:d],1,function(x){sum(is.na(x),na.rm = TRUE) == 0}))] <- dat$temp.labels;
  
  return(list(df = dat, Means = ClustMeans, Probs = probs))
}


#' Determine Closest Cluster to Sample
#'
#' Helper function for MyClustering which chooses the closest cluster for a given sample
#' @param nClust Expected total number of clusters (i.e the possible dosage values). Default is 5 for tetraploids.
#' @param val The specific sample for which the closest cluster is desired.
#' @param ClustMeans An \emph{nClust} x \emph{p} numeric Matrix of the mean values of each paramater for each cluster.
ChooseClosestCluster <- function(ClustMeans, val, nClust){
  scores <- rep(0,nClust);
  for (i in 1:nClust){
    for (j in 1:length(val)){
      scores[i] <- scores[i] + (val[j] - ClustMeans[i,j])^2
    }    
  }
  retVal <- which(scores == min(scores))
  if (is.na(retVal) || length(retVal) == 0)
    retVal <- sample(seq(nClust),1)
  else if (length(retVal) > 1)
    retVal <- sample(seq(retVal),1)
  return(retVal)
}


#' Measures the Confidence in Prediction
#'
#' Used to measure the confidence of each dosage call following standard multivariate normal confidence intervals.
#' @param index Which parameters to be used for measuring confidence. Default assumes only 2 features, but often this number is 1:\emph{p}.
#' @param limit The number of standard deviations from the mean that would warrant being called an outlier.
#' @param fitAllClust A 2 dimensional list containing a data frame (named \emph{df}), with the paramater values for each sample along with the predicted cluster, and  the ClustMeans matrix, an \emph{nClust} x \emph{p} numeric Matrix of the mean values of each paramater for each cluster.
GetConf <- function(fitAllClust, limit, index = 1:2, lbl.rng = 0:2){
  clust.distances  <- lapply(lbl.rng,function(z){
    apply(fitAllClust$df[which(fitAllClust$df$temp.labels==z),index],1,function(x){
      sapply(index,function(q){
        distance(x[q],fitAllClust$Means[z+1,q]);
      })
    })
  })
  clust.variances  <- sapply(lbl.rng,function(z){
    apply(fitAllClust$df[which(fitAllClust$df$temp.labels==z),index],2,function(x){
      var(x)
    })
  })
  clust.deviation <- lapply(lbl.rng+1,function(z){
    if(sapply(clust.distances,function(x){!is.null(dim(x))})[z]){
      apply(clust.distances[[z]],2,function(x){
        x/sqrt(clust.variances[,z])
      })
    }else{
      NA
    }
  })
  Confidence <- lapply(clust.deviation,function(q){
    apply(q, 2, function(x){
      mean(sapply(x, function(z){
        if (is.na(z)){
          NA
        }else
          if(z >= 0){
            (1-pnorm(z,0,1))/.5;
          }else{
            pnorm(z,0,1)/.5
          }
      }),na.rm = TRUE)
    })
  })
  
  clust.deviation <- lapply(clust.deviation,as.data.frame)
  for (i in index){
    names(clust.deviation[[i]]) <- rownames(fitAllClust$df)[which(fitAllClust$df$temp.labels==i-1)]
  }
  badChoices <- sort(unlist(lapply(lapply(clust.deviation,function(x){
    x[,union(which(x[1,] > limit),which(x[2,] > limit))]}),names)))
  badChoices <- badChoices[!duplicated(badChoices)]
  Data <- fitAllClust$df;
  Data$Prob <- NA;
  Probs <- Confidence
  prob <- unlist(sapply(lbl.rng,function(each.dose){
    Data$Prob[which(Data$temp.labels == each.dose)] <- Probs[[each.dose+1]];
  }))
  Data$Prob[order(Data$temp.labels, decreasing = FALSE)] <- prob
  return(list(Data = Data, Deviation = clust.deviation, Variance =  clust.variances, Confidence = Confidence, Worst.Samples = badChoices))
}

#' Distance Measure
#'
#' The Euclidean metric for measuring the \emph{p} dimensional distance between samples.
#' @param y A numeric vector of \emph{p} numbers.
#' @param x A numeric vector of \emph{p} numbers.
distance <- function(x,y){
  tot <- 0
  for (each.feature in seq(length(x))){
    tot <- tot + (x[each.feature] - y[each.feature])^2
  }
  return(sqrt(tot))
}


