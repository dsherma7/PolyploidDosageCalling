
#' Performs a Simulated Manual Call
#'
#' Following the pre-existing manual method of determining dosage from RQ values, this algorithm finds the obvious clusters among the RQ values and assigns each cluster a dosage.
#' @param Target The target gene for calling the dosage. Must follow the exact name in datasets.
#' @param AllData A serial data frame of all the relevant values including Job, Plate, RQ, CT, delta CT, and Sample Name.
GetManualCall <- function(AllData, Target = "Assay1 CN "){
  df <- NULL;
  Jobs <- unique(AllData$Job)
  AllData$Manual.Call <- NA;
  for (each.job in Jobs){
    Plates <- unique(AllData$Plate[which(AllData$Job == each.job)])
    for (each.plate in Plates){
      CurrentPlate <- intersect(intersect(which(AllData$Job == each.job),which(AllData$Plate == each.plate)),which(AllData$Target.Name == Target))
      RQ.Vals <- as.double(AllData$RQ[CurrentPlate]); 
      temp <- rep(NA,length(RQ.Vals));
      if (length(RQ.Vals) - sum(is.na(RQ.Vals)) > 3){
        temp[which(!is.na(RQ.Vals))] <- ManualCall.Trials(RQ.Vals)
      }else{
        print(paste("Warning:","Plate",each.plate,"Does Not Have Enough RQ Values to Perform Manual Call"))
      }
      AllData$Manual.Call[CurrentPlate][order(AllData$RQ[CurrentPlate])] <- temp;
      df <- c(df,temp);
    }
  }
  AllData <- AllData[which(AllData$Target.Name == Target),]
  
  return(AllData);
}

#' Runs Multiple Iterations of Manual Call
#'
#' This function runs multiple iterations of manual call with smaller and smaller refinement of clustering.
#' @param bin.threshold The minimum percentage that a bin (cluster) is allowed to be relative to all the values.
#' @param max The largest RQ value to be considered. NULL implies no limit.
#' @param min The smallest RQ value to be considered.
#' @param by The steps between from and to.
#' @param to The largest interval length for calling Manual Call with.
#' @param from The smallest interval length for calling Manual Call with.
#' @param ploidy The number of sets of chromosomes in a cell.
#' @param RQ.Vals A list of RQ values. Runs faster if this is already sorted.
ManualCall.Trials <- function(RQ.Vals, ploidy = 4, from = .01, to = 1, by = .01, min = 0, max = NULL, bin.threshold = .01){
  result <- result.broad <- result.narrow <- matrix(NA,nrow = length(seq(from,to,by)),ncol = sum(!is.na(RQ.Vals),na.rm = TRUE))
  
  range <- seq(from,to,by);
  for (i in 1:length(range)){
    result[i,] <- ManualCall(RQ.Vals, by.val = range[i], min = min, max = max, bin.threshold = bin.threshold,min.threshold = NULL)$Labels
    result.broad[i,] <- ManualCall(RQ.Vals, by.val = range[i], min = min, max = max, bin.threshold = bin.threshold, min.threshold = 4)$Labels
    result.narrow[i,] <- ManualCall(RQ.Vals, by.val = range[i], min = min, max = max, bin.threshold = bin.threshold, min.threshold = 1)$Labels
  }
  result <- rbind(result,result.broad,result.narrow)
  
  if (sum(apply(result,1,function(x){max(x,na.rm = TRUE) == ploidy}),na.rm = TRUE) > 0){
    result <- result[apply(result,1,function(x){max(x,na.rm = TRUE) == ploidy}),]
  }else if (sum(apply(result,1,function(x){max(x,na.rm = TRUE) %in% c(ploidy-1,ploidy,ploidy+1)}),na.rm = TRUE) > 0){
    result <- result[apply(result,1,function(x){max(x,na.rm = TRUE) %in% c(ploidy-1,ploidy,ploidy+1)}),]
  }else{
    result <- rep(NA,ncol(result))
  }
  
  if (!is.null(dim(result))){
    result <- apply(result,2,function(x){round(mean(x,na.rm = TRUE))})
  }
  return(result)
}

#' Simulated Manual Dosage Call
#'
#' Simulates pre-existing manual calling method by clustering sorted RQ values.
#' @param min.threshold Minimum number of samples allowed in each cluster.
#' @param bin.threshold The minimum percentage that a bin (cluster) is allowed to be relative to all the values.
#' @param max The largest RQ value to be considered. NULL implies no limit.
#' @param min The smallest RQ value to be considered.
#' @param by.val The steps between from and to.
#' @param RQ.Vals A list of RQ values. Runs faster if this is already sorted.
ManualCall <- function(RQ.Vals, by.val = .1, min = 0, max = NULL, bin.threshold = .01, min.threshold = 1){
  if (is.null(max)){max = max(RQ.Vals,na.rm = TRUE)+2*by.val}
  
  counts <- CountSections(RQ.Vals,by.val,min,max);
  index <- zero.index <- 1;
  splits <- zeros <- c(1,rep(NA,length(counts$Counts)-1));
  for (i in 2:length(counts$Counts)){
    if (counts$Counts[i] == 0 && counts$Counts[i-1] != 0){
      if (sum(counts$Counts[splits == index],na.rm = TRUE)/length(RQ.Vals) > bin.threshold){
        if (!is.null(min.threshold)){
          #  Labels[which(Labels == index)] <- index-1;
        }else{
          index <- index + 1;
        }
      }else{
        counts$Counts[i] <- -1; # So we will split at next 0
      }
    }
    if (counts$Counts[i] != 0 && counts$Counts[i-1] == 0){
      zero.index <- zero.index + 1;
    }
    
    splits[i] <- index;
    if (counts$Counts[i] == 0){
      zeros[i] <- zero.index;
    }else{
      zeros[i] <- NA;
    }
  }
  counts$Counts[which(counts$Counts == -1)] <- 0
  Zeros <- data.frame(Vals = unlist(lapply(seq(max(zeros,na.rm = TRUE)),function(x){sum(zeros == x,na.rm = TRUE)})))
  Chunks <- split(counts$Counts,splits)[1:(length(split(counts$Counts,splits))-1)]
  names(Chunks) <- 0:(length(Chunks)-1)
  Labels <- c(unlist(lapply(seq(length(Chunks)),function(x){rep(x-1,sum(Chunks[[x]],na.rm = TRUE))})),rep(length(Chunks),max(0,length(sort(RQ.Vals))-length(unlist(lapply(seq(length(Chunks)),function(x){rep(x-1,sum(Chunks[[x]],na.rm = TRUE))}))))))
  
  if (!is.null(min.threshold)){
    for (j in 1:max(Labels, na.rm = TRUE)){
      for (i in 1:max(Labels,na.rm = TRUE)){
        if (length(Labels[which(Labels == i)]) < min.threshold){
          Labels[which(Labels == i)] <- i-1;
        }
      }
    }
  }
  
  return(list(Chunks = Chunks,Labels = Labels, Zeros = Zeros))
}

#' Helper Function for Manual Call
#'
#' This function counts the the number of RQ values in a given range.
#' @param max The largest RQ value to be considered. NULL implies no limit.
#' @param min The smallest RQ value to be considered.
#' @param by.val The steps between from and to.
#' @param RQ.Vals A list of RQ values. Runs faster if this is already sorted.
CountSections <- function(RQ.Vals,by.val = .1,min = 0, max = 5){
  # Return RQ Values in a Given Interval
  return(data.frame(Counts = unlist(lapply(seq(min,max,by=by.val),
                                           function(x){
                                             length(RQ.Vals[intersect(which(RQ.Vals > x),which(RQ.Vals < x + by.val))])
                                           })),
                    Vals = seq(min,max,by.val)))
}
