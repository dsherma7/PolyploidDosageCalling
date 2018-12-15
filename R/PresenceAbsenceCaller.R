#' Runs CallPA() on all files in a directory
#'
#' This function determined the Presence/Absence of the appropriate assay based on the Quant Studio output file.
#' @param PAData A list of lists of each plate's Quant Studio Results file for Presence/Absence. Must contain Allele1.Delta.Rn and Allele2.Delta.Rn
#' @param sheetname Name of the sheet containing Allele Delta RN values
#' @param Controls File containing the known presence/absence information. In general this contains one column for each type of assay. 
#' @param output.dir The directory to write the output files to, Plots will be placed in Plots sub directory
#' @param outlier.threshold The number of standard deviations away from the mean a sample must be to be considered an outlier
#' @keywords Presence-Absence, Dosage Calling, Assay
CallAllPA <- function(PAData,sheetname,Controls,output.dir,outlier.threshold=3,
                      header = c('Job','#Samples','Plate','Rack','RackOrder','DNAWell','Well','Well.Position',
                                 'FBRecId','RecId','TissID','TISSID2','Sample.Name','PlantID','NewPlantID','Clone-plant_ID','ClonalID','NewCloneID',
                                 'Entry','Name','Variety','CrossID','Cross#','Pedigree','PedigreeLong',
                                 'SNP.Assay.Name','Allele1.Delta.Rn','Allele2.Delta.Rn','Quality(%)','Call','Control',
                                 'Assay1','Assay1_Pr','Assay1_QC','Assay3','Assay3_Pr','Assay3_QC','Assay2','Assay2_Pr','Assay2_QC'),console=FALSE){
  library(plyr)
  df.all <- NULL
  jobs <- names(PAData)
  assays <- c('Assay3','Assay2','Assay1')
  for (j in seq(length(PAData))){   #j=1; i=1
    each.job <- PAData[[j]]
    plates = names(each.job)
    df.plates <- NULL
    for (i in seq(length(plates))){
      df.new <- NULL
      tryCatch({
        df.new = CallSinglePA(df=each.job[[i]],Plate=plates[i],sheetname=sheetname,Controls=Controls,
                        output.dir=output.dir,outlier.threshold=outlier.threshold,console=console)
      },error = function(e){
        print(paste0('Plate ',plates[i], ' returned the error: ',e))
      },finally = {
        if (is.null(df.new)){
          tryCatch({
            df.new <- ParsePAComments(each.job[[i]])
          },error = function(e){
            print(paste0('Plate ',plates[i], ' has a bad comments column. This is necessary for the API'))
          },finally = {
            df.new <- each.job[[i]]
          })
          assay = assays[sapply(assays,function(x)length(grep(x,df.new$SNP.Assay.Name[1]))>0)]
          df.new[,assay] <- 'Undetermined'
          df.new$Plate <- substr(plates[i],0,8)
          df.new[,paste0(assay,'_Pr')] <- df.new$Control <- NA
        }
        df.plates <- rbind.fill(df.plates,df.new)
      })
    }
    df.plates$Job <- jobs[j]
    df.all <- rbind.fill(df.all,df.plates)
  }
  df.all <- df.all[,header[which(header %in% names(df.all))]]
  names(df.all)[which(names(df.all)=='Sample.Name')] <- 'Sample'
  return(df.all)
}

#' Determines if a Given Assay is Present in a Sample
#'
#' This function call the semi-supervised clustering algorithm for each plate/job of a given dataset.
#' @param sheetname Name of the sheet containing Allele Delta RN values
#' @param Controls File containing the known presence/absence information. In general this contains one column for each type of assay. 
#' @param output.dir The directory to write the output files to, Plots will be placed in Plots sub directory
#' @param outlier.threshold The number of standard deviations away from the mean a sample must be to be considered an outlier
#' @keywords Presence-Absence, Dosage Calling, Assay
CallSinglePA <- function(df,Plate,sheetname,Controls,output.dir,outlier.threshold=2,console=FALSE){
  df$Allele1.Delta.Rn <- as.numeric(df$Allele1.Delta.Rn)
  df$Allele2.Delta.Rn <- as.numeric(df$Allele2.Delta.Rn)
  
  tryCatch({
    # Parse the Comments string to obtain the associated categorical values.
    df <- ParsePAComments(df)
  },error = function(e){
    print(paste0('Plate ',plates[i], ' has a bad comments column. This is necessary for the API'))
  })
  
  # Remove any NTC samples
  if(length(grep('NTC_',df$Sample.Name))>0)
    df <- df[-grep('NTC_',df$Sample.Name),]
  
  if (sum(!is.na(df$`Clone-plant_ID`)) < 10){
    stop(paste0('Plate ',Plate, ' had no plant IDs '))
  }
  
  # Match the necessary controls based on the SNP.Assay.Name of the raw data file
  if (length(grep('Assay3',df$SNP.Assay.Name[1]))>0){
    df$Control = Controls[match(df$`Clone-plant_ID`,Controls$Plant.ID),grep('Assay3',names(Controls))]
    Assay = 'Assay3'
  }else if (length(grep('Assay2',df$SNP.Assay.Name[1]))>0){
    df$Control = Controls[match(df$`Clone-plant_ID`,Controls$Plant.ID),grep('Assay2',names(Controls))]
    Assay = 'Assay2'
  }else if (length(grep('Assay1',df$SNP.Assay.Name[1]))>0){
    df$Control = Controls[match(df$`Clone-plant_ID`,Controls$Plant.ID),grep('Assay1',names(Controls))]
    Assay = 'Assay1'
  }else{
    stop('ERROR on SNP Assay Name')
  }
  
  if(sum(!is.na(df$Control),na.rm=T)==0){
    stop(paste0('Plate ',Plate, ' had no matching control names'))
  }
  
  badPlot = T
  tryCatch({
    # Try to find the most obviously called samples to increase the number of controls
    Key = data.frame(Level = c('Positive','Negative','NTC'),ID = as.numeric(c(0,1,2)),stringsAsFactors = FALSE)
    df$Control[intersect(grep('Homo',df$Call),which(df$Allele1.Delta.Rn > max(df$Allele1.Delta.Rn[which(df$Control=='Negative')])))] <- 'Positive'
    df$Control[grep('Hetero',df$Call)] <- 'Negative'
    lbls = Key$ID[match(df$Control,Key$Level)]
    good.values <- intersect(which(!is.na(df$Allele1.Delta.Rn)),which(!is.na(df$Allele2.Delta.Rn)))
    RemoveOutliers <- function(df,lbls,Key,outlier.threshold){
      for (cat in unique(Key$ID)){
        temp = data.frame(A1 = df$Allele1.Delta.Rn,A2 = df$Allele2.Delta.Rn,Call = lbls,stringsAsFactors = FALSE);
        are.outliers = isOutlier(temp[which(temp$Call == cat),1:2],threshold = outlier.threshold);
        lbls[which(lbls==cat)][are.outliers] = NA;
      }
      return(lbls)
    }
    # Run multiple iterations of SS Clustering while adding/removing to the control values    
    df$Control.Trimmed = RemoveOutliers(df,lbls,Key,outlier.threshold)
    lbls <- df$Control.Trimmed
    Clusts <- lapply(1:10,function(none){
      Clust = Quadrophenia::MyClustering(data.frame(cbind(df$Allele1.Delta.Rn,df$Allele2.Delta.Rn),stringsAsFactors = FALSE),
                           nClust = 3,Labels = as.numeric(lbls),get.Prob = TRUE)
      lbls <- RemoveOutliers(df[good.values,],Clust$Probs$Data$temp.labels,Key,outlier.threshold=0.5)
      lbls[which(!is.na(df$Control.Trimmed))] <- df$Control.Trimmed[which(!is.na(df$Control.Trimmed))]
      return(Clust)
    })
    
    # Revert back to original controls
    if (Assay == 'Assay3')
      df$Control = Controls[match(df$`Clone-plant_ID`,Controls$Plant.ID),grep('Assay3',names(Controls))]
    if (Assay  == 'Assay2')
      df$Control = Controls[match(df$`Clone-plant_ID`,Controls$Plant.ID),grep('Assay2',names(Controls))]
    if (Assay == 'Assay1')
      df$Control = Controls[match(df$`Clone-plant_ID`,Controls$Plant.ID),grep('Assay1',names(Controls))]
    
    Clust <- NULL
    max.acc <- 0
    for (i in seq(length(Clusts))){
      # Keep the most accurate run of the clustering algorithm with varied controls
      temp.lbl <- Key$Level[match(Clusts[[i]]$df$temp.labels,Key$ID)]
      eq <- df$Control == temp.lbl
      acc <- sum(eq,na.rm=T)/length(which(!is.na(eq)))
      pos <- suppressWarnings(min(Clusts[[i]]$df$df.Allele1.Delta.Rn[which(Clusts[[i]]$df$temp.labels == 0)],na.rm=T))
      neg <- suppressWarnings(max(Clusts[[i]]$df$df.Allele1.Delta.Rn[which(Clusts[[i]]$df$temp.labels == 1)],na.rm=T))
      diff <- pos-neg;
      
      if (acc > max.acc){
        Clust <- Clusts[[i]]
        max.acc <- acc
      }
    }
    
    df$Assay3 <- df$Assay3_Pr <- df$Assay3_QC <- df$Assay2 <- df$Assay2_Pr <- df$Assay2_QC <- df$Assay1 <- df$Assay1_Pr <- df$Assay1_QC <- NA
    df[good.values,Assay] <- Clust$df$temp.labels
    
    # Adjust weird calls
    max.pos.A1 <- mean(df$Allele1.Delta.Rn,na.rm=T)+2*sqrt(var(df$Allele1.Delta.Rn,na.rm=T))
    min.pos.A1 <- mean(df$Allele1.Delta.Rn[which(df$Control=='Negative')],na.rm=T)-3*sqrt(var(df$Allele1.Delta.Rn[which(df$Control=='Negative')],na.rm=T))
    min.pos.A2 <- min(df$Allele2.Delta.Rn[which(df$Control!='Positive')],na.rm=T)+0.25*max(df$Allele2.Delta.Rn[which(df$Control!='Positive')],na.rm=T)
    temp <- df[,Assay]
    temp[intersect(which(df$Allele1.Delta.Rn>max.pos.A1),which(df$Allele2.Delta.Rn>min.pos.A2))] <- 0
    temp[intersect(which(df$Allele1.Delta.Rn<max.pos.A1),which(df$Allele2.Delta.Rn>min.pos.A2))] <- 1
    temp[which(df$Allele2.Delta.Rn<min.pos.A2)] <- 2
    df[,Assay] <- temp
    
    means <- sapply(Key$ID,function(x){cbind(mean(df$Allele1.Delta.Rn[which(df[,Assay]==x)],na.rm=T),mean(df$Allele2.Delta.Rn[which(df[,Assay]==x)],na.rm=T)) })
    means[,Key$ID[match('NTC',Key$Level)]+1] <- sapply(means[,3],function(x)mean(c(x,0),na.rm=T))
    sds <- sapply(Key$ID,function(x){sqrt(cbind(var(df$Allele1.Delta.Rn[which(df[,Assay]==x)],na.rm=T),var(df$Allele2.Delta.Rn[which(df[,Assay]==x)],na.rm=T))) })
    
    for (i in Key$ID){
      # Predict the likelihood of the calls
      temp <- df[which(df[,Assay]==i),]
      mu = c(mean(temp$Allele1.Delta.Rn,na.rm=T),mean(temp$Allele2.Delta.Rn,na.rm=T))
      sd = sqrt(c(var(temp$Allele1.Delta.Rn,na.rm=T),var(temp$Allele2.Delta.Rn,na.rm=T)))
      p1 <- dnorm(temp$Allele1.Delta.Rn,mu[1],sd[1])/dnorm(mu[1],mu[1],sd[1])
      p2 <- dnorm(temp$Allele2.Delta.Rn,mu[2],sd[2])/dnorm(mu[2],mu[2],sd[2])
      df[which(df[,Assay] == i),paste0(Assay,'_Pr')] <- p1*p2
  }
    df[good.values,paste0(Assay,'_Pr')] <- apply(cbind(Clust$Probs$Data$Prob,df[good.values,paste0(Assay,'_Pr')]),1,function(x)mean(c(x),na.rm=T))
    badPlot = F
  },error=function(e){
    print(paste0('Plate ',Plate,' for assay ',Assay,' returned the error \n',as.character(e)))
  },finally={
    # Add QC column for downstream Firewalling and plot the results.
    df[,paste0(Assay,'_QC')] <- NA
    PlotPA(df=df,output.dir=output.dir,Plate=Plate,Assay=Assay,console=console,badPlot=badPlot)
    dev.flush();
    df[,Assay] <- Key$Level[match(df[,Assay],Key$ID)]
    df[grep('NTC',df[,Assay]),Assay] <- 'UD'
    df$Plate = Plate
  })
  return(df)
}

#' Determines if a Given Assay is Present in a Sample
#'
#' This function call the semi-supervised clustering algorithm for each plate/job of a given dataset.
#' @param sheetname Name of the sheet containing Allele Delta RN values
#' @param control.file Location of the file containing the known Presence Absence values for each assay
#' @param output.dir The directory to write the output files to, Plots will be placed in Plots sub directory
#' @param outlier.threshold The number of standard deviations away from the mean a sample must be to be considered an outlier
#' @keywords Presence-Absence, Dosage Calling, Assay
ParsePAComments <- function(df){
  # Parse comment string  
  
  if (!("Comments"%in%names(df)))
    return(df)
  
  comment.names <- function(x){
    temp <- strsplit(x,'&|=')[[1]][-c(1)] 
    return(temp[seq(1,length(temp),2)])
  }
  split.comments <- function(x){
    if (!is.na(x)){
      temp <- strsplit(x,'&|=')[[1]][-c(1)] 
      temp.names <- comment.names(x)
      temp.vals  <- temp[seq(2,length(temp),2)]
      temp.vals[which(temp.vals=='NA')] <- NA
      matched <- cbind(all.names,temp.vals[match(all.names,temp.names)])
      return(c(matrix(t(matched),nrow=1)))
    }
    return(rep(NA,comment.sz))
  }
  comment.sz <- max(sapply(df$Comments,function(x)length(strsplit(x,'&|=')[[1]][-c(1)])))
  full.names <- sapply(df$Comments,comment.names)
  all.names <- unique(unlist(full.names[1:length(full.names)]))
  
  comments <- t(sapply(df$Comments,split.comments))
  c.names <- comments[1,seq(1,ncol(comments),2)]
  comments <- data.frame(cbind(comments[,seq(2,ncol(comments),2)]),stringsAsFactors=F,row.names=NULL)
  names(comments) <- c.names
  df <- as.data.frame(cbind(df,comments))
  return(df)
}

#' Print Output for Presence/Absence Calling
#'
#' Helper Function for CallSinglePA
PlotPA <- function(df,output.dir="",Plate="",Assay="",console=FALSE,badPlot=F){
  library(scales);
  if(!console)
    png(filename=paste0(output.dir,'PA.',Plate,'.png'));
  Key = data.frame(Level = c('Positive','Negative','NTC'),ID = as.numeric(c(0,1,2)),
                   Center=c(80,78,88),col=c('darkred','darkgreen','darkblue'),PCH=c(8,8,8),stringsAsFactors = FALSE)
  if (badPlot){
    plot(1,1,cex=0,axes=F,xlab='',ylab='',main=paste0('Plate: ',Plate, '\nPresence/Absence for ',Assay))
    text(1,1,'Bad Data',cex=3)
  }else{
    
    badControls <- intersect(which(!is.na(df$Control)),which(is.na(df$Control.Trimmed)))
    df$pch <- Key$PCH[match(df$Control,Key$Level)]
    df$pch[is.na(df$pch)] <- 16
    #df$pch[badControls] <- 13
    df$col = as.double(df[,Assay])+2
    #df$col[badControls] <- 5
    df$col[!is.na(df$Control)] <- Key$col[match(df$Control[!is.na(df$Control)],Key$Level)]
    df$cex <- 2*df[,paste0(Assay,'_Pr')]
    df$cex[!is.na(df$Control)] <- 2
    # Add Cluster Centers
    centers <- as.data.frame(matrix(NA,ncol=ncol(df),nrow=length(Key$ID)))
    names(centers) <- names(df)
    for (each.center in seq(length(Key$ID))){
      k <- Key$ID[each.center]
      centers$Allele1.Delta.Rn[each.center] <- mean(df$Allele1.Delta.Rn[which(df[,Assay]==k)],na.rm=T)
      centers$Allele2.Delta.Rn[each.center] <- mean(df$Allele2.Delta.Rn[which(df[,Assay]==k)],na.rm=T)
      centers$cex[each.center] <- 2
      centers$pch[each.center] <- Key$Center[which(Key$ID==k)]
      centers$col[each.center] <- Key$col[which(Key$ID==k)]
    }
    df <- rbind(df,centers)
    
    tryCatch({
      xlim = c(0.4,max(df$Allele1.Delta.Rn,na.rm=T))
      ylim = c(0,max(df$Allele2.Delta.Rn,na.rm=T))
      
      plot(df$Allele1.Delta.Rn,
           df$Allele2.Delta.Rn,
           col = alpha(df$col,df[,paste0(Assay,'_Pr')]),
           main = paste0('Plate: ',Plate, '\nPresence/Absence for ',Assay),
           xlab = 'Allele1',ylab = 'Allele2',
           pch = df$pch,cex=df$cex,
           xlim=xlim, ylim=ylim);
      
      Accuracy = round(sum(Key$ID[match(df$Control,Key$Level)]==df[,Assay],na.rm = TRUE)/sum(!is.na(df$Control)),3)
      ymin = max(c(0,min(df$Allele2.Delta.Rn,na.rm=T)),na.rm = TRUE)
      xmin = (min(df$Allele1.Delta.Rn,na.rm=T)+max(df$Allele1.Delta.Rn,na.rm=T))/2
      text(xmin,ymin,paste0("% Correct (Controls): ",Accuracy))
      #text(2.2,ymin,paste0('Likelihood:'),cex = .8)
      #for (t in seq(0,1,.1))
      #  text(2.2+1.3*t,ymin,t,cex=0.6*(1.5*t),col=alpha(3,t)) 
        
      legend("bottomright",legend = c('Positive','Negative','UD','Known','Center'),col = c(2:4,1,1),lwd = 1,pch = c(16,16,16,8,65));
    },error = function(e){
      print(paste0('Plotting ',Plate,' returned the error:\n',e))
    })
  }
  if (!console)
    dev.off();
}

#' Determines if Controls are Outliers
#'
#' Helper Function for CallPA
isOutlier <- function(x,threshold = 2){
  outliers <- rep(FALSE,nrow(x))
  tryCatch({
    outliers <- mahalanobis(x,colMeans(x),cov(x))>threshold
  },error = function(e){
    print(paste0('mahalanobis returned the error: ',e))
  })
  return(outliers)
}
