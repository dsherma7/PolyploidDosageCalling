#' Main Dosage Calling Function for the API
#'
#' Uses multiple models to fit the data, then uses each of these models to determine both the best fit, and liklihood of that fit. The models used in this algorithm are 
#' \describe{
#' \item{AllComboFit}{Following the \eqn{2^-\Delta\Delta Ct} method, this model considers all possible combinations of controls in order to generate a calibration \eqn{\Delta Ct} value. }
#' \item{Semi-Supervised K-means}{Uses standard unsupervised K-means clustering; however, by fixing controls values to their respective clusters it takes advantage of known information.}
#' \item{Linear Model}{Again following the \eqn{2^-\Delta\Delta Ct}, this method uses the fact that \eqn{log2(Dosage) ~ \Delta\Delta Ct}, which simplifies into \eqn{\log2(Dosage) ~ a + b\Delta Ct}. This way there is no dependence on finding a proper calibration sample, and we just use regression to determine \emph{a} and \emph{b}.}
#' \item{Simulated Manual Call}{Aligns the sorted RQ values and finds the natural clusters in the RQ values. It then labels this clusters from 0 to 4.}
#' }
#' @param header All the relevant variables for the final output.
#' @param swap.file Location of file containing any sample changes during tissue sampling.  
#' @param Ref.Vars Generate a new column containing the reference's value for any Ref.Vars. Default is CT to generate CT.PEF.
#' @param Output.File Name of output file containing each model's dosage/probability as well as the overall dosage/probability.
#' @param RQ.Threshold If RQ goes above this value we know something went wrong with the data. Default = 7, since RQ should never be above 5.
#' @param assays List of possible assays. Default is  c('Assay3','Assay2','Assay1','Assay1#','Assay3#'), but this could simply match the target on the raw data file.
#' @param Output.Dir Name of the directory to store the files. The path should be relative to Directory's path.
#' @param Target.Column Column name that lists the target assay. For standard Quant Studio output, this is "Target.Name".
#' @param Console Should plots should be printed to console instead of output directory? (default = False)
#' @param Skip.List List of Jobs to skip when searching for unupdated jobs. Should be a .csv file with a single column of Job names.
#' @param sheet Name of the sheet in the Raw Data where RQ, CT, and delta CT values are stored. Default is "Results" for an unchanged raw data file.
#' @param Job Name of the specific Job to run now during this call. 
#' @param Directory Name of the directory storing the raw data for the jobs (Output from Quant Studio)
#' @param Control.File If specified, these controls will be used for each job
#' @param Control.Dir Directory containing the Control Files (Ignored, if Control.File is specified)
#' @keywords Tetraploid, Polyploid, Dosage Calling, PCR, qPCR, Real Time
#' @export
CallDosageAPI <- function(Control.Dir, Control.File=NULL, Directory = NULL, Job = NULL, sheet = "Results", Skip.List = NULL, Console=F,
                       Target.Column = "Target.Name",Target = "Assay1", Output.Dir = NULL, assays = c('Assay3','Assay2','Assay1','Assay1#','Assay3#'), RQ.Threshold = 7,
                       Output.File = paste0(Output.Dir,"/",c(Job,basename(getwd()))[1],".Dosage.csv"), REF.Vars = c("CT"),swap.file = paste0(Output.Dir,'/src/Swaps.csv'),
                       header = c('...')){

  # Skip jobs in the Skip.List file. Put non-standard or old jobs here
  if (!is.null(Skip.List))  
    Skip <- sapply(read.csv(Skip.List,stringsAsFactors = FALSE,header = FALSE,quote = ""),as.character)
  if (Job %in% Skip){
    print(paste("Job",Job,"Has Been Skipped Per",Skip.List))
    return(paste("Job",Job,"Has Been Skipped Per",Skip.List))
  }
  
  # Switch to the provided directory where the Raw Data files are located
  curr.dir <- getwd();
  if (!is.null(Directory))   
    setwd(Directory)
  Directory <- getwd();
  if (is.null(Output.Dir))
    Output.Dir <- getwd();
  
  # Determine if their are new PCR Plates in the job
  PCRPlates <- list.files(Job)
  NotDone <- PCRPlates; 
  Completed <- list.files(paste0(c(Output.Dir,Job),collapse="/"))[grep('.csv',list.files(paste0(c(Output.Dir,Job),collapse="/")))]
  NotDone <- NULL
  for (each.plate in PCRPlates){
    if (length(grep(gsub('.xls','',each.plate),paste0(Completed,collapse=',')))==0)
      NotDone <- c(NotDone,each.plate)
  }
  if (is.null(NotDone)){
    print(paste("No New PCR Plates Exported for Job",Job))
    setwd(curr.dir);
    return(paste("No New PCR Plates Exported for Job",Job));
  }
  
  # Read in any data associated with target Job, then split into PA & Dosage jobs
  Dirs <- lapply(list.files(),dir); 
  names(Dirs) <- list.files();
  if (!is.null(Job))
    Dirs <- Dirs[which(names(Dirs) == Job)]
  Dirs[[Job]] <- NotDone
  Sub.Dirs <- lapply(lapply(seq(length(Dirs)),function(x){unlist(gsub(paste0(".xls|.xlsx|.txt|.csv|",paste0(names(Dirs)[x],"-")),"",Dirs[[x]]))}),as.character);

  BothData <- LoadAllData(Dirs = Dirs, sheet = sheet, want.list = TRUE)

  PAData <- lapply(BothData,function(each.plate){
    each.plate[which(sapply(seq(length(each.plate)),function(i)is.null(each.plate[[i]]$RQ)))]
    })
  PAData <- PAData[sapply(PAData,function(x)length(x)>0)]
  
  AllData <- lapply(BothData,function(each.plate){
    each.plate[which(sapply(seq(length(each.plate)),function(i)is.null(each.plate[[i]]$Allele1.Delta.Rn)))]
  })    
  AllData <- AllData[sapply(AllData,function(x)length(x)>0)]
  
  print(paste0('Job ',Job,' has ', sum(unlist(sapply(AllData,length)),na.rm=TRUE),' Dosage Calling Plates'))
  print(paste0('Job ',Job,' has ', sum(unlist(sapply(PAData,length)),na.rm=TRUE),' Presence/Absence Plates'))
  
  Fitted <- Manual.Dosage <- SS.Clusters <- CombinedDosageCalls <- CombinedPACalls <- NULL
  ## Dosage Calling
  if (length(AllData)>0){ 
  # Determine the target assay and load the correct 
  # controls. Add to this for new Assays.
  Target <- AllData[[1]][[1]][1,Target.Column]
  if (is.null(Control.File)){
    if (length(grep('Assay1',Target))>0)
      Control.File <- paste0(Control.Dir,'ControlsAssay1CN.csv')
    if (length(grep('Assay3',Target))>0)
      Control.File <- paste0(Control.Dir,'ControlsAssay3CN.csv')
    if (length(grep('Assay2',Target))>0)
      Control.File <- paste0(Control.Dir,'ControlsAssay2CN.csv')
    }
    Controls <- read.csv(Control.File)
    if (is.null(Controls)){
      print(paste0("Error:  Job ",Job," has bad controls. Either the target assay (",Target,") is unknown or the file ",Control.File," does not exist!"))
      return(paste0("Error:  Job ",Job," has bad controls. Either the target assay (",Target,") is unknown or the file ",Control.File," does not exist!"))
    }
    
    # Parse the comment string for any number of comments
    # This is how information gets passed through Quant Studio 
    Sub.Dirs <- lapply(seq(length(AllData)),function(x){
      Sub.Dirs[[x]][which(Sub.Dirs[[x]] %in% names(AllData[[x]]))]
    })
    comment.names <- function(x){
      temp <- strsplit(x,'&|=')[[1]][-c(1)] 
      return(temp[seq(1,length(temp),2)])
    }
    split.comments <- function(x,all.names){
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
    ParseComments <- function(each.job){
      
      return(lapply(each.job,function(df){
        comment.sz <- max(sapply(df$Comments,function(x)length(strsplit(x,'&|=')[[1]][-c(1)])))
        full.names <- sapply(df$Comments,comment.names)
        all.names <- unique(unlist(unique(full.names[1:length(full.names)])))
        
        comments <- t(sapply(df$Comments,function(x)split.comments(x,all.names)))
        c.names <- comments[1,seq(1,ncol(comments),2)]
        comments <- data.frame(cbind(comments[,seq(2,ncol(comments),2)]),stringsAsFactors=F,row.names=NULL)
        names(comments) <- c.names
        
        df <- as.data.frame(cbind(df,comments))
        return(df)
      }))
    }
    AllData <- lapply(AllData,ParseComments)
    
    # Compute RQ values from scratch after removing 
    # bad control values. Also, helps minimize batch effects.
    outlier <- 2
    badRQs <- NULL
    MaxZs <- NULL
    for (j in seq(length(AllData))){
      for (i in seq(length(AllData[[j]]))){
        temp <- AllData[[j]][[i]]
        temp$Control <- Controls$Dosage[match(temp$`Clone-plant_ID`,Controls$Name)]
        Ref <- -Inf
        for (dose in seq(4,0)){
          dCTs <- as.numeric(temp$Delta.Ct[which(temp$Control==dose)])
          dCTs[which(dCTs<Ref)] <- NA
          Ref <- mean(dCTs,na.rm = TRUE)
          sd <- sqrt(var(dCTs,na.rm=TRUE))
          Zs <- abs((dCTs-Ref)/sd)
          maxZ <- max(Zs,na.rm=TRUE)
          dCTs[which(Zs > outlier)] <- NA
          temp$Delta.Ct[which(temp$Control==dose)] <- dCTs
        }
        MaxZs <- c(MaxZs,maxZ)
        dCTs <- as.numeric(temp$Delta.Ct[which(temp$Control==1)])
        temp$RQ <- 2^(mean(dCTs,na.rm=T)-as.numeric(temp$Delta.Ct))
        Ref <- Inf
        for (dose in seq(4,0)){
          RQs <- as.numeric(temp$RQ[which(temp$Control==dose)])
          RQs[which(RQs>Ref)] <- NA
          Ref <- min(RQs,na.rm=T)
          temp$RQ[which(temp$Control==dose)] <- RQs
        }
        temp$Plate <- names(AllData[[j]])[i]
        badRQs <- rbind.fill(badRQs,temp[which(temp$RQ > RQ.Threshold),])
        temp$RQ[which(temp$RQ > RQ.Threshold)] <- NA
        AllData[[j]][[i]] = temp
      }
    }
    
    ### Basic data reformatting
    # Fix any duplicate names in raw data and controls
    Dirs <- lapply(AllData,names);
    names(Dirs) <- names(AllData);
    if (is.null(AllData)){
      print(paste("No Matching Controls for Job",Job))
      setwd(curr.dir)
      return(paste("No Matching Controls for Job",Job))
    }
    NoDupsLst <- RemoveDuplicateNames(AllData, Dirs, Controls,Target.Name = Target);
    AllData <- NoDupsLst$df;
    Controls <- NoDupsLst$Controls; rm(NoDupsLst);
    CombinedData <- ConvertList(AllData);
    
    # Separate function for each of the dosage calling methods. 
    tryCatch({Fitted <- GetDosage(AllData=AllData, Controls=Controls, Target=Target)},
             error = function(cmd){print(paste("Linear Fit Resulted in an Error: ",cmd))});
    tryCatch({Manual.Dosage <- GetManualCall(CombinedData, Target = Target)},
             error = function(cmd){print(paste("Manual Call Resulted in an Error: ",cmd))});
    tryCatch({SS.Clusters <- GetSemiSuperClusters(AllData, Controls = Controls, Target = Target)},
             error = function(cmd){print(paste("Semi-Supervised Fit Resulted in an Error: ",cmd))});
    dir.create(paste0(c(Output.Dir,Job),collapse="/"))
    # Combine each method into a single dataframe with the dosage calls for each sample
    Combined            <- CombineDosages(Dosages=Fitted, Manual=Manual.Dosage, SemiSuper=SS.Clusters,
                                          Target=Target, Jobs=names(AllData), Plates=Sub.Dirs)
    Combined <- Combined[,which(!(names(Combined)%in%assays))]
    
    # Determine the final dosage calls, create the plots, and swap any samples defined in the swap.file.
    CombinedDosageCalls <- GetPredicted(df=Combined, Controls = Controls, Output.Dir=paste0(c(Output.Dir,Job),collapse="/"), Target=Target, Console=Console, LL.min=0.0)
    CombinedDosageCalls <- CombinedDosageCalls[,header[which(header%in%names(CombinedDosageCalls))]]
    CombinedDosageCalls <- MakeSwaps(CombinedDosageCalls,swap.dir=swap.file,Job)
    
    # Create a separate output file for each PCR plate.
    Plates <- unique(CombinedDosageCalls$Plate)
    ByPlate <- lapply(Plates,function(x){CombinedDosageCalls[which(CombinedDosageCalls$Plate == x),]})
    names(ByPlate) <- Plates
    for (each.plate in Plates){
      print(paste0(each.plate,': ',!file.exists(paste0(c(Output.Dir,Job,paste0('Dosage.',each.plate,'.csv')),collapse="/"))))
      if (!file.exists(paste0(c(Output.Dir,Job,paste0('Dosage.',each.plate,'.csv')),collapse="/")))
        write.csv(x=ByPlate[[each.plate]],file=paste0(c(Output.Dir,Job,paste0('Dosage.',each.plate,'.csv')),collapse="/"),row.names=F,quote=F,na='');
      
    }
    dev.flush() # In case any plots failed
  }# Perform Copy Number Dosage Calling
  
  ## Presence/Absence
  if (length(PAData)>0){
    dir.create(paste0(c(Output.Dir,Job),collapse="/"))
    dir.create(paste0(c(Output.Dir,Job,'Plots/'),collapse="/"))
    
    # Read the controls and determine presence/absence for each of the plates in PAData. Swap samples as needed by swap.file
    ControlPAFile <- read.csv(file=paste0(Control.Dir,'ControlsPA.csv'), header = TRUE, stringsAsFactors = FALSE)
    CombinedPACalls <- CallAllPA(PAData=PAData, sheetname=sheet, Controls=ControlPAFile,
                                output.dir=paste0(c(Output.Dir,Job,'Plots/'),collapse="/"), outlier.threshold=2, console=Console)
    CombinedPACalls <- CombinedPACalls[,header[which(header%in%names(CombinedPACalls))]]
    CombinedPACalls <- MakeSwaps(CombinedPACalls,paste0(c(Output.Dir,'src','Swaps.csv'),collapse="/"),Job)
    
    # Create a separate output file for each PCR plate.
    Plates <- unique(CombinedPACalls$Plate)
    ByPlate <- lapply(Plates,function(x){CombinedPACalls[which(CombinedPACalls$Plate == x),]})
    names(ByPlate) <- Plates
    for (each.plate in Plates){
      if (!file.exists(paste0(c(Output.Dir,Job,paste0('PA.',each.plate,'.csv')),collapse="/")))
        write.csv(ByPlate[[each.plate]],paste0(c(Output.Dir,Job,paste0('PA.',each.plate,'.csv')),collapse="/"),row.names = FALSE,quote = FALSE,na = '');
    }
    dev.flush() # In case any plots failed
  }# Perform Presence/Absence Calling
  
  # Append the new calls to the Summary file
  Summary=BuildSummary(Job=Job,NotDone=NotDone)
  setwd(curr.dir);
  return(list(DosageCalls = CombinedDosageCalls,PACalls = CombinedPACalls,Summary=Summary))  
}


#' Calls dosage for a given input file
#'
#' Uses multiple models to fit the data, then uses each of these models to determine both the best fit, and liklihood of that fit.
#' 
#' @param REF.Vars REF.vars will be included in the output file with the reference assay values. (Default = CT)
#' @param Output.File The name of the output file containing the final dosage calls. Defaults to Dosage.csv.
#' @param RQ.Threshold If RQ goes above this value we know something went wrong with the data. (Default = 7), since RQ should never be above 5.
#' @param Target.Column Column containing the Target Assays. Default="Target.Name"
#' @param Console Should plots be printed to the console? Default is False which prints the plots to Output.Dir.
#' @param sheet The name in File where Clone-Plant_ID and RQ values are stored. Default = "Results".
#' @param Directory The directory containing the necessary files. Default to NULL which will use current directory.
#' @param Control.File Name of the file containing the controls. Must contain 2 columns, Name and Dosage, with Name matching the values in Clone-Plant_ID for control samples. (Default = Controls.csv)
#' @param File Name of the file containing the raw data from a Quant Studio qPCR run. The column Clone-plant_ID must be added.
#' @return A dataset containing the dosage calls of each of the methods and a final determination of the dosage and likelihood.
#' @seealso \code{\link{?CallPA}} for calling the presence/absence of provided data
#' @export
#' @examples
#' # Obtain the appropriate Control and Raw Data files. You can find these in the 
#' # data-raw/ folder of the Quadrophenia package files.They should match the data below.
#' 
#' # A 41 line header
#' # Followed by the output of QuantStudio's Results sheet
#' # Note that Clone-plant_ID must be added manually to this file
#' 
#' require(Quadrophenia)
#' RawData.Assay1CN[1:41,1] # Header
#' RawData.Assay1CN[42,]    # Columns of data
#' 
#' # The controls file has two columns: Name and Dosage. Name must 
#' # match the Clone-plant_ID column for the controls in the raw data file.
#' 
#' head(ControlsAssay1CN)
#' 
#' # Then run CallDosage() with these files as inputs
#' 
#' output = CallDosage(File='RawData.Assay1CN.xls',Control.File='ControlsAssay1CN.csv')
#' 
#' # The output directory Output/ should contain 1 .csv file and a Plots/ directory.
#' # The Plots/ directory contains two plots, one standard for manual calling, and 1
#' # plot showing the distribution of each class. The output file has columns for the
#' # Dosage calls as well as a likelihood for each call.
#' 
#' head(output)
#' 
#' ## Other ways to run CallDosage()
#' # You could also specifiy an output file or a directory to look for the raw data files
#' 
#' dir.create('Files/')
#' file.copy(c('RawData.Assay1CN.xls','ControlsAssay1CN.csv'),'Files/')
#' output = CallDosage(File='RawData.Assay1CN.xls',
#'                     Control.File='ControlsAssay1CN.csv',
#'                     Directory='Files/',Output.Dir='NewOutput')
#' 
#' # Or specify the files in Files/ wihout specifying the directory. You can
#' # also change the name of the output file to something more informative.
#' 
#' output = CallDosage(File='Files/RawData.Assay1CN.xls',
#'                     Control.File='Files/ControlsAssay1CN.csv',
#'                     Output.Dir='ThirdOutput',
#'                     Output.File='ThirdOutputDosage.csv')
#' 
#' # Duplicates will automatically be renamed to avoid conflicts
#' 
#' output = CallDosage(File='Files/RawData.Assay1CN.xls',
#'                     Control.File='Files/ControlsAssay1CN.csv',
#'                     Output.Dir='ThirdOutput',
#'                     Output.File='ThirdOutputDosage.csv')
#' 
#' # Clean up the test files
#' unlink(c('Files','Output','ThirdOutput'),recursive=T)
CallDosage <- function(File, Control.File="Controls.csv", Output.Dir="Output", Directory = NULL,sheet ="Results", Console=F,
                       Target.Column="Target.Name",RQ.Threshold=7,Output.File = "Dosage.csv", REF.Vars = c("CT")){
     
  # Switch to directory containing the files                     
  curr.dir <- getwd();
  if (!is.null(Directory))   
    setwd(Directory)
  Directory <- getwd();

  if(!is.null(Output.Dir)){
    Output.File <- paste0(Output.Dir,'/',Output.File)
    dir.create(Output.Dir)
  }
  
  # Fix name to avoid any duplicate name conflicts
  ctr = 1
  while(file.exists(Output.File)){
    Output.File <- paste0(gsub("-Copy|([0-9])*.csv$","",Output.File),"-Copy",ctr,".csv") 
    ctr = ctr+1;
  }
  Job <- strsplit(gsub('.xls|.xlsx|.csv','',File),'/')[[1]]
  Job <- Job[length(Job)]
  if (ctr > 1)
    Job <- paste0(Job[length(Job)],'-Copy',ctr-1)
  print(Output.File)
  NotDone = File;
  
  if (is.null(NotDone) | !file.exists(NotDone)){
    warning("Invalid raw data file!")
    setwd(curr.dir);
    return("Invalid raw data file!");
  }

  Dirs <- list(NotDone)
  names(Dirs) <- getwd()
  Sub.Dirs <- list(Dosage=Job)
  
  Fitted <- Manual.Dosage <- SS.Clusters <- NULL;
  
  BothData <- LoadAllData(Dirs = Dirs, sheet = sheet, want.list = TRUE)
  
  AllData <- lapply(BothData,function(each.plate){
    each.plate[which(sapply(seq(length(each.plate)),function(i)is.null(each.plate[[i]]$Allele1.Delta.Rn)))]
  })    
  AllData <- AllData[sapply(AllData,function(x)length(x)>0)]
  names(AllData[[1]]) <- Job
  if (sum(unlist(sapply(AllData,length)),na.rm=TRUE) <= 0){
    warning("Bad input file! Contains Allele1 Delta Rn, which is for Presence/Absence jobs.")
    setwd(curr.dir)
    return("Bad input file! Contains Allele1 Delta Rn, which is for Presence/Absence jobs.")
  }
  
  ## Dosage Calling
  Target <- AllData[[1]][[1]][1,Target.Column]
  
  if (!file.exists(Control.File)){
    warning("Bad Control File!")
    setwd(curr.dir)
    return("Bad Control File!")
  }
  Controls <- read.csv(Control.File)
  
  ### Compute RQ Values 
  outlier <- 2
  badRQs <- NULL
  MaxZs <- NULL
  for (j in seq(length(AllData))){
    for (i in seq(length(AllData[[j]]))){
      temp <- AllData[[j]][[i]]
      control.vals <- matrix(sapply(Controls$Dosage,function(x)rep(x,2)));
      temp$Control <- Controls$Dosage[match(temp$`Clone-plant_ID`,Controls$Name)]
      Ref <- -Inf
      for (dose in seq(4,0)){
        dCTs <- as.numeric(temp$Delta.Ct[which(temp$Control==dose)])
        dCTs[which(dCTs<Ref)] <- NA
        Ref <- mean(dCTs,na.rm = TRUE)
        sd <- sqrt(var(dCTs,na.rm=TRUE))
        Zs <- abs((dCTs-Ref)/sd)
        maxZ <- suppressWarnings(max(Zs,na.rm=TRUE))
        dCTs[which(Zs > outlier)] <- NA
        temp$Delta.Ct[which(temp$Control==dose)] <- dCTs
      }
      MaxZs <- c(MaxZs,maxZ)
      dCTs <- as.numeric(temp$Delta.Ct[which(temp$Control==1)])
      temp$RQ <- 2^(mean(dCTs,na.rm=T)-as.numeric(temp$Delta.Ct))
      Ref <- Inf
      for (dose in seq(4,0)){
        RQs <- as.numeric(temp$RQ[which(temp$Control==dose)])
        RQs[which(RQs>Ref)] <- NA
        Ref <- suppressWarnings(min(RQs,na.rm=T))
        temp$RQ[which(temp$Control==dose)] <- RQs
      }
      temp$Plate <- names(AllData[[j]])[i]
      badRQs <- rbind.fill(badRQs,temp[which(temp$RQ > RQ.Threshold),])
      temp$RQ[which(temp$RQ > RQ.Threshold)] <- NA
      AllData[[j]][[i]] = temp
    }
  }
  
  ### Basic data reformatting
  Dirs <- lapply(AllData,names);
  names(Dirs) <- names(AllData);
  if (is.null(AllData)){
    print(paste("No Matching Controls for Job",Job))
    setwd(curr.dir)
    return(paste("No Matching Controls for Job",Job))
  }
  NoDupsLst <- RemoveDuplicateNames(AllData, Dirs, Controls,Target.Name = Target);
  AllData <- NoDupsLst$df;
  Controls <- NoDupsLst$Controls; rm(NoDupsLst);
  
  CombinedData <- ConvertList(AllData);
  tryCatch({Fitted <- GetDosage(AllData=AllData, Controls=Controls, Target=Target);},
           error = function(cmd){print(paste("Linear Fit Resulted in an Error: ",cmd))});
  tryCatch({Manual.Dosage <- GetManualCall(CombinedData, Target = Target)},
           error = function(cmd){print(paste("Manual Call Resulted in an Error: ",cmd))});
  tryCatch({SS.Clusters <- GetSemiSuperClusters(AllData, Controls = Controls, Target = Target)},
           error = function(cmd){print(paste("Semi-Supervised Fit Resulted in an Error: ",cmd))});
  
  Combined            <- CombineDosages(Dosages=Fitted, Manual=Manual.Dosage, SemiSuper=SS.Clusters,
                                        Target=Target, Jobs=names(AllData), Plates=Sub.Dirs)
  
  CombinedDosageCalls <- GetPredicted(df=Combined, Controls = Controls, Output.Dir=Output.Dir, Target=Target, Console=Console, LL.min=0.0)
  write.csv(x=CombinedDosageCalls,file=Output.File,row.names=F,quote=F,na="")
  
  
  setwd(curr.dir)
  return(CombinedDosageCalls)  
}



#' Calls Presence/Absence for a given input file
#'
#' Uses a clustering algorithm to cluster the 3 categories: Positive, Negative, NTC on the file provided.
#' Requires known control values as specified in the Control.File.
#' @param Output.File The name of the output file containing the final dosage calls. Defaults to Dosage.csv.
#' @param Console Should plots be printed to the console? Default is False which prints the plots to Output.Dir.
#' @param sheet The name in File where Clone-Plant_ID and RQ values are stored. Default = "Results".
#' @param Directory The directory containing the necessary files. Default to NULL which will use current directory.
#' @param Control.File Name of the file containing the controls. Must contain 2 columns, Name and Presence, with Name matching the values in Clone-Plant_ID for control samples. (Default = Controls.csv)
#' @param File Name of the file containing the raw data from a Quant Studio qPCR run. The column Clone-plant_ID must be added.
#' @return A dataset containing the dosage calls of each of the methods and a final determination of the dosage and likelihood.
#' @seealso \code{\link{CallDosage}} for calling the dosage of provided data
#' @export
#' @examples
#' # Obtain the appropriate Control and Raw Data files. You can find these in the
#' # data-raw/ folder of the Quadrophenia package files.They should match the data below.
#' 
#' # A 41 line header
#' # Followed by the output of QuantStudio's Results sheet
#' # Note that Clone-plant_ID must be added manually to this file
#' 
#' require(Quadrophenia)
#' RawData.Assay2[1:41,1] # Header
#' RawData.Assay2[42,]    # Columns of data
#' 
#' # The controls file has two columns: Clone-plant_ID and Presence. Clone-plant_ID must
#' # match the Clone-plant_ID columns for the associated controls in the raw data file.
#' 
#' head(ControlsAssay2)
#' 
#' # Then run CallPA() with these files as inputs
#' 
#' output = CallPA(File='RawData.Assay2.xls',Control.File='ControlsAssay2.csv')
#' 
#' # The output directory Output/ should contain 1 .csv file and a Plots/ directory.
#' # The Plots/ directory contains one plot showing the clustering of Positive, Negative,
#' # and Undetermined calls for the expression of the target & reference assays.
#' 
#' head(output)
#' 
#' ## Other ways to run CallPA()
#' # You could also specifiy an output file or a directory to look for the raw data files
#' 
#' dir.create('Files/')
#' file.copy(c('RawData.Assay2.xls','ControlsAssay2.csv'),'Files/')
#' output = CallPA(File='RawData.Assay2.xls',
#'                 Control.File='ControlsAssay2.csv',
#'                 Directory='Files/',Output.Dir='NewOutput')
#' 
#' # Or specify the files in Files/ wihout specifying the directory. You can
#' # also change the name of the output file to something more informative.
#' 
#' output = CallPA(File='Files/RawData.Assay2.xls',
#'                 Control.File='Files/ControlsAssay2.csv',
#'                 Output.Dir='ThirdOutput',
#'                 Output.File='ThirdOutputPA.csv')
#' 
#' # Duplicates will automatically be renamed to avoid conflicts
#' 
#' output = CallPA(File='Files/RawData.Assay2.xls',
#'                 Control.File='Files/ControlsAssay2.csv',
#'                 Output.Dir='ThirdOutput',
#'                 Output.File='ThirdOutputPA.csv')
#' 
#' # Clean up the test files
#' unlink(c('Files','Output','ThirdOutput'),recursive=T)
CallPA   <- function(File, Control.File="Controls.csv",Output.Dir="Output",Directory=NULL,sheet="Results",Console=F,Output.File="PA.csv"){
  
  curr.dir <- getwd();
  if (!is.null(Directory))   
    setwd(Directory)
  
  Directory <- getwd();
  
  if(!is.null(Output.Dir)){
    Output.File <- paste0(Output.Dir,'/',Output.File)
    dir.create(Output.Dir)
  }
  
  ctr = 1
  while(file.exists(Output.File)){
    Output.File <- paste0(gsub("-Copy|([0-9])*.csv$","",Output.File),"-Copy",ctr,".csv") 
    ctr = ctr+1;
  }
  Job <- strsplit(gsub('.xls|.xlsx|.csv','',File),'/')[[1]]
  Job <- Job[length(Job)]
  if (ctr > 1)
    Job <- paste0(Job,'-Copy',ctr-1)
  print(Output.File)
  
  NotDone = File;
  
  if (is.null(NotDone) | !file.exists(NotDone)){
    warning("Invalid raw data file!")
    setwd(curr.dir);
    return("Invalid raw data file!");
  }
  
  Dirs <- list(NotDone)
  names(Dirs) <- getwd()
  
  #Sub.Dirs <- lapply(lapply(seq(length(Dirs)),function(x){unlist(gsub(paste0(".xls|.xlsx|.txt|.csv|",paste0(names(Dirs)[x],"-")),"",Dirs[[x]]))}),as.character);
  
  Fitted <- Manual.Dosage <- SS.Clusters <- NULL;
  
  BothData <- LoadAllData(Dirs = Dirs, sheet = sheet, want.list = TRUE)
  
  PAData <- lapply(BothData,function(each.plate){
    each.plate[which(sapply(seq(length(each.plate)),function(i)is.null(each.plate[[i]]$RQ)))]
  })
  PAData <- PAData[sapply(PAData,function(x)length(x)>0)]
  
  if (sum(unlist(sapply(PAData,length)),na.rm=TRUE) <= 0){
    warning("Bad input file! Contains RQ Values, which is for Dosage Calling jobs.")
    setwd(curr.dir)
    return("Bad input file! Contains RQ Values, which is for Dosage Calling jobs.")
  }
  
  dir.create(paste0(c(Output.Dir,'Plots/'),collapse="/"))
  Controls = read.csv(file=Control.File, header = TRUE, stringsAsFactors = FALSE)
  if (any(names(Controls) != c('Clone.plant_ID','Presence'))){
    warning("Bad Controls file! Must have only two columns: Clone-plant_ID & Presence")
    setwd(curr.dir)
    return("Bad Controls file! Must have only two columns: Clone-plant_ID & Presence")
  }
  names(Controls) <- c('Plant.ID',PAData[[1]][[1]]$SNP.Assay.Name[1])
  names(PAData[[1]]) <- Job
  CombinedPACalls <- CallAllPA(PAData=PAData, sheetname=sheet, Controls=Controls,
                               output.dir=paste0(c(Output.Dir,'Plots/'),collapse="/"), outlier.threshold=2, console=Console)
  
  write.csv(x=CombinedPACalls,file=Output.File,row.names=F,quote=F,na="");   
  
  setwd(curr.dir)
  dev.flush()
  return(CombinedPACalls)
}


#' Swaps any changed sample values
#'
#' Reads the Swaps.csv file in src and applies any substitions found
#' @param df The dataframe containing each of the models' values produced from CombineDosageCalls/CombinedPACalls
#' @param swap.dir The directory containing Swaps.csv
#' @param Job The name of the current job being called
#' @keywords Swap, Substitute, Swaps.csv, Tissue, Sampling
MakeSwaps <- function(df,swap.dir,Job){
  if (is.null(swap.dir))
    return(df)
  
  Swaps <- read.csv(swap.dir,stringsAsFactors=F,check.names=F)
  Swaps <- Swaps[which(Swaps$Job==Job),]
  if (nrow(Swaps)==0)
    return(df)
  for (each.row in seq(nrow(Swaps))){
    row <- Swaps[each.row,]
    vals <- row[which(names(row)%in%names(df))]
    vals <- vals[which(!is.na(vals))]
    df[intersect(which(df$Rack==row$Rack),which(df$DNAWell==row$DNAWell)),names(vals)] <- vals
  }
  return(df)
}

#' Gets the Predicted Dosage from Combined Dosage Calls
#'
#' Checks the different models' accuracy against the controls and then sets the predicted values to be the best model.
#' @param Summary.File Name of the current summary file
#' @param dose.dir Directory of where to find the Dosage.csv files. Should have sub-directory structure of dose.dir/<Job>/<PCR>.Dosage.csv
#' @param assays List of all possible assays. Default = c('Assay3','Assay2','Assay1','Assay3#','Assay1#')
#' @param PlateSz Size of the DNA plates excluding controls. For example, a 96 well DNA plate with 12 controls would have PlateSz = 84. Used to index the output.
#' @param NotDone Which PCR plates of the current job have yet to be added to the summary?
#' @param Job Name of the current job.
#' @keywords Combine, Output, All Models, Predicted
BuildSummary <- function(Job, NotDone, PlateSz=84, assays=c('Assay3','Assay2','Assay1','Assay3#','Assay1#'),
                         dose.dir = '/path/to/data', 
                         Summary.File = paste0(dose.dir,'../13_Summary/20',substr(Job,2,3),'_ProductionSummary.csv'),
                         header = c('RecId', 'Job', '#Samples','Name', 'Clone-plant_ID', 'Entry', 'TissID', 'Clone ID', 'Rack', 'RackOrder', 'Well Position', 'DNA Well', 'Pedigree','Pedigree Long',
                                    c(matrix(sapply(assays[-grep("#",assays)],function(x){c(paste0('(',x,')',c(' PCR',' RackOrder',' PCR Well',' A1 dRn', ' A2 dRn')),paste0(x,c('','_Pr','_QC')))}),nrow=1)),
                                    c(matrix(sapply(assays[grep("#",assays)],function(x)c(paste0('(',x,')',c(' PCR',' RackOrder',' PCR Well',' CT',' CT PEF',' dCT',' RQ')),paste0(x,c('','_Pr','_QC')))),nrow=1),'Action','uid')),
                         num.cols = c('CT', 'CT.PEF', 'dCT', 'RQ', 'Allele1.Delta.Rn', 'Allele2.Delta.Rn',matrix(rbind(assays,paste0(assays,'_Pr')),nrow=1),paste0(assays,'_QC')),
                         Assay.Cols = c('PCR','RackOrder','PCR Well','A1 dRn', 'A2 dRn', 'CT','CT PEF','dCT','RQ'),
                         consol.cols.num = c(c(matrix(sapply(assays[-grep("#",assays)],function(x){c(paste0('(',x,')',c(' A1 dRn', ' A2 dRn')),paste0(x,c('','_Pr')))}),nrow=1)),
                                             c(matrix(sapply(assays[grep("#",assays)],function(x)c(paste0('(',x,')',c(' CT',' CT PEF',' dCT',' RQ')),paste0(x,c('','_Pr')))),nrow=1))),
                         consol.cols.str = c(c(matrix(sapply(assays[-grep("#",assays)],function(x){c(paste0('(',x,')',c(' PCR',' RackOrder',' PCR Well')))}),nrow=1)),
                                             c(matrix(sapply(assays[grep("#",assays)],function(x)c(paste0('(',x,')',c(' PCR',' RackOrder',' PCR Well')))),nrow=1)))){
  
  library(gdata)
  print(paste0('Building Summary File: ',Summary.File,'..................................'))
  API.dir <- paste0(dose.dir,'../../API/')
  Job.dir <- paste0(API.dir,'Jobs')
  
  Key = data.frame(
    Level = c(rep(c('Positive','Negative','UD'),2)),
    ID = c(as.numeric(c(0,1,2)),'Positive','Negative','UD'),
    stringsAsFactors = FALSE)
  
  sub <- data.frame(
    old = c('Allele1.Delta.Rn','Allele2.Delta.Rn','Plant.ID','Well.Position','Well','CT.PEF','Delta.Ct','Plate','Clone.ID','DNA.Well','DNAWell'),
    new  = c('A1 dRn','A2 dRn','Clone-plant_ID','PCR Well','Well Position','CT PEF','dCT','PCR','CloneID','DNA Well','DNA Well'),
    stringsAsFactors = FALSE)
  
  padding <- function(x){
    if (str_length(x)>=3)
      return('')
    return(rep('0',3-str_length(x)))
  }
  
  JobRequest <- readxl::read_xlsx(paste0(API.dir,'JobRequest.xlsx'),sheet=1)
  JobRequest$JobID.Full <- paste0('J',sapply(JobRequest$Date,function(x)substr(x,nchar(x)-1,nchar(x))),'-',JobRequest$Station,sapply(JobRequest$JobID,padding),JobRequest$JobID)
  
  df <- df.temp <- NULL
  if (file.exists(Summary.File)){
    df <- read.csv(Summary.File,stringsAsFactors=F,check.names=F,na.strings=c(NA,''))
    df.temp <- df[which(df$Job == Job),]
    df <- df[which(df$Job != Job),]
  }
  
  calls <- list.files(paste0(dose.dir,Job))[grep('.csv',list.files(paste0(dose.dir, Job)))]
  calls <- calls[sapply(calls,function(x)gsub('^PA.|^Dosage.|.csv$','',x)%in%gsub('.xls(x)*$','',NotDone))]
  for (each.col in assays){
    cols <- c(names(df.temp)[grep(paste0('\\(',each.col,'\\)'),names(df.temp))],paste0(each.col,c('','_Pr')))       
    indx <- which(df.temp[,paste0('(',each.col,') PCR')]%in%gsub('.xls(x)*$','',NotDone))
    if(length(indx)>0)
      df.temp[indx,cols] <- NA
  }
  
  if (length(calls)>0){
      Current.Row <- JobRequest[which(JobRequest$JobID.Full == Job),]
      bulk <- Current.Row$Bulk
      SampleSheet <- read.csv(paste0(API.dir,'Jobs/',Job,'/',Current.Row$`Sample Sheet Path`),stringsAsFactors=F,check.names=F)
      SampleSheet$Rack <- NA
      for (i in 1:ceiling(nrow(SampleSheet)/PlateSz/bulk)){
        indx <- min(((i-1)*PlateSz*bulk+1),nrow(SampleSheet),na.rm=T):min(c(i*PlateSz*bulk,nrow(SampleSheet)),na.rm=T)
        current <- Current.Row$Plate_Start+i-1
        if (bulk > 1)
          current <- paste0(current,letters[1:bulk])
        SampleSheet$Rack[indx] <- current
      }
      names(SampleSheet)[which(names(SampleSheet)%in%sub$old)] <- sub$new[match(names(SampleSheet)[which(names(SampleSheet)%in%sub$old)],sub$old)]
      
      for (each.file in calls){
        print(paste0(Job,'/',each.file))
        tryCatch({
          temp <- read.csv(paste0(dose.dir, Job,'/',each.file),stringsAsFactors=F,check.names=F)
          temp$`Clone-plant_ID`[union(which(temp$`Clone-plant_ID`==''),which(is.na(temp$`Clone-plant_ID`)))] <- paste0('NTC_',seq(length(temp$`Clone-plant_ID`[union(which(temp$`Clone-plant_ID`==''),which(is.na(temp$`Clone-plant_ID`)))])))
          pa.assays <- assays[-grep('#',assays)][which(assays[-grep('#',assays)]%in%names(temp))]
          temp[pa.assays] <- apply(temp[pa.assays],2,function(x)Key$ID[match(x,Key$Level)])
          temp[num.cols[which(num.cols %in% names(temp))]] <- apply(temp[num.cols[which(num.cols %in% names(temp))]],2,as.numeric)
          names(temp)[which(names(temp)%in%sub$old)] <- sub$new[match(names(temp)[which(names(temp)%in%sub$old)],sub$old)]
          
          if (bulk > 1)
            temp <- interleave(temp,temp)
            
          job.dirs <- unlist(lapply(list.dirs(paste0(Job.dir,'/',Job)),function(x)paste0(x,'/',list.files(x))))
          pcr <- gsub('PA.|Dosage.|.csv','',each.file)
          info.file <- job.dirs[grep(pcr,job.dirs)]
          temp <- temp[!is.na(temp$Rack),]
          if (length(info.file)>0 && nrow(temp)>0){
            info.file <- info.file[-grep('Barcode',info.file)]
            info.file <- paste0(substr(info.file,0,max(str_locate_all(info.file,'/')[[1]][,1])),'Assignments.csv')
            assignments <-read.csv(info.file,stringsAsFactors=F,check.names=F)
            racks <- assignments[which(assignments$PCR == pcr),paste0('P',1:4)]
            if (bulk > 1)
              racks <- matrix(sapply(racks,paste0,letters[1:bulk]),ncol=1)
            info <- SampleSheet[which(SampleSheet$Rack %in% racks),]
            if (Current.Row$Field1 == "default"){
              non.controls <- 1:(min(which(as.numeric(temp$RecId)>80000000))-1)
              if (bulk > 1)
                temp[intersect(which(duplicated(temp$RecId)),non.controls),c('Entry','RecId','Clone-plant_ID')] <- info[match(temp$RecId[intersect(which(!duplicated(temp$RecId)),non.controls)],info$RecId)+1,c('Entry','RecId','Clone-plant_ID')] 
              temp$UID <- paste0(temp$Entry,temp$RecId,temp$`Clone-plant_ID`)
              info$UID <- paste0(info$Entry,info$RecId,info$`Clone-plant_ID`)
            }else{
              non.controls <- 1:(min(which(as.numeric(temp$TissID)>80000000))-1)
              if (bulk > 1)
                temp[intersect(which(duplicated(temp[,Current.Row$Field1])),non.controls),c(Current.Row$Field1,Current.Row$Field2)] <- info[match(temp[intersect(which(!duplicated(temp[,Current.Row$Field1])),non.controls),Current.Row$Field1],info[,Current.Row$Field1])+1,c(Current.Row$Field1,Current.Row$Field2)] 
              temp$UID <- paste0(temp[,Current.Row$Field1],temp[,Current.Row$Field2])
              info$UID <- paste0(info[,Current.Row$Field1],info[,Current.Row$Field2])
            }
            temp <- cbind(temp,info[match(temp$UID,info$UID),which(!(names(info)%in%names(temp)))])
          }
          
          present.assays <- assays[which(assays %in% names(temp))]
          present.assays <- present.assays[which(apply(temp[present.assays],2,function(x)sum(!is.na(x)))>0)]
          if (length(present.assays)>1){
            print(paste0('Too many assays for ',Job,'/',each.file,'!'))
          }
          names(temp)[which(names(temp)=='Plate')] <- 'PCR'
          names(temp)[which(names(temp)%in%Assay.Cols)] <- paste0('(',present.assays,') ',names(temp)[which(names(temp)%in%Assay.Cols)])
          if (length(grep('CN',present.assays))>0)
            temp[present.assays][is.na(temp[present.assays])] <- 'NTC'
          
          temp <- temp[header[which(header%in%names(temp))]]
          df.temp <- rbind.fill(df.temp,temp)
        },error = function(e){
          print(paste0('File ',each.file,' returned the error: ',as.character(e)))
          errfile <- paste0('/path/to/data../13_Summary/.src/ERROR_',gsub(' |:','_',Sys.time()))
          write.csv(x=temp,file=paste0(errfile,'_',each.file),quote=F,row.names=F)
          write.csv(x=as.character(e),file=paste0(errfile,'_output.csv'),quote=F,row.names=F)
        })
      }
      df.temp$TissID[which(as.numeric(df.temp$RecId)>80000000)] <- -seq(length(which(as.numeric(df.temp$RecId)>80000000)))
      df.temp <- df.temp[which(paste0(df.temp$RecId,df.temp$Entry,df.temp$`Clone-plant_ID`)!='NANANA'),]
      df.temp$uid <- paste0(df.temp$Job,',', df.temp$Rack,',',df.temp$`DNA Well`,',', df.temp$RecId,',',df.temp$`Clone-plant_ID`,',',df.temp$TissID,',',df.temp$Entry)
      cnts <- sapply(unique(df.temp$uid),function(x)length(which(df.temp$uid==x)))
      cnts <- names(cnts)[which(cnts > 3)]
      
      compact <- function(x,tag){
        if (length(x[!is.na(x)])==1)
          return(x[!is.na(x)])
        if (length(x[!is.na(x)])==0)
          return(NA)
        if (sum(as.numeric(x[!is.na(x)])-mean(as.numeric(x),na.rm=T),na.rm=T)!=0)
          if (length(x[!is.na(x)]) > 2)
            print(paste0('There are ',length(x[!is.na(x)]),' non-zero elements for ',tag))
        return(mean(x[!is.na(x)],na.rm=T))
      }
      consolidate <- function(y,tag){
        x <- compact(y,tag)
        #print(paste0('f(',paste0(y,collapse=', '),') :->',paste0(x,collapse=', ')))
        if (!is.na(as.numeric(x)))
          x <- as.numeric(x)
        return(rep(x,length(y)))
      }
      combine.row <- function(row,uid,ctr){
        row[which(names(row)%in%consol.cols.num)] <- apply(row[,which(names(row)%in%consol.cols.num)],2,function(each.col)consolidate(each.col,paste0(ctr,': ',uid)))
        row[which(names(row)%in%consol.cols.str)] <- apply(row[,which(names(row)%in%consol.cols.str)],2,function(each.col)consolidate(each.col,paste0(ctr,': ',uid)))
        #rbind(row.old,row[1,])
        return(row[1,])
      }
      
      df.compact <- NULL
      ctr <- 0

      df.temp <- as.data.frame(apply(df.temp,2,function(x){y<-x;y[which(x=='')]<-NA;return(y)}),stringsAsFactors = F)
      for (uid in unique(df.temp$uid)){
        ctr <- ctr+1
        row.all <- df.temp[which(df.temp$uid==uid),]
        col.cnts <- apply(row.all[,which(names(row.all)%in%consol.cols.str)],2,function(x)length(unique(x[!is.na(x)])))
        col.cnts <- col.cnts[which(col.cnts>1)]
        if (length(col.cnts)>0){
          for (each.col.name in names(col.cnts)){
            for (each.val in row.all[,each.col.name][!is.na(row.all[,each.col.name])]){
              row <- row.all[union(which(row.all[,each.col.name]==each.val),which(is.na(row.all[,each.col.name]))),]  
              compacted.row <- combine.row(row,uid,ctr)
              df.compact <- rbind.fill(df.compact,compacted.row)
            }
          }
        }else{
          compacted.row <- combine.row(row.all,uid,ctr)
          df.compact <- rbind.fill(df.compact,compacted.row)
        }
      }
      df.compact[,header[which(!(header%in%names(df.compact)))]] <- NA
      df.compact[,assays[-grep('#',assays)]] <- apply(df.compact[,assays[-grep('#',assays)]],2,function(x)Key$Level[match(x,Key$ID)])
      df <- rbind.fill(df,df.compact)
    }else{
    print(paste0('No calls for ',Job))
  }

  df.format <- as.data.frame(apply(df,2,function(x)gsub(',','|',x)),stringsAsFactors=F)
  df.format <- df.format[,header[which(header%in%names(df.format))]]
  
  year <- paste0('20',substr(Job,2,3))
  Summary.File <- paste0('/path/to/summary/',year,'_ProductionSummary.csv')
  write.csv(x=df.format,file=Summary.File,quote=F,na='',row.names=F)
  
  # Store a list of the bad controls
  tryCatch({
    assay.idx <- c(1,sapply(assays,function(x)which(header==x)),length(header))
    assay.idx <- sapply(seq(2,length(assay.idx)),function(i)(assay.idx[i-1]+1):assay.idx[i])
    bad.control.cols <- c(header[assay.idx[[1]]],unlist(sapply(seq(length(assays)),function(i)c(paste0(assays[i],'_Actual'),header[assay.idx[[i+1]]]))),'NumWrong')
    require(plyr)
    all_idx <- NULL
    bad_controls <- NULL
    temp <- df.format[which(as.numeric(df.format$RecId) < 0),]
    temp[,paste0(assays,'_Actual')] <- NA
    for (each.assay in assays){
      if (length(grep('#',each.assay))>0){
        Controls <- read.csv(paste0(Control.Dir,'Controls',gsub('#','CN',each.assay),'.csv'),stringsAsFactors=F)
        Controls <- Controls[,c('Name','Dosage')]
        Controls$Dosage[is.na(Controls$Dosage)] <- 'NTC'
        temp[,each.assay] <- as.numeric(temp[,each.assay])
      }else{
        Controls <- read.csv(paste0(Control.Dir,'ControlsPA.csv'),stringsAsFactors=F)
        cols <- c('Plant.ID',paste0(each.assay,'.1'),paste0(each.assay,'.2'))
        cols <- cols[which(cols%in%names(Controls))]
        Controls <- Controls[,cols]
        names(Controls) <- c('Name','Dosage')
        temp[,each.assay] <- Key$Level[match(temp[,each.assay],Key$ID)]
      }
      temp[,paste0(each.assay,'_Actual')] <- Controls$Dosage[match(temp$`Clone-plant_ID`, Controls$Name)]
      temp[intersect(which(temp[,paste0(each.assay,'_Actual')]=='NTC'),which(is.na(temp[,each.assay]))),each.assay] <- 'NTC'
      
    }
    
    temp$NumWrong <- apply(temp[,sapply(assays,function(x)c(paste0(x,c('','_Actual'))))],1,function(x)sum(x[assays]!=x[paste0(assays,'_Actual')],na.rm=T))
    bads <- temp[which(temp$NumWrong>0),bad.control.cols[which(bad.control.cols%in%names(temp))]]
    print(paste0('Printing ',nrow(bads),' bad controls to ',paste0(Control.Dir,'BadControls',year,'.csv')))
    write.csv(x=bads,file=paste0(Control.Dir,'BadControls',year,'.csv'),row.names=F,quote=F,na="")
  },error = function(e){
    print(as.character(e))
    write.table(x=Sys.time(),file=paste0(Control.Dir,'.BadCtrl.err'),append=T,row.names=F,quote=F,col.names=F,sep='\n')
    write.table(x=as.character(e),file=paste0(Control.Dir,'.BadCtrl.err'),append=T,row.names=F,quote=F,col.names=F,sep='\n')
    write.table(x="",file=paste0(Control.Dir,'.BadCtrl.err'),append=T,row.names=F,quote=F,col.names=F,sep='\n')
  })
    
  return(df.format)
}



#' Gets the Predicted Dosage from Combined Dosage Calls
#' 
#' Checks the different models' accuracy against the controls and then sets the predicted values to be the best model.
#' @param nulliplex.thres Final predictions with RQ < \code{nulliplex.thresh} will be set to nulliplex. Necessary the reaction fails for all nulliplex controls (very probable.). (Default = 0.2)
#' @param LL.min Any predictions with likelihood less than this will be ignored. Predictions with likelihood greater than this will be included in the final predictions weighted average. (Default=0.5)
#' @param Console Should the plots be printed to the console and not saved to \code{Output.Dir}? (Default=F).
#' @param Methods The different methods combined with \code{\link{CombineDosages}}. Must match \code{df}'s column names for these predictions.
#' @param Target The target assay. This is used primarily to rename the "Predicted" and "Probability" columns to the correct names. The assay names Assay1,Assay1#,Assay1CN, Assay1 CN, etc. will all generate a column named Assay1# and Assay1#_Pr.
#' @param Controls The Controls with Well information included
#' @param Output.Dir The directory to print the resulting plots to. Uneccesary if \code{Console=T}.
#' @param Controls the list of Controls for the dataset \code{df}.
#' @param df The dataset built from joining the results of each method specified in methods. Use \code{\link{CombineDosages}}
#' @keywords Combine, Output, All Models, Predicted
GetPredicted <- function(df, Controls, Output.Dir, Target = 'Assay1',Methods = c('Linear','Multinomial','Manual','SemiSupervised','QS.Dosage'), Console=F,LL.min = 0.5,nulliplex.thresh=0.02){
  #Methods = c('Fit.All','RoundQS','Multinom','Manual.Call','SS.Dosage','Predicted')
  Methods <- Methods[which(Methods %in% names(df))]
  df$Method <- df$Probability <- df$Predicted <- NA
  
  #df$Controls <- Controls$Dosage[match(df$`Clone-plant_ID`,Controls$Name)]
  GetAccuracies <- function(df,Controls){
    All_Accs <- rbind.fill(lapply(unique(df$Plate),function(each.plate){
      df.new <- df[which(df$Plate == each.plate),]
      Accuracies <- data.frame(t(sapply(Methods,function(col){
        sum(df.new[,col]==df.new$Controls,na.rm = TRUE)/sum(!is.na(df.new$Control),na.rm=TRUE)  
      })),stringsAsFactors = FALSE)
    }))
    All_Accs$Plate <- unique(df$Plate)
    All_Accs$Sum <- apply(All_Accs[,Methods],1,function(x)sum(x,na.rm = TRUE))
    return(All_Accs)
  }
  
  # Used to do this Multiple times
  pred <- NULL
  All_Accs <- GetAccuracies(df,Controls)
  for (each.plate in All_Accs$Plate){
      print(each.plate)
      df.new <- df[which(df$Plate == each.plate),Methods]
      df.new <- apply(df.new,2,as.numeric)
      df.new <- df.new+1
      df.new[which(is.na(df.new))] <- 0
      mods <- as.numeric(All_Accs[which(All_Accs$Plate == each.plate),Methods])
      if (length(which(mods>LL.min))>0){
        df.new <- matrix(df.new[,which(mods>LL.min)],ncol=length(which(mods>LL.min)))
        temp <- round((df.new %*% mods[which(mods>LL.min)])/sum(All_Accs[which(All_Accs$Plate == each.plate),Methods[which(mods>LL.min)]],na.rm=T),0)
      }else{
        temp <- df.new[,which(mods==max(mods,na.rm=T))[1]]
      }
      temp[which(temp==0)]<-NA
      temp <- temp-1
      df$Predicted[which(df$Plate == each.plate)] <- temp
      df$Method[which(df$Plate == each.plate)] <- Methods[which(All_Accs[which(All_Accs$Plate==each.plate),Methods]==max(All_Accs[which(All_Accs$Plate==each.plate),Methods],na.rm=T))[1]]
      
      df.new <- df[which(df$Plate==each.plate),]
      for (i in rep(c(1,2,0,3,4),2)){
        means <- sapply(0:4,function(x){mean(df.new$RQ[which(df.new$Predicted==x)],na.rm=TRUE) })
        sds <- sapply(0:4,function(x){sqrt(var(df.new$RQ[which(df.new$Predicted==x)],na.rm=TRUE)) })
        min.rng <- sapply(0:4,function(x){means[x+1]*sds[x+1]})
        max.rng <- sapply(0:4,function(x){means[x+1]+2*sds[x+1]})
        prob <- 1*abs(pnorm(df.new$RQ[which(df.new$Predicted == i)],means[i+1],sds[i+1])-0.5)
        prob2 <- dnorm(df.new$RQ[which(df.new$Predicted==i)],means[i+1],sds[i+1])/dnorm(means[i+1],means[i+1],sds[i+1])
        df.new$Probability[which(df.new$Predicted == i)] <- apply(cbind(prob,prob2,prob2),1,function(x)mean(x,na.rm=T))
        df.new$Predicted[intersect(which(df.new$RQ > min.rng[i+1]),which(df.new$RQ < max.rng[i+1]))] <- i
        df.new$Predicted[intersect(which(df.new$RQ > min.rng[i+2]),which(df.new$Predicted==i))] <- min(c(i+1,4))
        df.new$Predicted[intersect(which(df.new$RQ < max.rng[i]),which(df.new$Predicted==i))] <- max(c(i-1,0))
        df.new$Predicted[which(!is.na(df.new$Controls))] <- df.new$Controls[which(!is.na(df.new$Controls))]
      }
      df.new$Predicted[which(df.new$RQ < nulliplex.thresh)] <- 0
      df.new$Predicted[intersect(which(df.new$Predicted==0),which(df.new$RQ > mean(c(nulliplex.thresh,mean(df.new$RQ[which(df.new$Predicted==1)],na.rm=T)),na.rm=T)))] <- 1
      
      df$Probability[which(df$Plate == each.plate)] <- df.new$Probability
      df$Predicted[which(df$Plate == each.plate)] <- df.new$Predicted
  }
  print(sum(df$Predicted == df$Controls,na.rm = TRUE)/sum(!is.na(df$Controls),na.rm = TRUE))
  
  df$Probability[is.na(df$RQ)] <- .9
  df$Predicted <- as.double(df$Predicted)
  df$Predicted[which(df$Predicted > 4)] <- 4
  df$Predicted[which(df$Predicted < 0)] <- 0
  df$Predicted[which(is.na(df$CT))] <- 0
  df$Predicted[intersect(which(is.na(df$CT)),which(is.na(df$CT.PEF)))] <- 'UD'
  df$RQ[which(df$Predicted==0)] <- 0
  df$RQ[which(df$Predicted=='UD')] <- -0.0001
  
  # Construct Likelihood of Greater than 1
  df['> 1'] <- 0
  idx <- order(df$RQ)[which(df$RQ[order(df$RQ)] >= quantile(df$RQ[order(df$RQ)][which(df$Predicted[order(df$RQ)]==1)],0.85,na.rm=T))]
  df$`> 1`[idx] <- cumsum(df$Probability[idx])
  df$`> 1` <- df$`> 1`/quantile(df$`> 1`[which(df$Predicted==3)],0.8,na.rm=T)
  df$`> 1`[which(df$`> 1`>=1)] <- 0.99

  graph.df=PrintGraph(df.all = df,Controls = Controls,dir = paste0(Output.Dir,"/Plots"),console=Console);
  box.df=PrintBarGraph(df.all = df,Controls = Controls,dir = paste0(Output.Dir,"/Plots"),console=Console);
  
  
  if (length(grep('Assay1',Target)) > 0){
    names(df)[which(names(df)=='Predicted')] <- "Assay1#"
    names(df)[which(names(df)=='Probability')] <- 'Assay1#_Pr'
    df['Assay1#_QC'] <- NA
  }else if (length(grep('Assay3',Target)) > 0){
    names(df)[which(names(df)=='Predicted')] <- "Assay3#"
    names(df)[which(names(df)=='Probability')] <- "Assay3#_Pr"
    df['Assay3#_QC'] <- NA
  }else if (length(grep('Assay2',Target)) > 0){
    names(df)[which(names(df)=='Predicted')] <- "Assay2#"
    names(df)[which(names(df)=='Probability')] <- "Assay2#_Pr"
    df['Assay2#_QC'] <- NA
  }else{
    names(df)[which(names(df)=='Predicted')] <- Target
    names(df)[which(names(df)=='Probability')] <- paste0(Target,'_Pr')
    df[paste0(Target,'_QC')] <- NA
  }
  
  return(df)
}


#' Combine Each of the Models
#'
#' Combine each of the different models into a single dataframe. You could add more models and combine them here as desired.
#' @param Plates Names of the individual excel datasets within each Job. Usually <Plate Name>.xlsx
#' @param Jobs The names of each of the folders for each job. 
#' @param SemiSuper The semi-supervised clustering model's output. See \code{\link{MyClustering}} for more information
#' @param Manual The Manual calling simulated model's output. See \code{\link{ManualCall.Trials}} for more information
#' @param Dosages The output from \code{\link{GetDosages}} containing the SVM, Multinomial, and QS dosage predictions.
#' @keywords Combine, Output, All Models
CombineDosages <- function(Dosages, Manual, SemiSuper, Target, Jobs = names(Dirs),
                           Plates = lapply(lapply(seq(length(Dirs)),function(x){unlist(gsub(paste0(".xls|.xlsx|.txt|.csv",paste0(names(Dirs)[x],"-")),"",Dirs[[x]]))}),as.character)){
  
  ID.Cols <- c('Job','Plate','Sample.Name','Sample')
  get.uid <- function(df){
    if (is.null(df))
      return(df)
    cols <- ID.Cols[which(ID.Cols%in%names(df))]
    df$UID <- apply(df[,cols],1,function(x)paste0(x,collapse=''))
    return(df)
  }
  
  Dosages <- get.uid(Dosages)
  Manual <- get.uid(Manual)
  SemiSuper <- get.uid(SemiSuper)
  
  
  df <- Dosages
  for (each.col in names(Manual)[which(!(names(Manual)%in%names(df)))]){
    df[,each.col] <- Manual[match(df$UID,Manual$UID),each.col]
  }
  for (each.col in names(SemiSuper)[which(!(names(SemiSuper)%in%names(df)))]){
    df[,each.col] <- SemiSuper[match(df$UID,SemiSuper$UID),each.col]
  }
  Name.Key <- data.frame(Old=c('Dosage','Multinom','Manual.Call','Labels'),
                         New=c('Linear','Multinomial','Manual','SemiSupervised'),stringsAsFactors=F)
  names(df)[which(names(df)%in%Name.Key$Old)] <- Name.Key$New[match(names(df)[which(names(df)%in%Name.Key$Old)],Name.Key$Old)]
  return(df)
}

#' Driver Function for Calling Semi-Supervised Clustering Algorithm
#'
#' This function call the semi-supervised clustering algorithm for each plate/job of a given dataset.
#' @param Controls data.frame containing the Control values for this experiment. Must have column names of Actual.Name, Name, Dosage; where Actual.Name is what is found exactly in the data sets, Name is the desired output name (often the same as Actual.Name), and Dosage is the Dosage for each control. If NULL then this algorithm is no different than K-Means clustering, and will not return correct dosage numbers..
#' @param Target The target gene for calling the dosage. Must follow the exact name in datasets. 
#' @param VarIndices Index of the columns that contain variables to be used in the clustering process. 
#' @param VarNames Names of the columns that contain variables to be used in the clustering process. Using this will replace VarIndices with the indices of the provided VarNames.
#' @param AllData data.frame containing the relevent data. There must be a column for each of the VarNames/VarIndices, and the control samples must exactly match the Contols data.frame.
#' @keywords Semi-Supervised, K-Means
GetSemiSuperClusters <- function(AllData, Target = "Assay1",Controls = NULL){
  if (is.null(Controls)){
    warning("Semi-Supervised Clustering Works Best with Controls. However, none were provided.")
  }
  output <- NULL;
  
  
  
  for( each.job in 1:length(AllData)){
    for (each.plate in 1:length(AllData[[each.job]])){
      # The names of the relevant fields used for semi-supervised clustering
      VarNames = c("CT", "Delta.Ct", "RQ")
      VarIndices <- sapply(VarNames,function(x){which(names(AllData[[each.job]][[each.plate]]) == x)})
      
      tryCatch({
        temp <- suppressWarnings(data.frame(apply(AllData[[each.job]][[each.plate]][which(AllData[[each.job]][[each.plate]]$Target.Name == Target),VarIndices],2,as.double)));
        names(temp) <- names(AllData[[each.job]][[each.plate]])[VarIndices]
        Labels <- rep(NA,nrow(temp));
        Labels[which(AllData[[each.job]][[each.plate]]$Plant.ID[which(AllData[[each.job]][[each.plate]]$Target.Name == Target)] %in% Controls$Name)] <- 
          as.double(unlist(lapply(AllData[[each.job]][[each.plate]]$Plant.ID[intersect(which(AllData[[each.job]][[each.plate]]$Plant.ID %in% Controls$Name),
                                                                                          which(AllData[[each.job]][[each.plate]]$Target.Name == Target))],function(x){
                                                                                             Controls$Dosage[which(Controls$Name == x)]})))
        temp.clust <- MyClustering(temp,nClust = 5, Labels = Labels)
        GoodIndex <- apply(temp,1,function(x){sum(!is.na(x)) == length(x)})
        temp$Prob <- temp$Labels <- NA;
        temp$Prob[GoodIndex] <- temp.clust$Probs$Data$Prob;
        temp$Labels[GoodIndex] <- temp.clust$df$temp.labels
        temp$Job <- names(AllData)[each.job];
        temp$Plate <- names(AllData[[each.job]])[each.plate];
        temp$Sample <- AllData[[each.job]][[each.plate]]$Sample.Name[which(AllData[[each.job]][[each.plate]]$Target.Name == Target)];
        temp$Plant.ID <- AllData[[each.job]][[each.plate]]$Plant.ID[which(AllData[[each.job]][[each.plate]]$Target.Name == Target)];
        temp$Clone.ID <- AllData[[each.job]][[each.plate]]$Clone.ID[which(AllData[[each.job]][[each.plate]]$Target.Name == Target)];
        output <- rbind(output, temp)
      },error = function(e){
        print(paste0('Semi-Supervised Clustering: ',e))
      })
    }
  }
  return(output)
}


#' Loads All Data Sets
#'
#' Helper Function for Call Dosage
LoadAllData <- function(Dirs = NULL, sheet = NULL, want.list = FALSE, strings.factor = FALSE){
  library(stringr)
  library(gdata)
  if (is.null(Dirs)){
    Dirs <- lapply(ls(),dir);
    names(Dirs) <- ls();
  }
  
  AllData <- data.frame();
  AllList <- list();
  retList <- list();
  curr.dir <- getwd();
  for (each.dir in 1:length(Dirs)){
    print(paste0("Loading: ",sheet, " from ", names(Dirs)[each.dir]))
    temp.lst <- list();
    temp.names <- list();
    for (each.file in Dirs[[each.dir]]){
      setwd(names(Dirs)[each.dir]);
      filename <- each.file;
      #filename <- paste(names(Dirs)[each.dir],each.file,sep = "/");
      print(filename)
      tryCatch({
        wb <- suppressWarnings(read.xls(xls = filename, sheet = sheet, stringsAsFactors = strings.factor));
        wb <- wb[which(wb[,4] != ""),];
        names(wb) <- gsub(" ",".",wb[1,]);
        wb <- wb[-1,];
        temp.names[[length(temp.names)+1]] <- gsub(paste0(".xls|.xlsx|.txt|.csv|",paste0(names(Dirs)[each.dir],"-")),"",each.file)
        temp.lst[[length(temp.lst) + 1]] <- wb;
      },
      error = function(cmd){
        print(paste("Loading File",each.file,"Caused the Error: ",cmd));
      })
      setwd(curr.dir);
    }
    names(temp.lst) <- temp.names;
    retList[[length(retList) + 1]] <- temp.lst;
  }
  
  names(retList) <- names(Dirs)
  if (want.list){
    return(retList)
  }else{
    return(ConvertList(retList))
  }  
}

#' Converts List into data Frame
#'
#' Helper Function for Call Dosage
ConvertList <- function(AllList){
  library(plyr)
  AllData<-lapply(seq(length(AllList)),function(i){
    x <- AllList[[i]]; 
    df <- lapply(seq(length(x)),function(j){
      y <- x[[j]];
      y$Plate <- names(x)[j]
      y$Job <- names(AllList)[i];
      return(y)
    })
    names(df) <- names(x)
    return(df)
  })
  names(AllData) <- names(AllList)
  return(rbind.fill(lapply(AllData,rbind.fill)))
  
  Headers <- NULL;
  for (each.job in AllList){
    for (each.plate in each.job){
      Headers <- union(Headers, names(each.plate));
    }
  }
  
  AllData <- NULL;
  for (each.job in 1:length(AllList)){
    for (each.wb in 1:length(AllList[[each.job]])){
      temp <- as.data.frame(matrix(NA,nrow = nrow(AllList[[each.job]][[each.wb]]),ncol = (length(Headers))+2))
      names(temp) <- union(c("Job","Plate"),Headers);
      temp$Job <- rep(names(AllList)[each.job],nrow(AllList[[each.job]][[each.wb]]));
      temp$Plate <- rep(names(AllList[[each.job]])[each.wb],nrow(AllList[[each.job]][[each.wb]]));
      for (each.name in Headers){
        print(each.name)
        col = AllList[[each.job]][[each.wb]][,which(names(AllList[[each.job]][[each.wb]]) == each.name)]
        if (length(col)){
          temp[,which(Headers == each.name)+2] <- col
        }
      }
      AllData <- rbind(AllData,temp)
    }
  }
  names(AllData) <- union(c("Job","Plate"),Headers);
  
  return(AllData)
}

#' Removes Duplicated Names in List
#'
#' Helper Function for Call Dosage
RemoveDuplicateNames <- function(data, Dirs, Controls = NULL, Jobs = names(Dirs), Plates = lapply(lapply(seq(length(Dirs)),function(x){unlist(gsub(paste0(".xls|.xlsx|.txt|.csv|",paste0(names(Dirs)[x],"-")),"",Dirs[[x]]))}),as.character),
                                 SampleIndex = which(names(data[[1]][[1]]) == "Sample.Name"), Target.Name = "Assay1 CN "){
  df <- data
  for (each.job in 1:length(Jobs)){
    for(each.plate in 1:length(Plates[[each.job]])){
      target.index <- which(data[[each.job]][[each.plate]]$Target.Name == Target.Name);
      Samples <- data[[each.job]][[each.plate]][target.index,SampleIndex];
      for (each.sample in Samples){
        duplicates <- which(data[[each.job]][[each.plate]][target.index,SampleIndex] == each.sample);
        samples <- data[[each.job]][[each.plate]][target.index,][duplicates,SampleIndex]
        indices <- paste0(".",seq(length(samples)))
        if (length(data[[each.job]][[each.plate]][target.index,][duplicates,SampleIndex]) > 1){
          if (!is.null(Controls) && (samples[1] %in% Controls$Name) && !(paste0(samples,indices) %in% Controls$Name)){
            #Actual <- paste0(samples,indices);
            Name <- paste0(samples,indices);
            Dosage <- Controls$Dosage[which(Controls$Name == samples[1])][1];
            #Well <- Controls$Well[which(Controls$Name == samples[1])][1];
            Controls <- data.frame(rbind(Controls,data.frame(Name = Name,Dosage = Dosage)));
          }
          samples <- paste0(samples,indices);
          df[[each.job]][[each.plate]][target.index,][duplicates,SampleIndex] <- samples
        }
      }
      df[[each.job]][[each.plate]][-target.index,SampleIndex] <- df[[each.job]][[each.plate]][target.index,SampleIndex];
    }
  }
  
  if (!is.null(Controls))
    return(list(df = df, Controls = Controls))
  return(df);
}

#' Print Output Graph
#'
#' Helper Function for Call Dosage
GetPCH <- function(df.all,Controls){
  df.all$pch <- Controls$Dosage[match(df.all$`Clone-plant_ID`,Controls$Name)]+48
  df.all$pch[is.na(df.all$pch)] <- 16
  return(df.all)
}

#' Print Output Graph
#'
#' Helper Function for Call Dosage
PrintGraph <- function(df.all,Controls,dir = "",console = FALSE){
  library(scales);
  print(paste("Printing Standard Plots for",dir));
  Controls <- Controls[!duplicated(Controls$Name),];
  df.all <- GetPCH(df.all,Controls)
  Color.Key <- data.frame(dose=c('UD',NA,0,1,2,3,4),
                          col =c('yellow','orange','red','green','blue','cyan','magenta'),
                          col2=c('darkyellow','darkorange','darkred','darkgreen','darkblue','deepskyblue1','violet'),
                          stringsAsFactors = F)
  
  for (each.plate in unique(df.all$Plate)){
    bad <- T
    if(!console && !dir.exists(paste0(dir))){
      dir.create(paste0(dir));
    }
    if(!console)
      png(filename=paste0(dir,"/Dosage.",each.plate,".png"),width = 800,height = 600);

    tryCatch({
      df <- df.all[which(df.all$Plate == each.plate),];
      df.known <- df[which(df$pch != 16),];
      if(sum(is.na(df$RQ)) <= 0.9*length(df$RQ)){
        col <- Color.Key$col[match(df$Predicted,Color.Key$dose)]
        col[which(!is.na(df$Controls))] <- Color.Key$col2[match(df$Controls[which(!is.na(df$Controls))],Color.Key$dose)]
        col <- col[order(df$RQ)]
        cex <- rep(0.7,nrow(df))
        cex[which(!is.na(df$Controls))] <- 1
        cex <- cex[order(df$RQ)]
        plot(sort(df$RQ),col = alpha(col,df$Probability[order(df$RQ)]), cex=cex,
             pch = df$pch[order(df$RQ)], main = each.plate,xlab = "Index",ylab = "RQ",xlim = c(0,385),
             sub=paste0("Method: ",paste0(unique(df$Method),collapse=', ')));
        legend("bottomright",legend = c("UD",0:4,'Known'),col = c(7,0:4+2,1),lwd = 1,pch = c(rep(16,6),35));
        
        xlim <- c(0,nrow(df))
        ylim <- c(min(df$RQ,na.rm=T),max(df$RQ,na.rm=T))
        text(mean(xlim),ylim[1],  paste0("% Correct (Controls): ",round(sum(df.known$Predicted == df.known$pch-48,na.rm = TRUE)/sum(!is.na(df.known$pch)),4)))
        bad <- F
      }
    },error = function(e){
      print(paste0('Standard Plotting for ',each.plate,' Returned the error:\n ',as.character(e)))
    },finally = {
      if (bad){
        plot(1,1,cex=0,axes=F,xlab='',ylab='',
             main=each.plate,sub=paste0("Method: ",paste0(unique(df$Method),collapse=', ')));
        text(1,1,'Bad Data',cex=3)
      }
      if (!console)
        dev.off();
    })
  }
  #dev.off()
  return(df);
}

#' Print Output Bar/Histogram Graph
#'
#' Helper Function for Call Dosage
PrintBarGraph <- function(df.all,Controls,dir = "",console=FALSE){
  library(scales);
  print(paste("Printing Dosage Cluster Plots for",dir));
  Controls <- Controls[!duplicated(Controls$Name),];
  df.all <- GetPCH(df.all,Controls)
  
  for (each.plate in unique(df.all$Plate)){
    bad=T
    if(!console && !dir.exists(paste0(dir))){
      dir.create(paste0(dir));
    }
    if(!console)
      png(filename=paste0(dir,"/Dosage.box.",each.plate,".png"),width=800,height=700);
    tryCatch({
      # Prep one plate for plotting by making NTC = -1,
      # ordering the plate according to RQ, and setting
      # the 0 dosage RQ values to 0.
      df <- df.all[which(df.all$Plate == each.plate),]; 
      df <- df[order(df$RQ,na.last = FALSE),];
      df$Predicted[which(df$Predicted == "UD")] <- -1; 
      df$Predicted <-  as.double(df$Predicted);
      df$RQ[intersect(union(which(df$Predicted == 0),which(df$Predicted == -1)),which(is.na(df$RQ)))] <- 0
      df$col = df$Predicted+2; 
      df$col[which(df$pch != 16)] <- df$pch[which(df$pch != 16)]-46
      df$col[which(df$Predicted == -1)] = 7;
      df$cex = 0.7; df$cex[which(df$pch != 16)] <- 1.3; 
      df.known <- df[which(df$pch != 16),];
      df$Control[which(df$pch != 16)] <- 8;
      
      if(!console && !dir.exists(paste0(dir))){
        dir.create(paste0(dir));
      }
      plot(jitter(df$Predicted,0,.3),df$RQ,col = alpha(df$col,df$Probability),
           pch = df$pch,cex = df$cex, xlab = "Predicted Dosage", ylab ="RQ",xlim = c(-1.2,4.7),
           main = each.plate, sub=paste0("Method: ",paste0(unique(df$Method),collapse=', ')));
      df$Predicted[which(is.na(df$Predicted))] <- -1
      bins <- sort(unique(df$Predicted[which(!is.na(df$RQ))]))
      boxplot(RQ~Predicted,data=df,boxwex=0.2,border=6,add=T, varwidth=T, at=bins+.5, xaxt='n')
      legend("bottomright",legend = c("UD",0:4,'Known'),col = c(7,0:4+2,1),lwd = 1,pch = c(rep(16,6),35));
      
      xlim <- c(-1,4)
      ylim <- c(min(df$RQ,na.rm=T),max(df$RQ,na.rm=T))
      text(xlim[1]+0.1*(xlim[2]-xlim[1]), ylim[2],  paste0("% Correct (Controls): ",round(sum(df.known$Predicted == df.known$pch-48,na.rm = TRUE)/sum(!is.na(df.known$pch)),4)))
      bad=F
      if(!console)
        dev.off();
    },error = function(e){
      print(paste0('Plotting the Box Plot for ',each.plate,' Returned the Error:',e))
      #dev.off()
    },finally = {
      if (bad){
        plot(1,1,cex=0,axes=F,xlab='',ylab='',
             main=each.plate,sub=paste0("Method: ",paste0(unique(df$Method),collapse=', ')));
        text(1,1,'Bad Data',cex=3)
      }
    })
  }
  return(df);
}

