
#' Get the Dosage for All Combinations of Controls
#'
#' Uses all possible subsets of control values in order to generate a calibration delta CT value. 
#' @param max The largest subset of controls to use for the calibration sample. Usually this is the total number of controls.
#' @param min The smallest subset of controls to use for the calibration sample.
#' @param Controls data.frame containing the Control values for this experiment. Must have column names of Actual.Name, Name, Dosage; where Actual.Name is what is found exactly in the data sets, Name is the desired output name (often the same as Actual.Name), and Dosage is the Dosage for each control.
#' @param Results.Only FALSE means that the AllData is actually a three-dimensonal list with the third dimension being each sheet in the data sets. If AllData only contains the "Results' section then set this to TRUE.
#' @param all.control.combos All possible combinations of the control variables. 
#' @param AllData A two-dimensional list seperated first by Jobs, then by Plates, containing the data for each experiment.
DosageOfAllCombos <- function(AllData, all.control.combos = NULL, Results.Only = FALSE, Controls = NULL, min = 1, max = nrow(Controls)){
  df <- data.frame();
  if (is.null(all.control.combos)){
    all.control.combos <- AllPossibleCalSamples(AllData, Results.Only = TRUE, Controls = Controls, min = min, max = max);
  }
  indx <- 0;
  for (each.combo in all.control.combos){
    indx <- indx + 1;
    df <- rbind(df,DosageOfGivenCombo(AllData,each.combo, Results.Only = Results.Only))
    print(paste0(indx,"/",length(all.control.combos),"   Complete"))
  }
  return(df)
}

#' Find the Dosage Given a Subset of Controls
#'
#' This function determines the calibration sample based on the provided subset of controls. It then computes the RQ values with that samples delta CT value and finds the dosage from those.
#' @param Num.Iters The number of iterations to run on the dosage calling algorithm. Each iteration refines the "by" variable by 10. See ManualCall.Trials for more information.
#' @param Target The reference gene for calling the dosage. Must follow the exact name in datasets.
#' @param Target The target gene for calling the dosage. Must follow the exact name in datasets.
#' @param Results.Only FALSE means that the AllData is actually a three-dimensonal list with the third dimension being each sheet in the data sets. If AllData only contains the "Results' section then set this to TRUE.
#' @param Controls data.frame containing the Control values for this experiment. Must have column names of Actual.Name, Name, Dosage; where Actual.Name is what is found exactly in the data sets, Name is the desired output name (often the same as Actual.Name), and Dosage is the Dosage for each control.
#' @param AllData A two-dimensional list seperated first by Jobs, then by Plates, containing the data for each experiment.
DosageOfGivenCombo <- function(AllData, Controls, Results.Only = FALSE, Target = "Assay1 CN ", Reference = "PEF CN", Num.Iters = 5){
  df <- data.frame();
  for(each.job in 1:length(AllData)){
    for(each.plate in 1:length(AllData[[each.job]])){
      if (Results.Only){
        temp <- AllData[[each.job]][[each.plate]];      
      }else{
        temp <- AllData[[each.job]][[each.plate]]$Results;
      }
      Assay1 <- temp[which(temp$Target.Name == Target),];
      PEF   <- temp[which(temp$Target.Name == Reference),];
      CB.dCt <- mean(as.double(Assay1$Delta.Ct[which(Assay1$Sample.Name %in% Controls$Name)]),na.rm = TRUE)*Controls$Multiplier[1];
      temp.df <- suppressWarnings(data.frame(Job = rep(names(AllData)[each.job],nrow(Assay1)), Plate = rep(names(AllData[[each.job]])[each.plate],nrow(Assay1)), CB = gsub("10RL-68-| ","",paste(Controls$Name,sep = ".",collapse = "/")),
                            Sample.Name = Assay1$Sample.Name, CT = as.double(Assay1$CT), CT.PEF = as.double(PEF$CT), dCT = as.double(Assay1$Delta.Ct), 
                            RQ = 2^(CB.dCt - as.double(Assay1$Delta.Ct)), actRQ <- Assay1$RQ, Dosage = round(2^(CB.dCt - as.double(Assay1$Delta.Ct))), actual = round(as.double(Assay1$RQ)),CB.dCt = CB.dCt))
      for (i in 1:length(Assay1$CT)){
        if (is.na(Assay1$CT[i])){
          if (is.na(PEF$CT[i])){
            temp.df$Dosage[i] <- -1;
          }else{
            temp.df$Dosage[i] <- 0;
          }
        }
      }
      Manual <- NA;
      Refine <- .1;
      Count <- 1;
      while((sum(is.na(Manual)) == length(Manual)) && Count < Num.Iters){
        Manual <- ManualCall.Trials(temp.df$RQ, by = Refine)
        Refine <- Refine / 10;
        Count <- Count + 1;
      }
      temp.df$Manual <- NA;
      temp.df$Manual[!is.na(temp.df$RQ)] <- Manual
      #temp.df$Manual[is.na(temp.df$Manaul)] <- 0;
      df <- rbind(df,temp.df);
    }
  }
  #names(df) <- c("Job","Plate","CB","Sample.Name","CT","CT.PEF", "dCT","RQ","actRQ","Dosage","actual","CB.dCt","Manual")
  return(df)
}

#' Determines All Calibration Sample Subsets
#'
#' This function finds all possible subsets from min to max size of the Control variables.
#' @param Control.Range The row range describing which samples are part of the control list.
#' @param max The largest subset of controls to use for the calibration sample. Usually this is the total number of controls.
#' @param min The smallest subset of controls to use for the calibration sample.
#' @param Controls data.frame containing the Control values for this experiment. Must have column names of Actual.Name, Name, Dosage; where Actual.Name is what is found exactly in the data sets, Name is the desired output name (often the same as Actual.Name), and Dosage is the Dosage for each control.
#' @param AllData A two-dimensional list seperated first by Jobs, then by Plates, containing the data for each experiment.
AllPossibleCalSamples <- function(AllData, Controls = NULL, min = 8, max = nrow(Controls), Control.Range = 85:95){
  
  if (is.null(Controls)){
    warning("No Control Variables Provided, Operation Stopped Before Completion.")
    return(NA)
  }
  
  all.control.combos <- list();
  for (each.subset in min:max){
    temp.lst <- combn(seq(nrow(Controls)),each.subset)
    for (each.col in seq(ncol(temp.lst))){
      df <- suppressWarnings(data.frame(Controls[temp.lst[,each.col],] ))
      names(df) <- names(Controls)
      df$Multiplier <- sum(as.double(df$Dosage),na.rm = TRUE)/length(temp.lst[,each.col])
      all.control.combos[[length(all.control.combos)+1]] <- df
    }
  }
  return(all.control.combos)
}

#' Measures the Accuracy of a Model vs. Control Dosages
#'
#' This measures the percentage of controls that agree with the predicted values.
#' @param Controls data.frame containing the Control values for this experiment. Must have column names of Actual.Name, Name, Dosage; where Actual.Name is what is found exactly in the data sets, Name is the desired output name (often the same as Actual.Name), and Dosage is the Dosage for each control.
#' @param AllDosages A serial data frame of all the relevant values including Job, Plate, RQ, CT, delta CT, and Sample Name.
ComputeScores <- function(AllDosages, Controls = NULL){
  if (is.null(Controls)){
    warning("No Control Values Provided, Exiting Now");
    return(NA);
  }
  Jobs <- unique(AllDosages$Job)
  for (each.job in Jobs){
    Plates <- unique(AllDosages$Plate[which(AllDosages$Job == each.job)])
    print(paste0("Printing:  Job ",which(Jobs == each.job),"/",length(Jobs)))
    for (each.plate in Plates){
      print(paste0("Printing:  Plate ",which(Plates == each.plate),"/",length(Plates)))
      temp <- AllDosages[intersect(which(AllDosages$Job == each.job), which(AllDosages$Plate==each.plate)),]
      for (each.CB in unique(temp$CB)){
        temp2 <- temp[which(temp$CB == each.CB),]
        score <- sum(temp2$Dosage[which(temp2$Sample.Name %in% Controls$Name)] == Controls$Dosage,na.rm = TRUE)/length(Controls$Name)
        AllDosages$Control.Score[intersect(which(AllDosages$CB == each.CB),
                                           intersect(which(AllDosages$Job == each.job), which(AllDosages$Plate==each.plate)))] <- score
      }
    }
  }
  return(AllDosages)
}

#' Measures the Confidence of An AllComboFit Output
#'
#' Generates the probability of each dosage call by comparing the mean and standard deviation of all fits from an AllComboFit output. 
#' @param df An AllComboFit output from DosageOfAllCombos
MeasureConfidence <- function(df){
  Jobs <- unique(df$Job)
  summary.df <- data.frame();
  first <- TRUE;
  indx <- 0;
  total.loops <- length(unique(df$Plate))
  Header <- c("Job","Plate","Sample","CT","CT.PEF","dCT","CB.dCT","RQ","Dosage","SD","P0","P1","P2","P3","P4","MaxP","MaxP.Val","NTC","N0","N1","N2","N3","N4","MaxN","MaxN.Val","All3");
  
  for (each.job in Jobs){
    Plates <- unique(df[which(df$Job == each.job),]$Plate)
    for (each.plate in Plates){
      print(paste0("Job: ",which(Jobs == each.job),"/",length(Jobs)," Plate: ",which(Plates == each.plate), "/",length(Plates), " %: ",indx/total.loops))
      if (!first){
        names(summary.df) <- Header
        write.csv(summary.df, "ConfidenceOfManualDosageAll.Trimmed.csv")
      }
      first <- FALSE;
      indx <- indx + 1;
      Samples <- unique(df[intersect(which(df$Plate == each.plate),which(df$Job == each.job)),]$Sample.Name)
      for (each.sample in Samples){
        temp <- df[intersect(intersect(which(df$Plate == each.plate),which(df$Sample.Name == each.sample)),which(df$Job == each.job)),]
        CT <- mean(temp$CT,na.rm = TRUE)
        CT.PEF <- mean(temp$CT.PEF,na.rm = TRUE)
        dCT <- mean(temp$dCT, na.rm = TRUE)
        CB.dCT <- mean(temp$CB.dCt,na.rm = TRUE)
        RQ <- mean(temp$RQ,na.rm = TRUE)
        SD <- sqrt(var(temp$RQ,na.rm = TRUE))
        
        
        num.df <- NULL
        prob.df <- NULL
        
        num.df <- c(length(temp$Dosage[which(temp$Dosage == -1)]))
        for (i in 0:4){
          deviation <- abs((i-RQ))
          p <- (pnorm(deviation)-0.5)/0.5
          prob.df <- c(prob.df, (1-p))
          num.df <- c(num.df,length(temp$Dosage[which(temp$Dosage == i)]))
        }
        prob.df <- prob.df*(1-SD)
        
        prob.df <- data.frame(t(prob.df))
        num.df <- data.frame(t(num.df))
        names(prob.df) <- 0:4
        names(num.df) <- -1:4
        
        if(length(which.max(prob.df)) == 0){
          MaxP <- NA
        }else{
          MaxP = as.double(names(prob.df)[which.max(prob.df)])
        }
        
        if(is.na(MaxP)){
          MaxN <- NA
        }else{
          MaxN =  as.double(names(num.df)[which.max(num.df)])
        }
        MaxP.Val <- as.double(max(prob.df,na.rm = TRUE))
        MaxN.Val <- as.double(max(num.df,na.rm = TRUE))
        
        temp.df <- suppressWarnings(data.frame(
          Job = each.job, Plate = each.plate, Sample = each.sample, CT = CT, CT.PEF =  CT.PEF, dCT = dCT, CB.dCT = CB.dCT, RQ = RQ, Dosage = round(RQ), SD = SD,
          prob.df, MaxP = MaxP, MaxP.Val = MaxP.Val, num.df, MaxN = MaxN, MaxN.Val = MaxN.Val))
        if (is.na(temp.df$CT)){
          if (is.na(temp.df$CT.PEF)){
            temp.df$Dosage <- temp.df$MaxP <- -1;
          }else if (is.na(temp.df$CT)){
            temp.df$Dosage <- 0;
          }
        }
        
        if (is.na(temp.df$Dosage) || is.na(temp.df$MaxP) || is.na(temp.df$MaxN)){
          temp.df$All3 <- NA;
        }else
          if (temp.df$Dosage == temp.df$MaxP && temp.df$MaxP == temp.df$MaxN){
            temp.df$All3 <- temp.df$Dosage
          }else{
            temp.df$All3 <- NA;
          }
        
        names(temp.df) <- Header
        summary.df <- rbind(summary.df,temp.df)
      }
    }
  }
  names(summary.df) <- Header
  return(summary.df)
  
}

