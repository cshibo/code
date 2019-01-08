# Data analysis for directed study
# Shibo Cao

# Research question: RT change as a function of CP location? 
# comparison between RW and CP

# import packages
library(data.table)
library(matrixStats)
library(gmodels)

####################################################################################
## DATA PROCESSING

# compiled data for all participants (run ONLY ONCE!)
Compiled_data_for_cp <- list()
Compiled_data_for_rw <- list()
# Change-point(first trial) RT for all participants
Compiled_first_trial_cp <- list()
Compiled_first_trial_rw <- list()
# Change-point(second trial) RT for all participants
Compiled_second_trial_cp <- list()
Compiled_second_trial_rw <- list()

# participant number
Participant_num <- c("141","466","472","475","476","478","479","480","481","482","486","490","491","492","493","495")

## loop for all RW data
for (i in 1:(length(Participant_num))){
  working_directory_RW <- "/Users/cdlab-admin/Documents/LRB/Shibo/new data/RW_"
  current_wd <- paste(working_directory_RW, Participant_num[i], sep = "")
  setwd(current_wd)
  files <- list.files(recursive = TRUE)
  files <- files[ !grepl("code", files) ] # eliminate other folders
  # an empty list for storing all RTs (based on bins)
  result_data <- list()
  # an empty list for storing all Change-point (first trial) RT
  cp_rt_data <- list()
  cp_rt2_data <- list()
  # loop for bins 
  for (j in 1:(length(files))){
    filename <- files[j]
    data <- read.csv(file = filename, header = TRUE)
    # index of change point within a trial
    cp_index <- which(data$ChangePoints == '[1]')
    if (cp_index[length(cp_index)] != nrow(data)){
      cp_index <- append(cp_index, nrow(data))
    }
    # change-point data
    cp_rt_first_half_index <- cp_index[1:(length(cp_index)%/%2)]
    cp_rt_second_half_index <- cp_index[((length(cp_index)%/%2) + 1):(length(cp_index))]
    cp_rt_first_half <- data[cp_rt_first_half_index,]
    cp_rt_second_half <- data[cp_rt_second_half_index,]
    cp_rt_data <- c(cp_rt_data, subset(cp_rt_first_half, select = "RTs"), subset(cp_rt_second_half, select = "RTs"))
    # trial 2 after change-point data
    cp_rt2_first_half_index <- cp_index[1:(length(cp_index)%/%2)] + 1
    cp_rt2_second_half_index <- cp_index[((length(cp_index)%/%2) + 1):(length(cp_index))] + 1
    cp_rt2_first_half <- data[cp_rt2_first_half_index,]
    cp_rt2_second_half <- data[cp_rt2_second_half_index,]
    cp_rt2_data <- c(cp_rt2_data, subset(cp_rt2_first_half, select = "RTs"), subset(cp_rt2_second_half, select = "RTs"))
    # number of chunks we divided based on the change-point
    cp_num <- length(cp_index) - 1
    # name of each column
    cp_num_char <- as.character(1:cp_num)
    cp_num_char <- paste("Bin",cp_num_char, sep="_")
    # an empty list for storing RTs
    df_row <- list()
    # loop for separating different change-point bins
    for(k in seq(1,cp_num)){
      x <- cp_num_char[k]
      range <- cp_index[k]:(cp_index[k+1] - 1)
      df_row[[k]] <- data$RTs[range]
    }
    result_data <- c(result_data, df_row)
  }
  
  # add the value to Compiled_first_trial variable 
  new_first_trial_var <- list()
  new_second_trial_var <- list()
  for (n in 1:4){
    new_first_trial_var[n] <- mean(unlist(cp_rt_data[n]))
  }
  new_first_trial_var <- unlist(new_first_trial_var)
  Compiled_first_trial_rw[[i]] <- new_first_trial_var
  # add the value to Compiled_second_trial variable 
  for (n in 1:4){
    new_second_trial_var[n] <- mean(unlist(cp_rt2_data[n]),na.rm = TRUE)
  }
  new_second_trial_var <- unlist(new_second_trial_var)
  Compiled_second_trial_rw[[i]] <- new_second_trial_var
  
  # maximum length of chunk (the number of CP positions possible)
  max_df_l <- max(lengths(result_data))
  # name of each chunk
  df_num_char <- as.character(1:max_df_l)
  df_num_char <- paste("Position",df_num_char, sep="_")
  
  # loop for making all elements in the list into consistent length
  for(l in seq(1,length(result_data))){
    lg <- as.numeric(length(result_data[[l]]))
    if (lg != max_df_l){
      result_data[[l]][(lg+1):max_df_l] <- rep("NA",(max_df_l-lg))
    }
  }  
  
  # convert the data into matrix and group them based on position
  df_grouped <- matrix(unlist(result_data), ncol = max_df_l, byrow=TRUE)
  class(df_grouped) <- "numeric"
  # loop a vector for the length of amount of RT available for each position
  l_pos <- double()
  for (m in 1:max_df_l){
    temp_table <- table(df_grouped[,m])
    if (length(temp_table) != length(result_data)){
      l_pos[m] <- length(df_grouped[,m]) - length(df_grouped[,m][df_grouped[,m]=="NA"])
    } else{
      l_pos[m] <- length(df_grouped[,m])
    }
  }
  
  # average each position's RTs
  result <- colMeans(df_grouped, na.rm = TRUE)
  sd <- colSds(df_grouped,na.rm = TRUE)
  # margin of error
  MoE <- qnorm(0.975)*sd/sqrt(l_pos) 
  
  # combining trial 6 and after
  Average_six_or_more <- mean(result[6:max_df_l])   # combined average for trial 6 or further
  SD_six_or_more <- sd(result[6:max_df_l])  # combined sd for trial 6 or further
  MoE_six_or_more <- qnorm(0.975)*SD_six_or_more/sqrt(sum(l_pos[6:max_df_l]))
  # all together average and standard deviation (and margin of error)
  Average_RT <- c(result[1:5], Average_six_or_more) 
  SD <- c(sd[1:5], SD_six_or_more)
  MOE <- c(MoE[1:5], MoE_six_or_more)
  # save average RT to the compiled variable (don't forget to change index every time)!!!!!!!!!!
  Compiled_data_for_rw[[i]] <- Average_RT
}


## loop for all CP data
for (i in 1:(length(Participant_num))){
  working_directory_CP <- "/Users/cdlab-admin/Documents/LRB/Shibo/new data/CP_"
  current_wd <- paste(working_directory_CP, Participant_num[i], sep = "")
  setwd(current_wd)
  files <- list.files(recursive = TRUE)
  files <- files[ !grepl("code", files) ] # eliminate other folders
  # an empty list for storing all RTs (based on bins)
  result_data <- list()
  # an empty list for storing all Change-point (first trial) RT
  cp_rt_data <- list()
  cp_rt2_data <- list()
  # an empty list for storing all Change-point (first trial) RT
  cp_rt_data <- list()
  cp_rt2_data <- list()
  # loop for bins 
  for (j in 1:(length(files))){
    filename <- files[j]
    data <- read.csv(file = filename, header = TRUE)
    # index of change point within a trial
    cp_index <- which(data$ChangePoints == '[1]')
    if (cp_index[length(cp_index)] != nrow(data)){
      cp_index <- append(cp_index, nrow(data))
    }
    # change-point data
    cp_rt_first_half_index <- cp_index[1:(length(cp_index)%/%2)]
    cp_rt_second_half_index <- cp_index[((length(cp_index)%/%2) + 1):(length(cp_index))]
    cp_rt_first_half <- data[cp_rt_first_half_index,]
    cp_rt_second_half <- data[cp_rt_second_half_index,]
    cp_rt_data <- c(cp_rt_data, subset(cp_rt_first_half, select = "RTs"), subset(cp_rt_second_half, select = "RTs"))
    # trial 2 after change-point data
    cp_rt2_first_half_index <- cp_index[1:(length(cp_index)%/%2)] + 1
    cp_rt2_second_half_index <- cp_index[((length(cp_index)%/%2) + 1):(length(cp_index))] + 1
    cp_rt2_first_half <- data[cp_rt2_first_half_index,]
    cp_rt2_second_half <- data[cp_rt2_second_half_index,]
    cp_rt2_data <- c(cp_rt2_data, subset(cp_rt2_first_half, select = "RTs"), subset(cp_rt2_second_half, select = "RTs"))
    # number of chunks we divided based on the change-point
    cp_num <- length(cp_index) - 1
    # name of each column
    cp_num_char <- as.character(1:cp_num)
    cp_num_char <- paste("Bin",cp_num_char, sep="_")
    # an empty list for storing RTs
    df_row <- list()
    # loop for separating different change-point bins
    for(k in seq(1,cp_num)){
      x <- cp_num_char[k]
      range <- cp_index[k]:(cp_index[k+1] - 1)
      df_row[[k]] <- data$RTs[range]
    }
    result_data <- c(result_data, df_row)
  }
  
  # add the value to Compiled_first_trial variable 
  new_first_trial_var <- list()
  new_second_trial_var <- list()
  for (n in 1:4){
    new_first_trial_var[n] <- mean(unlist(cp_rt_data[n]))
  }
  new_first_trial_var <- unlist(new_first_trial_var)
  Compiled_first_trial_cp[[i]] <- new_first_trial_var
  # add the value to Compiled_second_trial variable 
  for (n in 1:4){
    new_second_trial_var[n] <- mean(unlist(cp_rt2_data[n]),na.rm = TRUE)
  }
  new_second_trial_var <- unlist(new_second_trial_var)
  Compiled_second_trial_cp[[i]] <- new_second_trial_var
  
  # maximum length of chunk (the number of CP positions possible)
  max_df_l <- max(lengths(result_data))
  # name of each chunk
  df_num_char <- as.character(1:max_df_l)
  df_num_char <- paste("Position",df_num_char, sep="_")
  
  # loop for making all elements in the list into consistent length
  for(l in seq(1,length(result_data))){
    lg <- as.numeric(length(result_data[[l]]))
    if (lg != max_df_l){
      result_data[[l]][(lg+1):max_df_l] <- rep("NA",(max_df_l-lg))
    }
  }  
  
  # convert the data into matrix and group them based on position
  df_grouped <- matrix(unlist(result_data), ncol = max_df_l, byrow=TRUE)
  class(df_grouped) <- "numeric"
  # loop a vector for the length of amount of RT available for each position
  l_pos <- double()
  for (m in 1:max_df_l){
    temp_table <- table(df_grouped[,m])
    if (length(temp_table) != length(result_data)){
      l_pos[m] <- length(df_grouped[,m]) - length(df_grouped[,m][df_grouped[,m]=="NA"])
    } else{
      l_pos[m] <- length(df_grouped[,m])
    }
  }
  
  # average each position's RTs
  result <- colMeans(df_grouped, na.rm = TRUE)
  sd <- colSds(df_grouped,na.rm = TRUE)
  # margin of error
  MoE <- qnorm(0.975)*sd/sqrt(l_pos) 
  
  # combining trial 6 and after
  Average_six_or_more <- mean(result[6:max_df_l])   # combined average for trial 6 or further
  SD_six_or_more <- sd(result[6:max_df_l])  # combined sd for trial 6 or further
  MoE_six_or_more <- qnorm(0.975)*SD_six_or_more/sqrt(sum(l_pos[6:max_df_l]))
  # all together average and standard deviation (and margin of error)
  Average_RT <- c(result[1:5], Average_six_or_more) 
  SD <- c(sd[1:5], SD_six_or_more)
  MOE <- c(MoE[1:5], MoE_six_or_more)
  # save average RT to the compiled variable (don't forget to change index every time)!!!!!!!!!!
  Compiled_data_for_cp[[i]] <- Average_RT
}



##########################################################################################################
## DATA ANALYSIS

## analysis 1: Trials after change-point vs. average reaction time
# compile CP 
Compiled_CP <- matrix(unlist(Compiled_data_for_cp), ncol = 6, nrow = 16, byrow = TRUE)
Compiled_RW <- matrix(unlist(Compiled_data_for_rw), ncol = 6, nrow = 16, byrow = TRUE)

# average RT
Compiled_averaged_up_CP <- colMeans(Compiled_CP)
Compiled_averaged_up_RW <- colMeans(Compiled_RW)
output_rt <- rbind(Compiled_averaged_up_CP,Compiled_averaged_up_RW) # compiled CP and RW into matrix

# SD of RT
Compiled_sd_CP <- (colSds(Compiled_CP)) / sqrt(16)
Compiled_sd_RW <- (colSds(Compiled_RW)) / sqrt(16)
output_sd <- rbind(Compiled_sd_CP, Compiled_sd_RW) # compiled sd into matrix

# plot overall graph
x_axis <- 1:(ncol(output_rt))
x_labels = c("1","2","3","4","5","6 or more")
plot(0,0,xlim = c(1,6),ylim = c(0.5,1.5),xaxt = "n",type = "o", main = "Average Reaction Time for Two conditions", xlab = "Trials after Change-Point/Outlier", ylab = "Average reaction time (s)")
axis(1, at=1:6, labels=x_labels)
cl <- rainbow(2)
for (i in 1:2){
  lines(x_axis,output_rt[i,],col = cl[i],type = 'o')
}
# error bar
for (i in 1:2){
  for (j in 1:6){
    arrows(x_axis[j], (output_rt[i,j] - (output_sd[i,j]/2)), x_axis[j], (output_rt[i,j] + (output_sd[i,j]/2)), length=0.05, angle=90, code=3)
  }
}

# add legend to the graph
legend(4.5, y = 1.5, legend=c("CP", "RW"), col=cl, lty=1:2, cex=0.8, text.font=4)

# Paired t-test between CP and RW on each point
Compiled_CP <- matrix(unlist(Compiled_data_for_cp), ncol = 6, nrow = 16, byrow = TRUE)
Compiled_RW <- matrix(unlist(Compiled_data_for_rw), ncol = 6, nrow = 16, byrow = TRUE)

p1_cp <- Compiled_CP[,1]
p1_rw <- Compiled_RW[,1]
t.test(p1_cp, p1_rw, paired = TRUE, alternative = "two.sided")
p2_cp <- Compiled_CP[,2]
p2_rw <- Compiled_RW[,2]
t.test(p2_cp, p2_rw, paired = TRUE, alternative = "two.sided")
p3_cp <- Compiled_CP[,3]
p3_rw <- Compiled_RW[,3]
t.test(p3_cp, p3_rw, paired = TRUE, alternative = "two.sided")
p4_cp <- Compiled_CP[,4]
p4_rw <- Compiled_RW[,4]
t.test(p4_cp, p4_rw, paired = TRUE, alternative = "two.sided")
p5_cp <- Compiled_CP[,5]
p5_rw <- Compiled_RW[,5]
t.test(p5_cp, p5_rw, paired = TRUE, alternative = "two.sided")
p6_cp <- Compiled_CP[,6]
p6_rw <- Compiled_RW[,6]
t.test(p6_cp, p6_rw, paired = TRUE, alternative = "two.sided")

## analysis 2: First & second half of round 1 and 2 with reaction time
# compiled first-trial data
Compiled_first_trial_CP <- matrix(unlist(Compiled_first_trial_cp), ncol = 4, nrow = 16, byrow = TRUE)
Compiled_first_trial_RW <- matrix(unlist(Compiled_first_trial_rw), ncol = 4, nrow = 16, byrow = TRUE)
colnames(Compiled_first_trial_CP) <- c("R1_1st_Half", "R1_2nd_Half","R2_1st_Half","R2_2nd_Half")
colnames(Compiled_first_trial_RW) <- c("R1_1st_Half", "R1_2nd_Half","R2_1st_Half","R2_2nd_Half")


# average RT across all first trial
Compiled_averaged_first_trial_CP <- colMeans(Compiled_first_trial_CP)
Compiled_averaged_first_trial_RW <- colMeans(Compiled_first_trial_RW)
output_first_trial_rt <- rbind(Compiled_averaged_first_trial_CP,Compiled_averaged_first_trial_RW) # compiled CP and RW into matrix

# SD of RT
Compiled_sd_first_trial_CP <- (colSds(Compiled_first_trial_CP)) / sqrt(16)
Compiled_sd_first_trial_RW <- (colSds(Compiled_first_trial_RW)) / sqrt(16)
output_first_trial_sd <- rbind(Compiled_sd_first_trial_CP, Compiled_sd_first_trial_RW) # compiled sd into matrix

# plot overall graph
x_axis <- 1:(ncol(output_first_trial_rt))
x_labels = c("R1_1st_Half","R1_2nd_Half","R2_1st_Half","R2_2nd_Half")
plot(0,0,xlim = c(1,4),ylim = c(0.5,1.5),xaxt = "n", type = "o", main = "Progression of Change-point Reaction time", xlab = "Section of the task", ylab = "Average reaction time (s)")
axis(1, at=1:4, labels=x_labels)
cl <- rainbow(2)
for (i in 1:2){
  lines(x_axis,output_first_trial_rt[i,],col = cl[i],type = 'o')
}
# error bar
for (i in 1:2){
  for (j in 1:4){
    arrows(x_axis[j], (output_first_trial_rt[i,j] - (output_first_trial_sd[i,j]/2)), x_axis[j], (output_first_trial_rt[i,j] + (output_first_trial_sd[i,j]/2)), length=0.05, angle=90, code=3)
  }
}

# add legend to the graph
legend(3.25, y = 1.5, legend=c("CP", "RW"), col=cl, lty=1:2, cex=0.8, text.font=4)


## analysis 2 (Ext.): comparison of R1_1st_half with R2_2nd_half
Compiled_R1_1st_CP <- Compiled_first_trial_CP[,colnames = "R1_1st_Half"]
Compiled_R2_2nd_CP <- Compiled_first_trial_CP[,colnames = "R2_2nd_Half"]
Compiled_R1_1st_RW <- Compiled_first_trial_RW[,colnames = "R1_1st_Half"]
Compiled_R2_2nd_RW <- Compiled_first_trial_RW[,colnames = "R2_2nd_Half"]

# paired t-test
t.test(Compiled_R1_1st_CP, Compiled_R2_2nd_CP, paired = TRUE, alternative = "two.sided")
t.test(Compiled_R1_1st_RW, Compiled_R2_2nd_RW, paired = TRUE, alternative = "two.sided")


## analysis 3: First & second half of round 1 and 2 with reaction time
# compiled second-trial data
Compiled_second_trial_CP <- matrix(unlist(Compiled_second_trial_cp), ncol = 4, nrow = 16, byrow = TRUE)
Compiled_second_trial_RW <- matrix(unlist(Compiled_second_trial_rw), ncol = 4, nrow = 16, byrow = TRUE)
colnames(Compiled_second_trial_CP) <- c("R1_1st_Half", "R1_2nd_Half","R2_1st_Half","R2_2nd_Half")
colnames(Compiled_second_trial_RW) <- c("R1_1st_Half", "R1_2nd_Half","R2_1st_Half","R2_2nd_Half")


# average RT across all second trial
Compiled_averaged_second_trial_CP <- colMeans(Compiled_second_trial_CP)
Compiled_averaged_second_trial_RW <- colMeans(Compiled_second_trial_RW)
output_second_trial_rt <- rbind(Compiled_averaged_second_trial_CP,Compiled_averaged_second_trial_RW) # compiled CP and RW into matrix

# SD of RT
Compiled_sd_second_trial_CP <- (colSds(Compiled_second_trial_CP)) / sqrt(16)
Compiled_sd_second_trial_RW <- (colSds(Compiled_second_trial_RW)) / sqrt(16)
output_second_trial_sd <- rbind(Compiled_sd_second_trial_CP, Compiled_sd_second_trial_RW) # compiled sd into matrix

# plot overall graph
x_axis <- 1:(ncol(output_second_trial_rt))
x_labels = c("R1_1st_Half","R1_2nd_Half","R2_1st_Half","R2_2nd_Half")
plot(0,0,xlim = c(1,4),ylim = c(0.5,1.5),xaxt = "n", type = "o", main = "Progression of Change-point Reaction time at Trial 2", xlab = "Section of the task", ylab = "Average reaction time (s)")
axis(1, at=1:4, labels=x_labels)
cl <- rainbow(2)
for (i in 1:2){
  lines(x_axis,output_second_trial_rt[i,],col = cl[i],type = 'o')
}
# error bar
for (i in 1:2){
  for (j in 1:4){
    arrows(x_axis[j], (output_second_trial_rt[i,j] - (output_second_trial_sd[i,j]/2)), x_axis[j], (output_second_trial_rt[i,j] + (output_second_trial_sd[i,j]/2)), length=0.05, angle=90, code=3)
  }
}

# add legend to the graph
legend(3.25, y = 1.5, legend=c("CP", "RW"), col=cl, lty=1:2, cex=0.8, text.font=4)


## analysis 3 (Ext.): comparison of R1_1st_half with R2_2nd_half
Compiled_R1_1st_CP <- Compiled_second_trial_CP[,colnames = "R1_1st_Half"]
Compiled_R2_2nd_CP <- Compiled_second_trial_CP[,colnames = "R2_2nd_Half"]
Compiled_R1_1st_RW <- Compiled_second_trial_RW[,colnames = "R1_1st_Half"]
Compiled_R2_2nd_RW <- Compiled_second_trial_RW[,colnames = "R2_2nd_Half"]

# paired t-test
t.test(Compiled_R1_1st_CP, Compiled_R2_2nd_CP, paired = TRUE, alternative = "two.sided")
t.test(Compiled_R1_1st_RW, Compiled_R2_2nd_RW, paired = TRUE, alternative = "two.sided")







