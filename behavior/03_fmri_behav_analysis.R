library(pacman)
library(viridisLite)
library(reshape2)
library(jsonlite)
library(tidyverse)
library(ggpubr)
library(hrbrthemes)
library(alluvial)
library(viridis)
library(patchwork)
library(wesanderson)
#library(factoextra)
library(BayesFactor)

#check OS ----
if (.Platform$OS.type == 'windows'){
  filesep <- .Platform$file.sep
  volume <- 'Z:'  #or the letter your windows mounted the research storage unter
} else if (Sys.info()["sysname"] == "Darwin") { #check if mac/OSX
  volume <- '/Volumes/psy_wfg_henke'
  filesep <- .Platform$file.sep
} else if (Sys.info()["sysname"] == "Linux") {
  volume <- 'storage/research/psy_wf_henke'
  filesep <- .Platform$file.sep
}

## create paths ----
basepath <- paste0(volume, "/s2019_twillems_fMRI_silent_engram/data/fMRI/")
raw <- "raw/"
bids <- "/bids/"
behav <- "/behav/"
output <- 'derivatives/behav/'
ses1 <- "ses-1"
func <- "func"
subject_tsv <- "participants.tsv"
specpath <- "sub-"
endpath <- "/behav/"

behavAnalysis <- function(df_detailed,group,output_dir){
  #turn strings into booleans
  df_detailed <- df_detailed %>% mutate(first.retrieval = if_else(first.retrieval == "correct", TRUE, FALSE),
                                        second.retrieval = if_else(second.retrieval=="correct",TRUE,FALSE),
                                        third.retrieval = if_else(third.retrieval=="correct",TRUE,FALSE),
                                        fourth.retrieval = if_else(fourth.retrieval=="correct",TRUE,FALSE),
                                        upperhalf        = if_else(upperhalf == "upper",TRUE,FALSE))

  #turn into long format
  df_detailed <- df_detailed[, c(1, 2, 19, 20, 3, 4, 5, 6 , 7, 8 , 9 , 10, 11, 12 , 13, 14, 15, 16, 17, 18 )]
  
  df_det_long <- df_detailed %>% rename(
    "Accuracy1" = "first.retrieval",
    "Accuracy2" = "second.retrieval",
    "Accuracy3" = "third.retrieval", 
    "Accuracy4" = "fourth.retrieval",
    "Button1" = "first.button",
    "Button2" = "second.button",
    "Button3" = "third.button",
    "Button4" = "fourth.button",
    "RT1" = "first.RT",
    "RT2" = "second.RT",
    "RT3" = "third.RT",
    "RT4" = "fourth.RT",
    "consc1" = "first.rating",
    "consc2" = "second.rating",
    "consc3" = "third.rating",
    "consc4" = "fourth.rating") %>% 
    
    pivot_longer(
      cols = 5:20,
      names_to=c(".value","retrieval"),
      names_sep = -1)
  
  #aggregation to count different trial types
  df_aggr2 <- df_detailed %>% filter(!is.na(second.rating) & !is.na(second.retrieval)) %>%
    mutate(ConsciousnessRatingRetrievalBinary = if_else(second.rating > 2,"high","low"))%>%
    group_by(id) %>%
    mutate(RT_mean = mean(second.RT),
           RT_sd = sd(second.RT),
           RTsplit = if_else(second.RT>=RT_mean,1,0),
           unclear_trials = if_else(second.retrieval ==0 & second.rating == 2, 1,0),
           consc_forgotten  = if_else(second.retrieval  == 0 & second.rating < 3, 1,0),
           forgotten_but_fast = if_else(consc_forgotten == 1 & second.RT < RT_mean,1,0)
    ) %>%
    summarise(.groups='drop_last',
              Accuracy_imm = mean(second.retrieval,na.rm = TRUE),
              Accuracy_imm_sem = sd(second.rating,na.rm = TRUE)/sqrt(n()),
              Accuracy_del = mean(third.retrieval,na.rm = TRUE),
              Accuracy_del_sem = sd(third.rating,na.rm = TRUE)/sqrt(n()),
              Accuracy_del_recog = mean(fourth.retrieval,na.rm = TRUE),
              Accuracy_del_recog_sem = sd(fourth.rating,na.rm = TRUE)/sqrt(n()),
              n = n())
  
  # Normalize the reaction times, ------
  df_det_long <- df_det_long %>%
    mutate(invRT = 1/RT, id = as.character(id))  %>%
    group_by(id, retrieval) %>%
    mutate(madRT = (invRT - median(invRT))/mad(invRT)) %>%
    ungroup() %>%
    filter(abs(madRT) < 5) %>% 
    group_by(id, retrieval) %>%
    mutate(madRT = (invRT - median(invRT))/mad(invRT)) %>%
    ungroup()
  
  df_agg_all <- df_det_long %>% group_by(id,retrieval,consc) %>% summarise(.groups = 'drop_last',Acc=mean(Accuracy),n=n()) %>% mutate(freq=n/sum(n))
  write.csv(df_agg_all,file = paste0(output_dir,"accuracies_per_sub_ret_consc",group,".csv"))
  df_agg_all_noconsc <- df_det_long %>% group_by(id,retrieval) %>% summarise(.groups = 'drop_last',Acc=mean(Accuracy),n=n()) %>% mutate(freq=n/sum(n))
  write.csv(df_agg_all_noconsc,file=paste0(output_dir,"accuracies_per_sub_ret",group,".csv"))
  

  # Section 2: CONTINUATION - creation of dataframes etc. -----
  df_det_long <- df_det_long %>% 
    mutate(consc = case_when(
      consc == 2 ~ "Unconscious",
      consc == 3 | consc == 0 ~ "Neither",
      consc == 4 ~ "Conscious"),
           combined_categories = case_when(
             Accuracy == 1 & consc == "Unconscious" ~ "Correct\nUnconscious",
             Accuracy == 1 & consc == "Neither" ~ "Correct\nInBetween",
             Accuracy == 1 & consc == "Conscious" ~ "Correct\nConscious",
             Accuracy == 0 & consc == "Unconscious" ~ "Incorrect\nUnconscious",
             Accuracy == 0 & consc == "Neither"   ~ "Incorrect\nInBetween",
             Accuracy == 0 & consc == "Conscious"  ~ "Incorrect\nConscious",),
    )
  
  buttonBias <- df_det_long %>%
    group_by(id, retrieval) %>% 
    filter(consc=="Unconscious") %>% 
    mutate(BBias1 = Button-2.5) %>% #2.5 should be mean if 2 and 3 was pressed equally
    summarise(BBias = mean(BBias1)) %>% 
    mutate(z_Bias = abs(BBias-mean(BBias))/sd(BBias))
  
  t_tests = buttonBias %>% 
    group_by(retrieval) %>% 
    summarise(P = t.test(BBias,mu=0)$p.value,
              Sig = ifelse(P < 0.05, "*"," "),
              MaxWidth = max(BBias))
  biasplot <- ggplot(data=buttonBias, aes(x=retrieval, y=BBias)) + 
    geom_boxplot() + 
    geom_text(aes(label=Sig, y = MaxWidth + 0.05), size = 15, data = t_tests) +
    geom_hline(yintercept = 0, linetype = 4,color="red")
  print(biasplot)
  print("no button Bias in recognition")
  
  buttonBias2 <- buttonBias %>%
    ungroup() %>%
    group_by(id) %>%
    summarise(BBiasAcrossSubj = mean(BBias)) %>%
    mutate(z_Bias = abs(BBiasAcrossSubj-mean(BBiasAcrossSubj))/sd(BBiasAcrossSubj))
  print(paste0("number of subjects with a buttonBias = z > 3: ",sum(buttonBias2$z_Bias >= 3)))
        
  
  by_id_ret <- df_det_long %>%  group_by(id, retrieval)  %>%
    select(id,Accuracy,retrieval,consc) %>% 
    summarise(.groups='drop_last',
              num_consc = sum(consc == "Conscious"),
              trues = mean(Accuracy),
              sd = sd(Accuracy)) 
  
  combined_cats <- df_det_long %>% group_by(retrieval, combined_categories) %>% 
    summarise(.groups='drop_last',N = n()) # %>% mutate(freq = N/sum(n()))
  
  subj_cats <- df_det_long %>%  group_by(retrieval, combined_categories, id,) %>%
    summarise(.groups='drop_last',N = n()) 
  
  df_trajectory <- df_det_long %>%
    select(retrieval,combined_categories,item,id ) %>%
    pivot_wider(id_cols = c(id,item),names_from = retrieval, names_prefix = 'retrieval',values_from = combined_categories) %>% 
    drop_na()
  
  df_traj_consc <- df_det_long %>% 
    select(retrieval,consc,item,id) %>% 
    pivot_wider(id_cols=c(id,item),names_from = retrieval, names_prefix = 'retrieval', values_from = consc) 
  df_alluv_consc <- df_traj_consc %>%
    group_by(retrieval1,retrieval2,retrieval3,retrieval4) %>% 
    summarise(.groups='drop_last',n=n()) %>% 
    drop_na() %>% 
    ungroup()
  df_alluv_consc <- df_alluv_consc[df_alluv_consc$n >= 15,]
  
  
  png(file=paste0(output_dir, "alluvial_",group,".png"),
      width=1118, height=781)
  alluv_plot <- alluvial(... = df_alluv_consc[c("retrieval2","retrieval3","retrieval4")], freq=df_alluv_consc$n, 
                         col = viridis(nrow(df_alluv_consc)),axis_labels =c("30min","24h Ret","24h Recog"),
                         cex=0.6,cex.axis = 1.25) 
  
  mtext(paste0("Alluvial Plot: Results for group: ", group), 3, line=3, font=2)
  dev.off()

  alluv_plot <- alluvial(... = df_alluv_consc[c("retrieval2","retrieval3","retrieval4")], freq=df_alluv_consc$n, 
                         col = viridis(nrow(df_alluv_consc)),axis_labels =c("30min","24h Ret","24h Recog"),
                         cex=0.6,cex.axis = 1.25) 
  mtext(paste0("Alluvial Plot: Results for group: ", group), 3, line=3, font=2)
  
  df_alluvial <- df_trajectory %>%  group_by(retrieval2,retrieval3,retrieval4) %>% 
    summarise(.groups='drop_last',n=n()) %>% 
    ungroup() 
  
  df_alluvial2 <- df_trajectory %>%  group_by(retrieval1,retrieval2,retrieval3,retrieval4) %>% 
    summarise(.groups='drop_last',n=n()) %>% 
    ungroup() %>% 
    filter(retrieval1 == "Incorrect\nConscious")
  
  df_alluvial3 = df_trajectory %>%  group_by(retrieval1,retrieval2,retrieval3,retrieval4) %>% 
    summarise(.groups='drop_last',n=n()) %>% 
    ungroup() %>% 
    filter(retrieval1 == "Correct\nConscious")

  
  if (nrow(df_alluvial2) != 0){
    png(file=paste0(output_dir, "alluvial_immediateIncor",group,".png"),
        width=1118, height=781)
    alluv_plot <- alluvial(... = df_alluvial2[c("retrieval1","retrieval2","retrieval3","retrieval4")], freq=df_alluvial2$n, 
                           col = viridis(nrow(df_alluvial2)),axis_labels =c("imm","30min","24h Ret","24h Recog"),
                           cex=0.6,cex.axis = 1.25) 
    
    mtext(paste0("Alluvial Plot: Results including immediate ret for group: ", group), 3, line=3, font=2)
    dev.off()
  }
    
  png(file=paste0(output_dir, "alluvial_immediateCorr",group,".png"),
      width=1118, height=781)
  alluv_plot <- alluvial(... = df_alluvial3[c("retrieval1","retrieval2","retrieval3","retrieval4")], freq=df_alluvial3$n, 
                         col = viridis(nrow(df_alluvial3)),axis_labels =c("imm","30min","24h Ret","24h Recog"),
                         cex=0.6,cex.axis = 1.25) 
  
  mtext(paste0("Alluvial Plot: Results including immediate ret for group: ", group), 3, line=3, font=2)
  dev.off()
  
  
  df_alluvial<-df_alluvial[!df_alluvial$n <= ceiling(sum(df_alluvial$n)/(6^3)),]
  df_alluvial<- df_alluvial %>% mutate(
    trajectories = case_when(
      retrieval2 == "Correct\nConscious" & 
        retrieval3 == "Correct\nConscious"  & 
        retrieval4 == "Correct\nConscious" ~ "A: Conscious\nEpisodic\nMemory",
      
      (retrieval2 == "Incorrect\nUnconscious"  |retrieval2 == "Correct\nUnconscious" | retrieval2 == "Correct\nInBetween") & 
        (retrieval3 == "Correct\nConscious") &
        retrieval4 == "Correct\nConscious" ~ "B: Conscious\nEpisodic\nMemory\nneeded\nConsolidation",
      
      (retrieval2 == "Incorrect\nConscious" | retrieval2=="Incorrect\nInBetween") & 
        (retrieval3 == "Correct\nConscious") &
        retrieval4 == "Correct\nConscious" ~ "C: CEMaR", #Conscious\nEpisodic\nafter\nReconsolidation
      
      (retrieval2 == "Correct\nConscious" | retrieval2 == "Correct\nInBetween")& 
        (retrieval3 == "Incorrect\nUnconscious" |retrieval3 == "Correct\nUnconscious" |retrieval3 == "Correct\nInBetween"|retrieval3 == "Incorrect\nInBetween") & 
        retrieval4 == "Correct\nConscious" ~ "D: Conscious\nEpisodic\nMemory\nneeded Cue",
      
      (retrieval2 == "Incorrect\nUnconscious" | retrieval2 == "Correct\nUnconscious" | 
         retrieval2 == "Correct\nInBetween" | retrieval2 == "Incorrect\nInBetween" |
       retrieval2 == "Incorrect\nConscious" | retrieval2 == "Correct\nUnconscious") & 
        (retrieval3 == "Incorrect\nUnconscious" | retrieval3 == "Correct\nUnconscious" |retrieval3 == "Correct\nInBetween" |retrieval3 == "Incorrect\nInBetween") & 
        (retrieval4 == "Correct\nConscious" | retrieval4 == "Correct\nInBetween")   ~ "E: Conscious\nFamiliarity",
      
      (retrieval2 == "Correct\nConscious")  & 
        (retrieval3 == "Incorrect\nUnconscious" | retrieval3 == "Correct\nUnconscious" |
           retrieval3 == "Incorrect\nConscious"| retrieval3 == "Inorrect\nInBetween" | retrieval3 == "Correct\nInBetween") & 
        (retrieval4 == "Incorrect\nUnconscious" |retrieval4 == "Correct\nUnconscious"|retrieval4 == "Correct\nInBetween"|retrieval4 == "Incorrect\nInBetween") ~ "F: FaiC ",#Forgotten\nafter\ninitial\nConsolidation",
      
      retrieval2 == "Correct\nUnconscious" & 
        retrieval3 == "Correct\nUnconscious" & 
        retrieval4 == "Correct\nUnconscious"  ~ "G: UEM",
      
      (retrieval2 == "Incorrect\nUnconscious" | retrieval2 == "Correct\nUnconscious" | retrieval2 == "Incorrect\nConscious") & 
        (retrieval3 == "Incorrect\nUnconscious" | retrieval3 == "Correct\nUnconscious"| retrieval3 == "Incorrect\nConscious") &
        retrieval4 == "Incorrect\nConscious" ~ "H: Forgotten with\nInterference",
      
      (retrieval2 == "Correct\nConscious") &
        (retrieval3=="Correct\nConscious") &
        (retrieval4=="Correct\nInBetween") ~ "I: Category Memory",
      
      (retrieval2 ==   "Incorrect\nUnconscious" | retrieval2 == "Correct\nUnconscious" | retrieval2 == "Incorrect\nConscious"| retrieval2 == "Incorrect\nInBetween" | retrieval2 == "Correct\nInBetween")  & 
        (retrieval3 == "Incorrect\nUnconscious" | retrieval3 == "Correct\nUnconscious" | retrieval3 == "Incorrect\nConscious"| retrieval3 == "Incorrect\nInBetween" | retrieval3 == "Correct\nInBetween") & 
        (retrieval4 == "Incorrect\nUnconscious" |retrieval4 == "Correct\nUnconscious"|retrieval4 == "Correct\nInBetween"|retrieval4 == "Incorrect\nInBetween") ~ "J: Forgotten",
      
      retrieval2 == "Incorrect\nConscious" & 
        (retrieval3 == "Correct\nInBetween" | retrieval3 == "InCorrect\nInBetween" | retrieval3 == "Incorrect\nConscious") &
           retrieval4 == "Correct\nInBetween" ~ "K: False\nMemory",
      
      TRUE ~ "undefined"
    )
  )
  

  ###### Section 2.1 Dataframes that divide by subject ----
  
  df_traj_by_subject <- df_det_long %>%
    select(retrieval,combined_categories,id, item) %>%
    pivot_wider(id_cols = c(id,item),names_from = retrieval, names_prefix = 'retrieval', values_from = combined_categories)
  df_traj_by_subject <- na.omit(df_traj_by_subject)
  df_alluvial_by_subj <- df_traj_by_subject %>% 
    group_by(id,retrieval2,retrieval3,retrieval4) %>%
    summarise(.groups='drop_last',
              n=n()) %>%
    ungroup()
  
  #Welche sollte man ausschliessen? Wir haben 6^3 mögliche Trajectories, also 216. Wenn eine bestimmte Trajectory zufällig auftreten sollte 
  df_alluvial_by_subj<-df_alluvial_by_subj[!df_alluvial_by_subj$n <= ceiling(sum(df_alluvial_by_subj$n)/(6^3)),]
  df_alluvial_by_subj<- df_alluvial_by_subj %>% 
    mutate(
      trajectories = case_when(
        retrieval2 == "Correct\nConscious" & 
          retrieval3 == "Correct\nConscious"  & 
          retrieval4 == "Correct\nConscious" ~ "A: Conscious\nEpisodic\nMemory",
        
        (retrieval2 == "Incorrect\nUnconscious"  |retrieval2 == "Correct\nUnconscious" | retrieval2 == "Correct\nInBetween") & 
          (retrieval3 == "Correct\nConscious") &
          retrieval4 == "Correct\nConscious" ~ "B: Conscious\nEpisodic\nMemory\nneeded\nConsolidation",
        
        (retrieval2 == "Incorrect\nConscious" | retrieval2=="Incorrect\nInBetween") & 
          (retrieval3 == "Correct\nConscious") &
          retrieval4 == "Correct\nConscious" ~ "C: CEMaR", #Conscious\nEpisodic\nafter\nReconsolidation
        
        (retrieval2 == "Correct\nConscious" | retrieval2 == "Correct\nInBetween")& 
          (retrieval3 == "Incorrect\nUnconscious" |retrieval3 == "Correct\nUnconscious" |retrieval3 == "Correct\nInBetween"|retrieval3 == "Incorrect\nInBetween") & 
          retrieval4 == "Correct\nConscious" ~ "D: Conscious\nEpisodic\nMemory\nneeded Cue",
        
        (retrieval2 == "Incorrect\nUnconscious" | retrieval2 == "Correct\nUnconscious" | 
           retrieval2 == "Correct\nInBetween" | retrieval2 == "Incorrect\nInBetween" |
           retrieval2 == "Incorrect\nConscious" | retrieval2 == "Correct\nUnconscious") & 
          (retrieval3 == "Incorrect\nUnconscious" | retrieval3 == "Correct\nUnconscious" |retrieval3 == "Correct\nInBetween" |retrieval3 == "Incorrect\nInBetween") & 
          (retrieval4 == "Correct\nConscious" | retrieval4 == "Correct\nInBetween")   ~ "E: Conscious\nFamiliarity",
        
        (retrieval2 == "Correct\nConscious")  & 
          (retrieval3 == "Incorrect\nUnconscious" | retrieval3 == "Correct\nUnconscious" |
             retrieval3 == "Incorrect\nConscious"| retrieval3 == "Inorrect\nInBetween" | retrieval3 == "Correct\nInBetween") & 
          (retrieval4 == "Incorrect\nUnconscious" |retrieval4 == "Correct\nUnconscious"|retrieval4 == "Correct\nInBetween"|retrieval4 == "Incorrect\nInBetween") ~ "F: FaiC ",#Forgotten\nafter\ninitial\nConsolidation",
        
        retrieval2 == "Correct\nUnconscious" & 
          retrieval3 == "Correct\nUnconscious" & 
          retrieval4 == "Correct\nUnconscious"  ~ "G: UEM",
        
        (retrieval2 == "Incorrect\nUnconscious" | retrieval2 == "Correct\nUnconscious" | retrieval2 == "Incorrect\nConscious") & 
          (retrieval3 == "Incorrect\nUnconscious" | retrieval3 == "Correct\nUnconscious"| retrieval3 == "Incorrect\nConscious") &
          retrieval4 == "Incorrect\nConscious" ~ "H: Forgotten with\nInterference",
        
        (retrieval2 == "Correct\nConscious") &
          (retrieval3=="Correct\nConscious") &
          (retrieval4=="Correct\nInBetween") ~ "I: Category Memory",
        
        (retrieval2 ==   "Incorrect\nUnconscious" | retrieval2 == "Correct\nUnconscious" | retrieval2 == "Incorrect\nConscious"| retrieval2 == "Incorrect\nInBetween" | retrieval2 == "Correct\nInBetween")  & 
          (retrieval3 == "Incorrect\nUnconscious" | retrieval3 == "Correct\nUnconscious" | retrieval3 == "Incorrect\nConscious"| retrieval3 == "Incorrect\nInBetween" | retrieval3 == "Correct\nInBetween") & 
          (retrieval4 == "Incorrect\nUnconscious" |retrieval4 == "Correct\nUnconscious"|retrieval4 == "Correct\nInBetween"|retrieval4 == "Incorrect\nInBetween") ~ "J: Forgotten",
        
        retrieval2 == "Incorrect\nConscious" & 
          (retrieval3 == "Correct\nInBetween" | retrieval3 == "InCorrect\nInBetween" | retrieval3 == "Incorrect\nConscious") &
          retrieval4 == "Correct\nInBetween" ~ "K: False\nMemory",
        
        TRUE ~ "undefined"
      )
    )
  
  
  unconscious_RTs_bysub <- df_det_long %>% filter(consc=='Unconscious') %>% 
    group_by(id,retrieval,Accuracy) %>% 
    summarise(.groups='drop_last',
      median_rt = median(madRT),
      mean_rt = mean(madRT),
      sem_rt  = sd(madRT)/sqrt(n()))
  
  unconscious_RTs_deltas <- unconscious_RTs_bysub %>% ungroup() %>% group_by(id,retrieval) %>% summarise(.groups='drop_last',
    delta_RTs_mean   =mean_rt[Accuracy==FALSE]-mean_rt[Accuracy==TRUE],
    delta_RTs_median =median_rt[Accuracy==FALSE]-median_rt[Accuracy==TRUE],)
  
  
  
  consc_accuracies <- df_det_long %>% 
    select(consc,Accuracy,retrieval, id) %>% 
    filter(retrieval>1) %>% 
    group_by(consc, id) %>% 
    summarise(.groups='drop_last',
              acc_by_consc = mean(Accuracy)) 
  
  ### Section3: Descriptive PLOTS ----
  cor1to2 <- ggscatter(df_aggr2, x = "Accuracy_imm", y = "Accuracy_del", 
                       add = "reg.line", conf.int = TRUE, 
                       cor.coef = TRUE, cor.method = "pearson",
                       xlab = "Accuracy Day 1 Retrieval", ylab = "Accuracy Day 2 Retrieval") +
    theme(text=element_text(size=18)) +
    ggtitle(paste0("Correlation of retrieval accuracies"))#, Results for group: ", group))
  print(cor1to2)
  
  
  cor1to3 <- ggscatter(df_aggr2, x = "Accuracy_imm", y = "Accuracy_del_recog", 
                       add = "reg.line", conf.int = TRUE, 
                       cor.coef = TRUE, cor.method = "pearson",
                       xlab = "Accuracy 30 min Retrieval", ylab = "Accuracy 24h Recognition") +
    ggtitle(paste0("Correlation of retrieval Accuracies, Results for group: ", group))
  print(cor1to3)
  
  cor2to3 <- ggscatter(df_aggr2, x = "Accuracy_del", y = "Accuracy_del_recog", 
                       add = "reg.line", conf.int = TRUE, 
                       cor.coef = TRUE, cor.method = "pearson",
                       xlab = "Accuracy 24h min Retrieval", ylab = "Accuracy 24h Recognition") +
    ggtitle(paste0("Correlation of retrieval Accuracies, Results for group: ", group))
  print(cor2to3)

  ## Unconscious Accuracy
  df_mean_acc_by_sub_unconscious_halfed <- df_det_long %>% filter(consc=='Unconscious') %>%  
    group_by(id, retrieval, upperhalf) %>% 
    summarise(.groups='drop_last',
              mean_acc_bysub = mean(Accuracy),
              sem_mean_acc_bysub = sd(Accuracy)/sqrt(n()),
              n = n())
  
  recogDf <- df_det_long %>% 
    mutate(pairs = substr(item,start=1,stop=3)) %>%
    filter(retrieval==4) %>%
    select(-c(group,retrieval,Button,RT,madRT,combined_categories,invRT,item)) %>% 
    pivot_wider( names_from = upperhalf,values_from = c(consc,Accuracy), names_prefix = "upperhalf_") %>% 
    mutate(dangerousPair = case_when((consc_upperhalf_TRUE == "Unconscious" & consc_upperhalf_FALSE == "Conscious") | (consc_upperhalf_TRUE == "Conscious" & consc_upperhalf_FALSE == "Unconscious")  ~ "dangerous",
                                                 TRUE ~ "no")) 
  upperHalfCheck <- recogDf %>% group_by(id,consc_upperhalf_TRUE,dangerousPair) %>% summarise(acc=mean(Accuracy_upperhalf_TRUE),n=n()) %>% filter(consc_upperhalf_TRUE == "Unconscious")
  lowerHalfCheck <- recogDf %>% group_by(id,consc_upperhalf_FALSE,dangerousPair) %>% summarise(acc=mean(Accuracy_upperhalf_FALSE),n=n()) %>% filter(consc_upperhalf_FALSE == "Unconscious")
  ttestUpperHalf <- t.test(x=upperHalfCheck[upperHalfCheck$dangerousPair=='dangerous',]$acc,y=upperHalfCheck[upperHalfCheck$dangerousPair=='no',]$acc)
  ttestLowerHalf <- t.test(x=lowerHalfCheck[lowerHalfCheck$dangerousPair=='dangerous',]$acc,y=lowerHalfCheck[lowerHalfCheck$dangerousPair=='no',]$acc)
  lowerHalfCheck$whichHalf ="lowerhalf"
  lowerHalfCheck = lowerHalfCheck %>% subset(select=(-consc_upperhalf_FALSE))
  upperHalfCheck$whichHalf = "upperhalf"
  upperHalfCheck = upperHalfCheck %>% subset(select=(-consc_upperhalf_TRUE))
  RecogCheck = rbind(upperHalfCheck,lowerHalfCheck)
  NewRecogCheck = RecogCheck %>% 
    mutate(dangerousPair=if_else(dangerousPair == "dangerous","pairConsc","pairNotConsc")) %>% 
    pivot_wider(names_from = dangerousPair,values_from = c(n,acc)) %>% 
    na.omit() %>% 
    write_csv("recogCheck.csv")
  
  
  df_mean_acc_by_sub_unconscious <- df_det_long %>% 
    filter(consc=='Unconscious') %>%  
    group_by(id, retrieval) %>% 
    summarise(.groups='drop_last',
              mean_acc_bysub = mean(Accuracy),
              sem_mean_acc_bysub = sd(Accuracy)/sqrt(n()),
              n = n()) 
  write.csv(df_mean_acc_by_sub_unconscious, paste0(output_dir,"mean_unconsc_accuracy_longFormat_",group,".csv"))
  
  df_mean_acc_by_sub <- df_det_long %>% 
    group_by(retrieval,consc) %>% 
    mutate(bign = n()) %>% 
    ungroup() %>% 
    group_by(id,retrieval,consc) %>% 
    summarise(.groups='drop_last',
              mean_acc_bysub = mean(Accuracy),
              sem_mean_acc_bySub = sd(Accuracy)/sqrt(n()),
              n_all = n())

  df_mean_acc_by_sub_all_plready <- df_mean_acc_by_sub %>% 
    group_by(retrieval,consc) %>% 
    mutate(zscores = (mean_acc_bysub-mean(mean_acc_bysub))/sd(mean_acc_bysub),
           consc = case_when(consc=='Conscious' ~ "Sure", 
                                   consc=='Neither'   ~ "Unsure",
                                   consc=='Unconscious'~"Guess",
                                   TRUE ~ "undefined")) %>%
    filter(zscores<abs(2)) %>% 
    ungroup() %>% 
    na.omit() %>% 
    filter(retrieval != "1") %>% 
    mutate(Task=case_when(retrieval=="2" ~ '30min ret', retrieval=="3" ~ '24h ret', retrieval==4 ~ 'recog'))
  
  df_mean_acc_by_sub_unconscious_plready <- df_mean_acc_by_sub_unconscious %>% 
    group_by(retrieval) %>% 
    mutate(zscores = (mean_acc_bysub-mean(mean_acc_bysub))/sd(mean_acc_bysub)) %>%
    filter(zscores<abs(2)) %>% 
    ungroup() %>% 
    mutate(Task=case_when(retrieval=="2" ~ '30min ret', retrieval=="3" ~ '24h ret', retrieval==4 ~ 'recog'))

  t_tests_mean_all = df_mean_acc_by_sub_all_plready %>% 
    group_by(Task,consc) %>% 
    summarise(P = t.test(mean_acc_bysub,mu=0.5)$p.value,
              df = t.test(mean_acc_bysub,mu=0.5)$parameter,
              M = mean(mean_acc_bysub),
              SD = sd(mean_acc_bysub),
              Tval = as.double(t.test(mean_acc_bysub,mu=0.5)$statistic),
              Sig = ifelse(P < 0.05, "*"," "),
              MaxWidth = max(mean_acc_bysub),
              n=n()) %>% 
    group_by(Task) %>% 
    mutate(percentage = n/sum(n),
           acc=round(M*100))

  t_tests_mean_unconsc = df_mean_acc_by_sub_unconscious_plready %>% 
    group_by(Task) %>% 
    summarise(P = t.test(mean_acc_bysub,mu=0.5)$p.value,
              Sig = ifelse(P < 0.05, "*"," "),
              MaxWidth = max(mean_acc_bysub))
  
  recog_unconsc_acc = df_mean_acc_by_sub_unconscious_plready %>% 
    filter(Task == 'recog') %>% 
    select(mean_acc_bysub)
  bf_anchor = ttestBF(recog_unconsc_acc$mean_acc_bysub,mu=0.5)
  print(group)
  print(bf_anchor)
  print(t_tests_mean_unconsc$P)
  
  df_mean_acc_by_sub_unconscious_plot <- ggboxplot(df_mean_acc_by_sub_unconscious_plready, x="Task",y="mean_acc_bysub",
                                                   color="Task",#title=paste0("Mean accuracy by subject in all unconscious retrievals, Results for group: ", group),
                                                   #palette = "viridis",
                                                   ylab = 'Guessing accuracy',
                                                   xlab = "",#task 
                                                   size = 1.25,width=0.9,) + 
                                                   scale_color_viridis(alpha=1,begin=0,end=0,discrete=TRUE) + 
                                                   #ggtheme = theme_pubclean(),)+
                                                   #ggtheme = theme_pubr() ) +
    geom_hline(yintercept = 0.5, linetype = 4) + 
    ylim(0.3,0.9) + 
    geom_jitter(size=2.25,alpha=0.6,aes(color=Task),) + 
    geom_text(aes(label=Sig, y = MaxWidth + 0.025), size = 15, data = t_tests_mean_unconsc) + 
    theme(legend.position = 'none',text = element_text(size = 18))#, text = element_text(size = 27.5))
  
  
  print(df_mean_acc_by_sub_unconscious_plot)
  ggsave(plot = df_mean_acc_by_sub_unconscious_plot,filename = paste0(output_dir, "mean_unconsc_accuracy_bySubj_",group,".pdf"))
  
  if (grepl("corrected",group)){
    cats_plot_ready <- df_det_long  %>% mutate(consc_paper = case_when(consc=='Conscious' ~ "Sure", 
                                                                       consc=='Neither'   ~ "Unsure",
                                                                       consc=='Unconscious'~"Guess",
                                                                       TRUE ~ "undefined")) %>% 
      mutate(consc_paper = factor(consc_paper, levels=c("Guess", "Unsure", "Sure"))) %>%
      na.omit() %>% 
      group_by(retrieval,consc_paper) %>% 
      summarise(n = n(), acc = mean(Accuracy), sdAcc = sd(Accuracy)) %>% 
      filter(retrieval!=1) %>% 
      group_by(retrieval) %>% 
      mutate(percentage = n/sum(n),nacc=round(acc*100))
  } else {
    cats_plot_ready <- df_det_long  %>% mutate(consc_paper = case_when(consc=='Conscious' ~ "Sure", 
                                                                       consc=='Neither'   ~ "Unsure",
                                                                       consc=='Unconscious'~"Guess",
                                                                       TRUE ~ "undefined")) %>% 
      mutate(consc_paper = factor(consc_paper, levels=c("Guess", "Unsure", "Sure"))) %>%
      group_by(id,retrieval,consc_paper) %>% 
      mutate(zscores = (Accuracy-mean(Accuracy))/sd(Accuracy)) %>%
      filter(zscores<abs(2)) %>% 
      ungroup() %>% 
      na.omit() %>% 
      group_by(retrieval,consc_paper) %>% 
      summarise(n = n(), acc = mean(Accuracy)) %>% 
      filter(retrieval!=1) %>% 
      group_by(retrieval) %>% 
      mutate(percentage = n/sum(n),nacc=round(acc*100))
  }
  
  freqs_ratings = df_det_long %>% 
    group_by(retrieval,id,consc) %>% 
    summarise(n = n(), .groups = 'drop_last') %>% 
    mutate(freq = n / sum(n)) %>% 
    group_by(retrieval,consc) %>% 
    summarise(meanFreq = mean(freq), sdFreq=sd(freq), .groups="drop_last")
    

  cats <- cats_plot_ready %>%  ggplot(aes(x=retrieval,y=percentage,fill=consc_paper)) + 
    geom_bar(stat ='identity') + 
    labs(x="",y="Relative frequency of ratings",fill='Rating:  ')  +
    scale_x_discrete(labels=c("1" = "Immediate", "2" = "30min ret","3" = "24h ret", "4"="recog"))  + 
    theme_pubr() +    
    theme(text=element_text(size=18)) +
    scale_fill_viridis(alpha=.75,discrete=TRUE,begin = 0 ,end = 0.75, option="viridis") #
  print(cats)
  ggsave(plot = cats,filename = paste0(output_dir, "categories_plotN_",group,".png"))

  
  write.csv(df_mean_acc_by_sub_unconscious_halfed,paste0(output_dir,'mean_unconsc_accuracy_long_halfed',group,'.csv'))
  df_mean_acc_by_sub_unconsc_by_ret_halfed <- df_mean_acc_by_sub_unconscious_halfed %>% 
    select(-sem_mean_acc_bysub) %>% 
    pivot_wider(id_cols = c(id,upperhalf),values_from = mean_acc_bysub,names_from = retrieval,names_prefix = "ret")
  write.csv(df_mean_acc_by_sub_unconsc_by_ret_halfed,paste0(output_dir, "mean_unconsc_accuravy_wideFormat_HALFED_",group,".csv"))
  
  df_mean_acc_by_sub_unconsc_by_ret <- df_mean_acc_by_sub_unconscious %>% 
    select(-sem_mean_acc_bysub) %>% 
    pivot_wider(id_cols = id,values_from = mean_acc_bysub,names_from = retrieval,names_prefix = "ret")
  write.csv(df_mean_acc_by_sub_unconsc_by_ret,paste0(output_dir, "mean_unconsc_accuravy_wideFormat_",group,".csv"))
  
  bbias_df = df_mean_acc_by_sub_unconsc_by_ret %>% 
    pivot_longer(cols=starts_with("ret"),names_prefix="ret",names_to="retrieval",values_to="acc") %>%  
    right_join(y=buttonBias,by=c("id","retrieval"))
  
  

  ## REACTION TIME 
  RT_deltaplot_ggboxplot <- ggboxplot(unconscious_RTs_deltas, x="retrieval",y="delta_RTs_mean",
                                   color="retrieval",palette = "jco",
                                   title= paste0("Unconscious RT deltas by subj, Results for group: ", group),
                                   add= "jitter") + geom_hline(yintercept = 0, linetype = 2)
  print(RT_deltaplot_ggboxplot) + stat_compare_means()

  rt_densities <- ggplot(na.omit(df_det_long), aes(madRT,group=Accuracy, color=Accuracy, fill=Accuracy)) + geom_density(alpha=0.4) + facet_grid(vars(consc)) +
    labs(subtitle = paste0("Densities of reaction times over all retrievals, separated by consciousness, results for group: ", group))
  plot(rt_densities)
  t_test_RT = t.test(x = unconscious_RTs_deltas[unconscious_RTs_deltas$retrieval==4,]$delta_RTs_median)
  
  myFontStyle = list(size = 18,
                     face = "bold", 
                     family = NULL)
  
  arrangedPlot <- ggarrange(cats,df_mean_acc_by_sub_unconscious_plot,ncol=2,nrow=1,labels=c("A","B"),font.label = myFontStyle)
  print(arrangedPlot)
  #ggsave(plot = arrangedPlot, width = 12, height = 6, dpi = 300,filename = paste0(output_dir, "behav_plot_combined_small_",group,".pdf"))
  ggsave(plot = arrangedPlot, width = 12, height = 5, dpi = 300,filename = paste0(output_dir, "behav_plot_combined_",group,".png"))
  ggsave(plot = arrangedPlot,filename = paste0(output_dir, "behav_plot_combined_small_",group,".png"))
  
}

behavVersions <- c("_corrected.csv","_maxTrials.csv") #  b <- "_withgroup.csv"
b <- "_withgroup.csv"


sequences = c("hipp","wb")
for (behavVersion in behavVersions){
    #load different versions (with more or less data cleanup
      df_hipp <- read.csv(paste0(basepath, sequences[1], "/bids/derivatives/behav/df_", sequences[1], behavVersion))
      df_wb <- read.csv(paste0(basepath, sequences[2], "/bids/derivatives/behav/df_", sequences[2], behavVersion))
      df_detailed <- rbind(df_hipp, df_wb)
      output_dir <- paste0(basepath, "/behav/")
      write.csv(df_detailed,paste0(output_dir,"df_bothgroups",behavVersion))
    group <- paste0(sequence,tools::file_path_sans_ext(behavVersion))
    behavAnalysis(df_detailed,group,output_dir)
  }

