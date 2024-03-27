library(pacman)
p_load(eyelinker, intervals, stringr, readr, dplyr, ggplot2, tidyverse, stringr, emmeans, sjPlot, nlme, patchwork, psych, zoo, magrittr, gazer, ggridges, cowplot)
# https://github.com/dmirman/gazer/blob/master/vignettes/Pupil-vignette.Rmd
# https://github.com/dmirman/gazer/blob/master/vignettes/blink_detection.Rmd

subjects <- commandArgs(trailingOnly = TRUE)

for (s in subjects) {
  base_dir <- "/Users/bai/Library/CloudStorage/Box-Box/CAMSLab_Projects/Studies/AlcGen/Data/In_Lab"
  behav_dir <- paste0(base_dir, "/processed_behav")
  subj_dir <- paste0(base_dir, "/pupillometry/restructured/trial_split/", s)
  preproc_dir <- paste0(base_dir, "/pupillometry/preprocessed/", s)
  
  check_and_create <- function(dir) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
      cat("Created directory:", dir, "\n")
    } else {
      cat("Directory already exists:", dir, "\n")
    }
  }
  
  check_and_create(subj_dir)
  check_and_create(preproc_dir)
  
  preproc_eg <- list()
  
  s_rate_500 <- c('L003', 'L004', 'L006', 'L007', 'L009', 'L020', 'L021', 'H001', 'L022', 'L024', 'L026', 'L027', 'L028', 'L029')
  s_rate_250 <- c('L013', 'L014', 'L015', 'L017', 'L019', 'H002', 'H003')
  s_rate_1000 <- c('L025', 'L031')
  
  if (s %in% s_rate_500) {
    time_adj <- 500
    amp_adj <- 1
  } else if (s %in% s_rate_250) {
    time_adj <- 250
    amp_adj <- 2
  } else {
    time_adj <- 1000
    amp_adj <- 0.5
  }
  
  load(paste0(subj_dir, "/", s, "_cond_split.RData"))
  
  task <- "conditioning"
  target_names <- c("cs_onset", "us_onset")
  for (name in target_names) {
    if (name %in% names(subject_phases_list)) {
      cat("  Starting processing phase ", name, "\n")
      
      # Get raw data in shape
      df <- subject_phases_list[[name]] %>%
        rename(trial = trial_num) %>%
        group_by(trial) %>%
        mutate(time = (row_number() - time_adj - 1) *  amp_adj / 1000) %>%
        ungroup() %>%
        select(subject, trial, time, xp, yp, ps)
      
      # Preprocessing starts here...
      pupil_formatted <- make_gazer(df,
                                    subject="subject",
                                    trial="trial",
                                    time = "time",
                                    x = "xp",
                                    y = "yp",
                                    pupil = "ps")
      
      # Discard trials or subjects in which too many pupils are missing (currently considering >.75)
      pup_missing <- count_missing_pupil(pupil_formatted,
                                         pupil = "pupil",
                                         missingthresh = 0.75)
      
      pup_missing_stat <- pup_missing %>%
        arrange(subject, trial) %>%
        group_by(subject, trial) %>%
        slice(1) %>%
        ungroup() %>%
        select(subject, trial, averageMissingSub, averageMissingTrial)
      
      preproc_eg[["miss"]] <- pup_missing_stat %>% filter(trial==1)
      preproc_eg[["raw"]] <- pup_missing %>% filter(trial==1)
      
      
      # Blink Detection
      # Eye links are characterized by a pronounced drop in the pupillary signal, followed by a full loss of signal, and then usually a recovery artifact when the signal comes back online. 
      
      # To efficiently detect blinks, smooth the signal using a weighted moving window average of 40 samples (20ms before and after; 1000hz) and calculate a velocity profile for each trial. 
      pupil_blink_algo <-  pup_missing %>%
        mutate(smooth_pupil = moving_average_pupil(pupil, n = 20))
      
      # Subtracting each sample from the immediately preceding sample in the signal
      pupil_blink_algo1 <- pupil_blink_algo %>%
        mutate(velocity_pupil=c(diff(smooth_pupil), NA))
    
      preproc_eg[["vel"]] <- pupil_blink_algo1 %>% filter(trial==1)
      
      # Blink onsets were subsequently identified as occurring when the velocity crossed a predetermined negative threshold (I selected -5 based on Mathot's (2013) recommendation. This rapid decrease in pupil diameter corresponds to the apparent decrease in size of the pupil as the eyelid closes. Likewise, when the eyelid reopens there is a recovery artifact wherein pupil size rapidly gets larger. Thus, the algorithm detected the recovery period by indexing the time since onset that velocity exceed some positive threshold (I selected -5, again based off Matĥot (2013) recommendation). Finally, the offset was detected as the time at which velocity fell back down to 0. In this way, a blink corresponds to an onset, recovery, and offset index. According to Mathôt (2013), this detection algorithm underestimates the blink period by several milliseconds, thus I selected a margin value (10 ms) which was subtracted from the onset and added to the offset.
      pupil_blink_algo2 <-  pupil_blink_algo1 %>%
        mutate(blinks_onset_offset = ifelse(velocity_pupil <= -5 | velocity_pupil >= 5, 1, 0)) %>%
        mutate(blinks_pupil = ifelse(blinks_onset_offset == 1, pupil == NA, pupil)) %>%
        mutate(extendpupil = extend_blinks(blinks_pupil, fillback = 10, fillforward = 10, hz = 500)) %>%
        mutate(interp = na.approx(extendpupil, na.rm = FALSE, rule=2))
      
      preproc_eg[["blink"]] <- pupil_blink_algo2 %>% filter(trial==1)
      
      # Detect artifacts which can arise from quick changes in pupil size (Kret & Sjak-Shie, in press). The max_dilation function calculates the normalized dilation speed, which is the max absolute change between samples divided by the temporal separation for each sample, preceding or succeeding the sample. To detect out liters, the median absolute deviation is calculated from the speed dilation variable, multiplied by a constant, and added to the median dilation speed variable--values above this threshold are then removed.
        
      pupil_artifact <- pupil_blink_algo2 %>% 
        group_by(subject, trial) %>% 
        mutate(speed = speed_pupil(interp, time)) %>%
        mutate(MAD = calc_mad(speed)) %>%
        filter(speed < MAD)
      
      preproc_eg[["artifact"]] <- pupil_artifact %>% filter(trial==1)
        
      # Smoothing with 20-point moving average and linear interpolation. Point can cary.
      pupil_smoothed <- smooth_interpolate_pupil(pupil_artifact,
                                                 pupil = interp,
                                                 extendblinks = TRUE,
                                                 step.first = "interp",
                                                 filter = "moving",
                                                 maxgap = Inf,
                                                 type = "linear",
                                                 hz = 500,
                                                 n = 20)
      
      # Baseline Correction
      baseline_median <- pupil_smoothed %>%
        group_by(subject, trial) %>%
        filter(time <= 0) %>%
        summarize(baseline_median = median(pup_interp)) %>%
        ungroup()
      
      preproc_eg[["baseline"]] <- baseline_median
      
      pupil_corrected <- pupil_smoothed %>%
        left_join(baseline_median, by = c("subject", "trial")) %>%
        mutate(pup_corrected = pup_interp - baseline_median) %>%
        ungroup()
      preproc_eg[["corrected"]] <- pupil_corrected %>% filter(trial==1)
      
      # Combine with behavior
      behav_data <- read.csv(file.path(paste0(behav_dir, "/combined_", task, "_data.csv"))) %>%
        filter(subject == s) %>%
        rename(trial = trial_num) %>%
        select(subject, phase, trial, stim, incong_loc, us_loc, us_img, iter, drink_2grp, trial_cat)
      
      pupil_preproc <- pupil_corrected %>%
        left_join(behav_data, by = c("subject", "trial"))
      
      # Save data
      save(pupil_preproc, file = file.path(paste0(preproc_dir, "/", s, "_preproc_", name, ".RData")))
      save(preproc_eg, file = file.path(paste0(preproc_dir, "/", s, "_preproc_interim_", name, ".RData")))
  
      cat("    Saving data... ","\n")
      cat("    Completed pupil preprocessing for subject ", s, " phase ", name, "\n")
      
    } # if statement
  } # for loop
} # subject loop