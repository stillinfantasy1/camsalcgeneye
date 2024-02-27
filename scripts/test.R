base_path <- "/Users/bai/Library/CloudStorage/Box-Box/CAMSLab_Projects/Studies/AlcGen/Data/In_Lab"


# Input
trial_split_path <- file.path(base_path, "pupillometry", "restructured", "trial_split")
# Output
preproc_path <- file.path(base_path, "pupillometry", "preprocessed")
if (!dir.exists(preproc_path)) {
  dir.create(preproc_path)
}

# Which subjects?
subject <- c('L007')

input_path <- file.path(trial_split_path, subject, paste0(subject, "_cond_split.RData"))
data <- load(input_path)
cs_onset <- subject_phases_list$cs_onset %>%
  rename(trial = trial_num) %>%
  group_by(trial) %>%
  mutate(time = row_number() - 500 - 1) %>%
  ungroup()


stimulu_onset <- cs_onset %>%
  filter(time >=0)

Sdata <- make_pupillometryr_data(data = stimulu_onset,
                                 subject = subject,
                                 trial = trial,
                                 time = time,
                                 condition = stim)

downsamp_data <- downsample_time_data(data = Sdata,
                                      pupil = ps,
                                      timebin_size = 10,
                                      option = 'mean')

# 4. Exclude trials with >75% missing values
missing <- calculate_missing_data(downsamp_data, pupil = ps)
excl_data <- clean_missing_data(downsamp_data,
                                pupil = ps,
                                trial_threshold = .75,
                                subject_trial_threshold =.75)

int_data <- interpolate_data(data = downsamp_data,
                             pupil = ps,
                             type = 'linear')
int_data <- int_data %>%
  mutate(time = Timebin / 100 * 2) %>%
  select(-Timebin)



ggplot(avg_data, aes(x=time, y=avg_ps, color = stim)) + 
  geom_line()+ 
  theme_classic()

ggplot(int_data, aes(x=time, y=ps, color = stim)) + 
  geom_smooth()+ 
  theme_classic()

base_data <- cs_onset %>%
  filter(time>=-500 & time<=0) %>%
  group_by(trial) %>%
  summarise(baseline = mean(ps, na.rm =TRUE))

processed_data <- int_data %>%
  left_join(base_data, by = "trial") %>%
  mutate(ps_cor = ps-baseline)
  
avg_data <- processed_data %>%
  group_by(subject, time, stim) %>%
  summarize(avg_ps= mean(ps, na.rm = TRUE))
ggplot(avg_data, aes(x=time, y=avg_ps, color = stim)) + 
  geom_line()+ 
  theme_classic()
