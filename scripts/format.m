% Read the dataset from a CSV file
data = readtable('/Users/bai/Library/CloudStorage/Box-Box/CAMSLab_Projects/Studies/AlcGen/Data/In_Lab/pupillometry/restructured/trial_split/L006/L006_cond_cs_trial2.csv');

% Create the structure S with required fields
S.data.sample = data.ps; % Assuming 'ps' column has the pupil size values
S.data.smp_timestamp = data.time; % Assuming 'time' column has the timestamp values

% Save the structure to a .mat file
save('pupil_data.mat', 'S');