dat_dir='/Users/bai/Library/CloudStorage/Box-Box/CAMSLab_Projects/Studies/AlcGen/Data/In_Lab/pupillometry/restructured/trial_split/L006/';
fname = 'L006_cond_cs.csv';
samp_rate = 500;

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 6);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["subject", "time", "ps", "stim", "trial"];
opts.VariableTypes = ["double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
% opts = setvaropts(opts, "text", "EmptyFieldRule", "auto");

% Import the data
L006cond = readtable(fullfile([dat_dir, fname]), opts);


%% Reformat

S1.data.sample = L006cond.ps;

save_file_path = fullfile(dat_dir, 'S1.mat'); % Define the full file path for the .mat file
save(save_file_path, 'S1'); % Save the S1 variable to a .mat file

%% find NAs
function na_idx = nacheck(df) 
    nrow = size(df,1);
    na_idx = [];
    for n = 1:nrow
        if contains(df(n,:), 'NA')
            na_idx = cat(1,na_idx, n);
        end
    end
end

