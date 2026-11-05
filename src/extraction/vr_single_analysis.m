%% This code takes one summary file and analyzes it %%
% created by CAS: 28 July 2025
% you can either say summaryFileName = 'FILL' or
% it will be pulled into the vr group analysis script and get data like
% that

function out = vr_single_analysis(summary_file_name, save_fig_yn)
% Append .txt to filename if needed
if ~contains(summary_file_name, '.txt')
    summary_file_name = [summary_file_name '.txt'];
end

summary_file = importdata(summary_file_name);

% Read file as table
summary_table = readtable(summary_file_name); % puts into column and sees that there is a header

% Get subject number
subj_num_str = summary_file_name(2:3); % e.g., 01
subj_num = str2double(subj_num_str);

% Determine side order: odd = L->R, even = R->L
if mod(subj_num, 2) == 1
    exp_side_order = 'left_right';
    majority_side = 1;
    minority_side = 2;
else
    exp_side_order = 'right_left';
    majority_side = 2;
    minority_side = 1;
end

% Extract relevant data
mcom_at_lift_onset = summary_table{:, 32};  % column 32 mcom, 24 cop
com_condition = summary_table{:, 2};        % column 2
trial_nums = summary_table{:, 1};           % column 1
num_trials = height(summary_table);

% Find trial indices when CoM switched to minority side
minority_indices = find(com_condition == minority_side);
minority_indices = [0; minority_indices];  % prepend zero
seq_diff_vec = diff(minority_indices) - 1;  % subtract 1 for true sequence length
order_of_seqs = seq_diff_vec;
max_seq = max(seq_diff_vec);

% Initialize output variables
mcom_data_vec = [];
mean_mcom_data = [];
corr_seq_mean_data = [];
trial_num_vec = [];
seq_vec = [];

% Loop through sequences
for each_seq = [1,5]
    seq_label = ['seq' num2str(each_seq)];
    seq_occurrences = find(seq_diff_vec == each_seq);
    num_repeats = length(seq_occurrences);

    current_seq_mcom_vals = NaN(num_repeats, 1);
    current_seq_trial_nums = NaN(num_repeats, 1);

    for each_repeat = 1:num_repeats
        current_seq = seq_occurrences(each_repeat);
        trial_idx = minority_indices(current_seq + 1);  % +1 to offset prepended 0
        current_trial_num = trial_nums(trial_idx);
        current_mcom = mcom_at_lift_onset(trial_idx-1); % add -1 to see pre switch

        current_seq_trial_nums(each_repeat) = current_trial_num;
        current_seq_mcom_vals(each_repeat) = current_mcom;
    end

    % Store current sequence data
    trial_num_vec = [trial_num_vec, current_seq_trial_nums'];
    mcom_data_vec = [mcom_data_vec, current_seq_mcom_vals'];
    seq_nums = ones(1, num_repeats) * each_seq;
    seq_vec = [seq_vec, seq_nums];

    mean_seq_mcom = mean(mcom_data_vec);
    mean_mcom_data = [mean_mcom_data, mean_seq_mcom];
    mean_seq_val = mean(seq_nums);
    corr_seq_mean_data = [corr_seq_mean_data, mean_seq_val];

    % Assign to output struct
    out.trial_nums.(seq_label) = current_seq_trial_nums;
    out.seq_mcoms.(seq_label) = current_seq_mcom_vals;
    out.seq_nums.(seq_label) = seq_nums;
end

% Add remaining output fields
out.mcom_data_vec = mcom_data_vec;
out.trial_num_vec = trial_num_vec;
out.seq_vec = seq_vec;
out.mean_mcom_data = mean_mcom_data;
out.corr_mean_com_data = corr_seq_mean_data;
out.order_of_seqs = order_of_seqs;

%%%%%%
% Extract first and last 5 Mcom values and compute means
if length(mcom_data_vec) >= 5
    first5 = mcom_data_vec(1:5);
    last5 = mcom_data_vec(end-4:end);
else
    first5 = NaN(1, 5);
    last5 = NaN(1, 5);
end

% Save only the means as scalar values in output
out.mcom_seq_1 = mean(first5);
out.mcom_seq_5 = mean(last5);
