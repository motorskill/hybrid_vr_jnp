%% VR group data analysis %%
% created by CAS: 28 July 2025
% gets data by looping through vr_single_analysis

function out = vr_group_analysis()
% Load all summary .txt files in the directory
all_summary_files = struct2cell(dir('*summary_*.txt'));
file_names = all_summary_files(1, :);
num_files = length(file_names);

all_seq_nums = [1,5];
num_seqs = length(all_seq_nums);

subject_ids = {};
subject_seq1 = containers.Map;
subject_seq5 = containers.Map;

% Data storage
mean_mcoms = struct;
combined_mcoms = [];
combined_seqs = [];
all_mcom_seq_1 = [];
all_mcom_seq_5 = [];
results_table = []; % for easier data extraction
mean_table =[]; 

% Loop through each summary file
for i = 1:num_files
    file_name = file_names{i};
    subject_id = file_name(1:3);  % e.g., 'S01_1'
    subject_num = str2double(file_name(2:3));  % e.g., 01

    current_data = vr_single_analysis(file_name);

    % multiplying by -1 for EVEN subjects
    if mod(subject_num, 2) == 0
        current_data.mcom_seq_1 = -current_data.mcom_seq_1;
        current_data.mcom_seq_5 = -current_data.mcom_seq_5;
    end
    
        % Add to subject-wise list
    if ~ismember(subject_id, subject_ids)
        subject_ids{end+1} = subject_id;
        subject_seq1(subject_id) = current_data.mcom_seq_1;
        subject_seq5(subject_id) = current_data.mcom_seq_5;
    else
        % Append values to existing entry
        subject_seq1(subject_id) = [subject_seq1(subject_id), current_data.mcom_seq_1];
        subject_seq5(subject_id) = [subject_seq5(subject_id), current_data.mcom_seq_5];
    end
    
    % Store data
    mean_mcoms.(subject_id) = current_data.mean_mcom_data;
    all_mcom_seq_1.(subject_id)  = current_data.mcom_seq_1; 
    all_mcom_seq_5.(subject_id)  = current_data.mcom_seq_5;        
    combined_mcoms = [combined_mcoms, current_data.mean_mcom_data];
    combined_seqs = [combined_seqs, current_data.corr_mean_com_data];
    % Create row for this subject
    new_row = table(...
        string(subject_id), ...
        current_data.mcom_seq_1, ...
        current_data.mcom_seq_5, ...
        'VariableNames', {'SubjectID', 'MeanSeq1', 'MeanSeq5'});

    % Append to summary table
    results_table = [results_table; new_row];
end

% Create summary table by averaging across runs
mean_table = table();
for i = 1:length(subject_ids)
    sid = subject_ids{i};
    seq1_mean = mean(subject_seq1(sid));
    seq5_mean = mean(subject_seq5(sid));

    new_row = table(string(sid), seq1_mean, seq5_mean, ...
        'VariableNames', {'SubjectID', 'MeanSeq1', 'MeanSeq5'});
    mean_table = [mean_table; new_row];
end
%% Extra: seq 1 v seq 5 mcoms

all_mcom_seq_1 = [];
all_mcom_seq_5 = [];

for i = 1:num_files
    file_name = file_names{i};
    subject_num = str2double(file_name(2:3));
    current_data = vr_single_analysis(file_name);

    % Flip for even-numbered subjects (R→L)
    if mod(subject_num, 2) == 0
        current_data.mcom_seq_1 = current_data.mcom_seq_1 * -1;
        current_data.mcom_seq_5 = current_data.mcom_seq_5 * -1;
    end

    all_mcom_seq_1 = [all_mcom_seq_1; current_data.mcom_seq_1];
    all_mcom_seq_5 = [all_mcom_seq_5; current_data.mcom_seq_5];
end

x = 1:5;

figure;
hold on;
plot(x, mean(all_mcom_seq_1, 1), '-o', 'DisplayName', 'First 5 Trials');
plot(x, mean(all_mcom_seq_5, 1), '-o', 'DisplayName', 'Last 5 Trials');
xlabel('Post-Switch Block (1–5)', 'FontSize', 20);
ylabel('Mcom at Lift Onset (Nmm)', 'FontSize', 20);
title('Group Average Post-Switch (LR & RL*-1)', 'FontSize', 22);
ylim([-150 150]);
xlim([1 5]);
legend('Location', 'best');
grid on;
set(gca, 'FontSize', 18);
set(gcf, 'color', 'w');

end
