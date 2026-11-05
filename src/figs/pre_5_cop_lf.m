%% VR single-subject analysis %%
% created by CAS: 14 Jan 2026
%
% PURPOSE:
% For seq = 5 only:
%   - extract the 5 lifts immediately BEFORE the switch
%   - get COP and LF and graph is in a scatterplot of LF on the Y and COP
%   on the X

function out = pre_5_cop_lf(summary_file_name)
%% --- Load data ---
if ~contains(summary_file_name,'.txt')
    summary_file_name = [summary_file_name '.txt'];
end

T = readtable(summary_file_name);

com_condition = T{:,2};
cop_diff     = T{:,24};
lf_diff      = T{:,34};
ntrials      = height(T);

%% --- Subject info ---
subj_num = str2double(summary_file_name(2:3));

% Minority side definition
if mod(subj_num,2)==1
    minority_side = 2;
else
    minority_side = 1;
end

%% --- Find switch structure ---
minority_idx = find(com_condition == minority_side);
minority_idx = [0; minority_idx];
seq_len = diff(minority_idx) - 1;

%% --- Collect ALL 5→1 blocks ---
cop_blocks = [];   % rows = blocks, cols = trials (1→5)
lf_blocks  = [];

for i = 1:length(seq_len)

    if seq_len(i) == 5
        switch_idx = minority_idx(i+1);
        pre_idx = (switch_idx-5):(switch_idx-1);

        % Safety check
        if any(pre_idx <= 0) || any(pre_idx > ntrials)
            continue
        end

        % Preserve order
        cop_blocks(end+1,:) = cop_diff(pre_idx)';
        lf_blocks(end+1,:)  = lf_diff(pre_idx)';
    end
end

%% --- Flatten to 25×1 (block1→block5, trial1→trial5) ---
cop_25 = reshape(cop_blocks', [], 1);
lf_25  = reshape(lf_blocks',  [], 1);

%% --- Output ---
out.cop_diff = cop_25;
out.lf_diff  = lf_25;
out.subj_num = subj_num;
fprintf('S%02d: cop_25 size = %d, lf_25 size = %d\n', subj_num, length(cop_25), length(lf_25));
end