%% VR single-subject analysis %%
% created by CAS: 28 July 2025
% updated: 6 Jan 2026 
%
% PURPOSE:
%   For seq = 5 only:
%     - extract the 5 lifts immediately BEFORE the switch
%     - preserve trial order WITHIN each block (1→5)
%     - preserve block order ACROSS repetitions (block 1→5)
%     - flatten into a single 1×25 vector
%
% OUTPUT:
%   out.pre5_all_row: 1×25 vector (mcom)
%   out.pre5_all_row_gf: 1×25 vector (gf)
%   out.pre5_blocks: 5×5 (mcom)
%   out.pre5_blocks_gf: 5×5 (gf)
%   out.subj_num
function out = pre_5_row_single(summary_file_name)
    %% --- Load data ---
    if ~contains(summary_file_name,'.txt')
        summary_file_name = [summary_file_name '.txt'];
    end
    
    T = readtable(summary_file_name);
    
    com_condition = T{:,2};
    mcom          = T{:,32};
    gf            = T{:,27};
    
    ntrials = height(T);
    
    %% Subject info
    subj_num = str2double(summary_file_name(2:3));
    
    % Minority side definition
    if mod(subj_num,2)==1
        minority_side = 2;
    else
        minority_side = 1;
    end
    
    %% Find switch
    minority_idx = find(com_condition == minority_side);
    minority_idx = [0; minority_idx];
    seq_len = diff(minority_idx) - 1;
    
    %% collect ALL 5-trial pre-switch blocks
    pre5_blocks    = [];
    pre5_blocks_gf = [];
    
    for i = 1:length(seq_len)
    
        if seq_len(i) == 5
            switch_idx = minority_idx(i+1);
            pre_idx = (switch_idx-5):(switch_idx-1);
    
            if any(pre_idx <= 0) || any(pre_idx > ntrials)
                continue
            end
    
            pre5_blocks(end+1,:)    = mcom(pre_idx)';
            pre5_blocks_gf(end+1,:) = gf(pre_idx)';
        end
    end
    
    %% Flatten
    pre5_all_row    = reshape(pre5_blocks', 1, []);
    pre5_all_row_gf = reshape(pre5_blocks_gf', 1, []);
    
    %% Output
    out.pre5_blocks    = pre5_blocks;
    out.pre5_blocks_gf = pre5_blocks_gf;
    
    out.pre5_all_row    = pre5_all_row;
    out.pre5_all_row_gf = pre5_all_row_gf;
    
    out.subj_num = subj_num;
end