%% VR group analysis %%
% created by CAS
% updated: 6 Jan 2026 


%% Find all summary files
files = dir('*summary_*.txt');
file_names = {files.name};
num_subs = length(file_names);

%% Containers
all_pre5    = [];% MCom
all_pre5_gf = []; % GF
subject_ids = [];

%% loop through subjects
for i = 1:num_subs

    fname = file_names{i};
    subj_num = str2double(fname(2:3));

    data = pre_5_row_single(fname);

    % Flip sign for R→L starters (Mcom only)
    if mod(subj_num,2) == 0
        data.pre5_all_row = -data.pre5_all_row;
        data.pre5_all_row_gf = -data.pre5_all_row_gf;  % IMPORTANT: keep alignment !
    end

    all_pre5(end+1,:)    = data.pre5_all_row;
    all_pre5_gf(end+1,:) = data.pre5_all_row_gf;
    subject_ids(end+1)   = subj_num;

end

%% group means
out.group_mean_25    = mean(all_pre5, 1, 'omitnan');
out.group_mean_25_gf = mean(all_pre5_gf, 1, 'omitnan'); 

out.all_pre5    = all_pre5;
out.all_pre5_gf = all_pre5_gf; 

out.subject_ids = subject_ids;

%% block structure (5×5 per subject > group)
nBlocks = 5;
nTrials = 5;

pre5_reshaped = reshape(all_pre5, size(all_pre5,1), nTrials, nBlocks);
pre5_reshaped = permute(pre5_reshaped, [1 3 2]); % subj × block × trial

pre5_reshaped_gf = reshape(all_pre5_gf, size(all_pre5_gf,1), nTrials, nBlocks);
pre5_reshaped_gf = permute(pre5_reshaped_gf, [1 3 2]);

%% --- Block means ---
out.block_means    = squeeze(mean(pre5_reshaped,3)); % MCOM subj × block
out.block_means_gf = squeeze(mean(pre5_reshaped_gf,3)); % GF subj × block

out.group_block_mean    = mean(out.block_means,1,'omitnan');
out.group_block_mean_gf = mean(out.block_means_gf,1,'omitnan');

