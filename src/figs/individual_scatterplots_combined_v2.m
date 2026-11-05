%% VR group + 15 individual scatter plots in one figure %%
% created by CAS, 15 Jan 2026
%
% PURPOSE:
%   - Plot all pre-switch seq=5 COP vs LF points for the group
%   - Overlay individual-subject plots in the same figure (subplots)
%   - Show regression line per subject
%   - Compute r and Fisher z per subject

%% --- Find all summary files ---
files = dir('*summary_*.txt');
file_names = {files.name};
num_subs = length(file_names);

if num_subs > 15
    warning('More than 15 subjects detected; only first 15 will be plotted.');
    file_names = file_names(1:15);
    num_subs = 15;
end

%% Colors for subjects
colors = lines(num_subs); % unique color per subject

%% Participant-specific axis limits

x_limits = {
[-5 10] % Participant 1
[0 40] % Participant 2
[0 20] % Participant 3
[-5 50] % Participant 4
[-40 0] % Participant 5
[-20 45] % Participant 6
[-15 5] % Participant 7
[-10 20] % Participant 8
[0 40] % Participant 9
[-5 50] % Participant 10
[0 30] % Participant 11
[0 50] % Participant 12
[-5 10] % Participant 13
[-5 10] % Participant 14
[-10 30] % Participant 15
};

y_limits = {
[0 10] % Participant 1
[-4 4] % Participant 2
[0 6] % Participant 3
[-4 6] % Participant 4
[-6 15] % Participant 5
[-10 6] % Participant 6
[0 10] % Participant 7
[-6 4] % Participant 8
[-8 1] % Participant 9
[-2 6] % Participant 10
[-6 6] % Participant 11
[-3 5] % Participant 12
[0 8] % Participant 13
[0 15] % Participant 14
[-10 15] % Participant 15
};
%% Containers 
all_cop = [];
all_lf  = [];
subject_id = [];

subject_r = zeros(num_subs,1);
subject_z = zeros(num_subs,1);

%% One figure for all subjects 
figure('Name','All Subjects Scatter');
nCols = 5;
nRows = ceil(num_subs/nCols);
for i = 1:num_subs
    fname = file_names{i};
    subj_num = str2double(fname(2:3));
    data = pre_5_cop_lf(fname);

    % Keep only last 3 trials of every block
    idx = [3:5, 8:10, 13:15, 18:20, 23:25];
    data.cop_diff = data.cop_diff(idx);
    data.lf_diff  = data.lf_diff(idx);

    % Flip for R->L starters
    if mod(subj_num,2)==0
        data.cop_diff = -data.cop_diff;
        data.lf_diff  = -data.lf_diff;
    end

    % Exclude extremes
    valid_idx = data.cop_diff ~= -999 & data.lf_diff ~= -999; % remove -999 if present
    cop_clean = data.cop_diff(valid_idx);
    lf_clean  = data.lf_diff(valid_idx);

    % Store for group analysis
    all_cop = [all_cop; cop_clean];
    all_lf  = [all_lf; lf_clean];
    subject_id = [subject_id; repmat(subj_num,length(cop_clean),1)];

    % Compute correlation r per subject
    if length(cop_clean) >= 3
        r = corr(cop_clean, lf_clean);
        z = atanh(r);
    else
        r = NaN;
        z = NaN;
    end

    subject_r(i) = r;
    subject_z(i) = z;

    % Subplot for this subject
    subplot(nRows, nCols, i); hold on;
    scatter(cop_clean, lf_clean, 40, colors(i,:), ...
        'filled', 'MarkerFaceAlpha',0.6);

    % Regression line
    if length(cop_clean) >= 2
        p = polyfit(cop_clean, lf_clean, 1);
        x_fit = linspace(min(cop_clean), max(cop_clean), 100);
        y_fit = polyval(p, x_fit);
        plot(x_fit, y_fit, 'k--', 'LineWidth', 1.2);
    end

    xline(0, ':k', 'LineWidth',1);

    % Participant-specific axis limits
    xlim(x_limits{subj_num});
    ylim(y_limits{subj_num});

    % xlabel('COP DIFF (mm)');
    % ylabel('\it{LF DIFFERENCE (N)}');
    % title(sprintf('S%02d: r=%.2f', subj_num, r));
    grid on;
    set(gca,'FontSize',10);
end

%% Output
out.all_cop = all_cop;
out.all_lf  = all_lf;
out.subject_id = subject_id;
out.subject_r = subject_r;
out.subject_z = subject_z;
