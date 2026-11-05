%% VR group scatter analysis with unique participant colors %%
% created by CAS 14 Jan 2026

%% Find all summary files 
files = dir('*summary_*.txt');
file_names = {files.name};
num_subs = length(file_names);

%% Containers 
all_cop = [];
all_lf  = [];
cop_clean = [];
lf_clean = [];
subject_id = [];

%% Generate unique colors
colors = lines(num_subs);  % lines colormap, one color per subject

%% Scatter plot setup
figure; hold on;

%% Loop through subjects 
for i = 1:num_subs

    fname = file_names{i};
    subj_num = str2double(fname(2:3));

    data = pre_5_cop_lf(fname);

    % last 15 trials (comment out for all 25)
    % data.cop_diff = data.cop_diff(end-14:end);
    % data.lf_diff  = data.lf_diff(end-14:end);

    % last three trials of every block
    idx = [3 4 5 8 9 10 13 14 15 18 19 20 23 24 25];
    data.cop_diff = data.cop_diff(idx);
    data.lf_diff  = data.lf_diff(idx);

    % Flip sign for R→L starters
    if mod(subj_num,2) == 0
        data.cop_diff = -data.cop_diff;
        data.lf_diff = -data.lf_diff;
    end

    %% Exclude extreme values 
    valid_idx = data.cop_diff ~= -999 & data.lf_diff ~= -999; % removing -999 if present

    cop_clean = data.cop_diff(valid_idx);
    lf_clean  = data.lf_diff(valid_idx);

    %% Append for output
    all_cop = [all_cop; cop_clean];
    all_lf  = [all_lf; lf_clean];
    subject_id = [subject_id; repmat(subj_num, length(cop_clean), 1)];

    %% Plot subject points with unique color
    %scatter(cop_clean, lf_clean, 40, colors(i,:),
    %'filled','MarkerFaceAlpha', 0.6); %colors & change to 'o' instead of
    %filled for vr
    %% for graphical abstract
    %scatter(cop_clean, lf_clean, 40, [0.494 0.184 0.556], 'filled','MarkerFaceAlpha', 0.6); %purple
    scatter(cop_clean, lf_clean, 40, [0.929 0.694 0.125], 'filled', 'MarkerFaceAlpha', 0.6); %yellow
    %% Axis formatting 
    xlabel('DIGIT POSITIONING (mm)');
    ylabel('LF DIFFERENCE (N)');

    xlim([-50 50]);
    ylim([-10 15]);

    % Dotted vertical line at COP & LF = 0
    xline(0, ':k', 'LineWidth', 2.5);
    yline(0, ':k', 'LineWidth', 2.5);

    ax = gca;
    ax.Box = 'off';
    ax.LineWidth = 1.5;
    ax.FontSize = 16;
    ax.XAxisLocation = 'bottom';
    ax.YAxisLocation = 'left';

    grid on;
    set(gcf,'color','w');
end

%% Mixed-effects regression

% Make sure the variables are column vectors
all_cop = all_cop(:);
all_lf = all_lf(:);
subject_id = subject_id(:);

% Remove invalid observations
good = ~isnan(all_cop) & ~isnan(all_lf) & ...
       ~isinf(all_cop) & ~isinf(all_lf);

model_cop = all_cop(good);
model_lf = all_lf(good);
model_subject = subject_id(good);

% Create table
tbl = table(model_lf, model_cop, categorical(model_subject), ...
    'VariableNames', {'LFdiff','COPdiff','Subject'});

%% Fit mixed-effects model 
lme_model = fitlme(tbl, 'LFdiff ~ COPdiff + (1|Subject)');

%% Generate fitted regression line

x_fit = linspace(min(all_cop), max(all_cop), 100)';

% Fixed-effect estimates
beta = fixedEffects(lme_model);
beta_cov = lme_model.CoefficientCovariance;

% Design matrix for fixed effects
X = [ones(length(x_fit),1), x_fit];

% Population-level fitted values
y_fit = X * beta;

% Standard error of fitted mean
y_se = sqrt(diag(X * beta_cov * X'));

% 95% confidence interval
t_crit = tinv(0.975, lme_model.DFE);
y_lower = y_fit - t_crit .* y_se;
y_upper = y_fit + t_crit .* y_se;

%% Plot 95% confidence interval

fill([x_fit; flipud(x_fit)], ...
    [y_lower; flipud(y_upper)], ...
    [0.929 0.694 0.125], ...
    'FaceAlpha', 0.20, ...
    'EdgeColor', 'none');

%% Plot fitted regression line
plot(x_fit, y_fit, 'k-', 'LineWidth', 2);
