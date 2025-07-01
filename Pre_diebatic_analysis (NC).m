%% Whole body PET reveals inter-organ metabolic connectivity disruption across glycemic states

%%  Run GLM for each organ with FBG_Category
organ_cols = {'liver','adrenalGlandLeft_norm', 'adrenalGlandRight_norm', ...
              'brain_norm', 'colon_norm','gluteusMaximusLeft_norm',	'gluteusMaximusRight_norm',...
          	'heartMyocardium_norm',	'iliopsoasLeft_norm',	'iliopsoasRight_norm',	'kidneyLeft_norm',...
        	'kidneyRight_norm', 'lungLeft_norm'	,'lungRight_norm',...
         	'pancreas_norm',	'smallBowel_norm','spleen_norm', 'subcutaneousFat_norm',...
        	'thyroidLeft_norm',	'thyroidRight_norm','visceralFat_norm'};


SUVTable = readtable('/Volumes/Extreme Pro/SZNormal/LookupTable/Mean_SUL_1000.xlsx');

SUVTable.FBG_Category = categorical(SUVTable.FBG_Category);
SUVTable.FBG_Category = reordercats(SUVTable.FBG_Category, ["Normal", "Pre-diabetic", "Diabetic"]);

models = cell(length(organ_cols), 1);

p_values = table('Size', [length(organ_cols), 3], ...
    'VariableTypes', {'string', 'double', 'double'}, ...
    'VariableNames', {'Organ', 'p_pre_vs_normal', 'p_diab_vs_normal'});

for i = 1:length(organ_cols)
    organ = organ_cols{i};
    formula = sprintf('%s ~ FBG_Category + Age + Gender + BMI', organ);
    models{i} = fitglm(SUVTable, formula);
    
    % Extract coefficient table
    coeffs = models{i}.Coefficients;
    Rsq(i) = models{i}.Rsquared.Ordinary;% extract R square

    % Store organ name
    p_values.Organ(i) = string(organ);

    % Find and store p-values for each comparison
    idx_pre = strcmp(coeffs.Row, 'FBG_Category_Pre-diabetic');
    idx_diab = strcmp(coeffs.Row, 'FBG_Category_Diabetic');

    p_values.p_pre_vs_normal(i) = coeffs.pValue(idx_pre);
    p_values.p_diab_vs_normal(i) = coeffs.pValue(idx_diab);
end

% Apply FDR Correction (Benjamini-Hochberg)
p_values.p_pre_vs_normal_FDR = bh_fdr(p_values.p_pre_vs_normal);
p_values.p_diab_vs_normal_FDR = bh_fdr(p_values.p_diab_vs_normal);

output_folder = '/Volumes/Extreme Pro/SZNormal/Pre_diab/p_values';
filenames = 'Mean_SUL_1000.xlsx';
path = fullfile(output_folder,filenames);
% writetable(p_values,path);


%% Group metabolic network construction : LAYOUT 
% Figure layout
organ_cols = {'liver','adrenalGlandLeft_norm', 'adrenalGlandRight_norm','thyroidLeft_norm', ...
              'brain_norm', ...
          	'heartMyocardium_norm',	...
            'gluteusMaximusLeft_norm',	'gluteusMaximusRight_norm','iliopsoasLeft_norm',	'iliopsoasRight_norm',...
        	'kidneyLeft_norm','kidneyRight_norm',...
            'lungLeft_norm'	,'lungRight_norm',...
         	'pancreas_norm',	'smallBowel_norm','colon_norm','spleen_norm',...
            'subcutaneousFat_norm',...
        	'visceralFat_norm'};
% Change display name 
roi_display_names = containers.Map;
roi_display_names('brain_norm') = 'Brain';
roi_display_names('pancreas_norm') = 'Pancreas';
roi_display_names('spleen_norm') = 'Spleen';
roi_display_names('liver') = 'Liver';
roi_display_names('heartMyocardium_norm') = 'Heart';
roi_display_names('lungLeft_norm') = 'Lung L';
roi_display_names('lungRight_norm') = 'Lung R';
roi_display_names('kidneyLeft_norm') = 'Kidney L';
roi_display_names('kidneyRight_norm') = 'Kidney R';
roi_display_names('adrenalGlandLeft_norm') = 'Adrenal L';
roi_display_names('adrenalGlandRight_norm') = 'Adrenal R';
roi_display_names('colon_norm') = 'Colon';
roi_display_names('gluteusMaximusLeft_norm') = 'Gluteus Maximus L';
roi_display_names('gluteusMaximusRight_norm') = 'Gluteus Maximus R';
roi_display_names('iliopsoasLeft_norm') = 'Iliopsoas L';
roi_display_names('iliopsoasRight_norm') = 'Iliopsoas R';
roi_display_names('smallBowel_norm') = 'Small Bowel';
roi_display_names('subcutaneousFat_norm') = 'Subcutaneous Fat';
roi_display_names('thyroidLeft_norm') = 'Thyroid';
roi_display_names('visceralFat_norm') = 'Visceral Fat';


clean_labels = cell(size(organ_cols));
for i = 1:length(organ_cols)
    if isKey(roi_display_names, organ_cols{i})
        clean_labels{i} = roi_display_names(organ_cols{i});
    else
        clean_labels{i} = organ_cols{i};  
    end
end

% change display color 
organ_groups = containers.Map;
organ_groups('brain_norm') = 'Brain';
organ_groups('liver') = 'Digestive';
organ_groups('pancreas_norm') = 'Digestive';
organ_groups('smallBowel_norm') = 'Digestive';
organ_groups('colon_norm') = 'Digestive';
organ_groups('spleen_norm') = 'Immune';
organ_groups('heartMyocardium_norm') = 'Cardiac';
organ_groups('lungLeft_norm') = 'Respiratory';
organ_groups('lungRight_norm') = 'Respiratory';
organ_groups('kidneyLeft_norm') = 'Renal';
organ_groups('kidneyRight_norm') = 'Renal';
organ_groups('adrenalGlandLeft_norm') = 'Endocrine';
organ_groups('adrenalGlandRight_norm') = 'Endocrine';
organ_groups('thyroidLeft_norm') = 'Endocrine';
% organ_groups('thyroidRight_norm') = 'Endocrine';

organ_groups('iliopsoasLeft_norm') = 'Muscle';
organ_groups('iliopsoasRight_norm') = 'Muscle';
organ_groups('gluteusMaximusLeft_norm') = 'Muscle';
organ_groups('gluteusMaximusRight_norm') = 'Muscle';
organ_groups('visceralFat_norm') = 'Fat';
organ_groups('subcutaneousFat_norm') = 'Fat';

group_labels = values(organ_groups, organ_cols); % Assign a color per group
group_labels = categorical(group_labels);
group_colors = lines(numel(categories(group_labels)));  % distinct colors
node_colors = group_colors(double(group_labels), :);


%% 7b.1 : Create covariate-adjusted SUL + FDR 

SUV_table = readtable('/Volumes/Extreme Pro/SZNormal/LookupTable/Mean_SUL_1000.xlsx');

organ_cols = {'liver','adrenalGlandLeft_norm', 'adrenalGlandRight_norm','thyroidLeft_norm', ...
              'brain_norm', ...
          	'heartMyocardium_norm',	...
            'gluteusMaximusLeft_norm',	'gluteusMaximusRight_norm','iliopsoasLeft_norm',	'iliopsoasRight_norm',...
        	'kidneyLeft_norm','kidneyRight_norm',...
            'lungLeft_norm'	,'lungRight_norm',...
         	'pancreas_norm',	'smallBowel_norm','colon_norm','spleen_norm',...
            'subcutaneousFat_norm',...
        	'visceralFat_norm'};
% Create residual table 
residuals = table();
residuals.Patient = SUV_table.Patient;

for i = 1:numel(organ_cols)
    roi = organ_cols{i};
    mdl = fitlm(SUV_table, sprintf('%s ~ Age + Gender + BMI', roi));  % Create regression model
    residuals.(roi) = mdl.Residuals.Raw;  % Store residuals
end  % Extracts the raw residuals, which represent the part of SUV not explained by covariates.
     % These residuals are treated as the "adjusted metabolic signal".

% Stratify by group 
group_label = SUV_table.FBG_Category;

res_norm = residuals(strcmp(group_label, 'Normal'), organ_cols);
res_predm = residuals(strcmp(group_label, 'Pre-diabetic'), organ_cols);
res_dm = residuals(strcmp(group_label, 'Diabetic'), organ_cols);


% Define bootstrap parameters
n_boot = 1000;
n_rois = numel(organ_cols);


% Apply to each group
[R_norm,z_norm]    = bootstrap_network(table2array(res_norm), n_boot, 0.05);
[R_predm,z_predm]  = bootstrap_network(table2array(res_predm), n_boot, 0.05);
[R_dm,z_dm]        = bootstrap_network(table2array(res_dm), n_boot, 0.05);


%% Center B
SUV_table = readtable('/Volumes/Extreme Pro/SZNormal/LookupTable/Mean_SUL_54.xlsx');

% Create residual table 
residuals = table();
residuals.Patient = SUV_table.Patient;
for i = 1:numel(organ_cols)
    roi = organ_cols{i};
    mdl = fitlm(SUV_table, sprintf('%s ~ Age + Gender + BMI', roi));  % Create regression model
    residuals.(roi) = mdl.Residuals.Raw;  % Store residuals
end  % Extracts the raw residuals, which represent the part of SUV not explained by covariates.
     % These residuals are treated as the "adjusted metabolic signal".

% Stratify by group 
group_label = SUV_table.FBG_Category;
res_norm = residuals(strcmp(group_label, 'Normal'), organ_cols);
res_predm = residuals(strcmp(group_label, 'Pre-diabetic'), organ_cols);
res_dm = residuals(strcmp(group_label, 'Diabetic'), organ_cols);

% Define bootstrap parameters
n_boot = 1000;
n_rois = numel(organ_cols);

% Apply to each group
[R_norm_b,z_norm_b]    = bootstrap_network(table2array(res_norm), n_boot, 0.05);
[R_predm_b,z_predm_b]  = bootstrap_network(table2array(res_predm), n_boot, 0.05);
[R_dm_b,z_dm_b]        = bootstrap_network(table2array(res_dm), n_boot, 0.05);

%% Mantel test
% ========1=======cross center 
% comparing norm networks between Center A and Center B
[r_mantel_norm, p_norm] = mantel_test(R_norm, R_norm_b, 100000);

% You can repeat this for other pairs:
[r_mantel_predm, p_predm] = mantel_test(R_predm, R_predm_b, 100000);

[r_mantel_dm, p_dm] = mantel_test(R_dm, R_dm_b, 100000);

%%  ploting group network
close all

R_all = {R_norm, R_predm, R_dm};

group_names = {'Normal', 'Pre-diabetic', 'Diabetic'};

for g = 1:3
    R = R_all{g};
    R(isnan(R)) = 0; 
    % R(1:size(R,1)+1:end) = 0;  % zero diagonal
    R(abs(R) < 0.5) = 0;      % Zero out weak correlations with |r| < 0.5
  
    % Build graph
    R_abs = abs(R);
    R_abs_sym = (R_abs + R_abs') / 2;
   
    fig = figure; 
    circularGraph(R_abs_sym,'Label',clean_labels,'ColorMap',node_colors);
 
    % fig_names = {'normal055centerb.fig','pre_055centerb.fig','diabetic_055centerb.fig'};
    % file_base = fig_names{g}; 
    % 
    % figure_path = fullfile('/Volumes/Extreme Pro/SZNormal/Pre_diab/network_plot', fig_names{g});
    % saveas(fig, figure_path);
    % saveas(fig, fullfile('/Volumes/Extreme Pro/SZNormal/Pre_diab/network_plot', [file_base '.png'])); 
end

%% ploting heatmap
close all
R_all = {R_norm, R_predm, R_dm};

group_names = {'Normal', 'Pre-diabetic', 'Diabetic'};

for g = 1:3
    R = R_all{g};
    R(isnan(R)) = 0; 
    R(abs(R) < 0) = NaN;      % NaN out weak correlations with |r| < 0.5

    figure('Color', 'w');
    h = heatmap(clean_labels, clean_labels, R, ...
    'Colormap', bluewhitered(), ...
    'ColorLimits', [-1 1], ...
    'MissingDataColor', [0.95 0.95 0.95], ...
    'MissingDataLabel', '');

    % Aesthetics
    h.CellLabelColor = 'none';  % hide numbers
    h.FontName = 'Arial';
    h.FontSize = 10;
    h.GridVisible = 'off';

end


%% topology
SUV_table = readtable('/Volumes/Extreme Pro/SZNormal/LookupTable/Mean_SUL_1000.xlsx');

organ_cols =  {'liver','adrenalGlandLeft', 'adrenalGlandRight','thyroidLeft', ...
              'brain', ...
          	'heartMyocardium',	...
            'gluteusMaximusLeft',	'gluteusMaximusRight','iliopsoasLeft',	'iliopsoasRight',...
        	'kidneyLeft','kidneyRight',...
            'lungLeft'	,'lungRight',...
         	'pancreas',	'smallBowel','colon','spleen',...
            'subcutaneousFat',...
        	'visceralFat'};


% Create residual table 
residuals = table();
residuals.Patient = SUV_table.Patient;
for i = 1:numel(organ_cols)
    roi = organ_cols{i};
    mdl = fitlm(SUV_table, sprintf('%s ~ Age + Gender + BMI', roi));  % Create regression model
    residuals.(roi) = mdl.Residuals.Raw;  % Store residuals
end  

% -------------用到建立网络时存的residual变量----------------------------------------------------------
data_norm = table2array(residuals(strcmp(group_label, 'Normal'), organ_cols));
data_predm = table2array(residuals(strcmp(group_label, 'Pre-diabetic'), organ_cols));
data_dm = table2array(residuals(strcmp(group_label, 'Diabetic'), organ_cols));

% -------------得到每一次boot，对应的网络， 存在R_boot_xx 和A_boot_xx 里----------------------------------------------------------
% Bootstrap parameters
n_boot = 1000;
n_roi = numel(organ_cols);
alpha = 0.05;
% Magnitude
r_thresh = 0.3; % THRESHOLD

% Initialize output
R_boot_norm   = NaN(n_roi, n_roi, n_boot);
A_boot_norm   = false(n_roi, n_roi, n_boot);
R_boot_predm  = NaN(n_roi, n_roi, n_boot);
A_boot_predm  = false(n_roi, n_roi, n_boot);
R_boot_dm     = NaN(n_roi, n_roi, n_boot);
A_boot_dm     = false(n_roi, n_roi, n_boot);

% Bootstrap function
bootstrap_group = @(data, R_out, A_out) deal( ...
    arrayfun(@(b) ...
        compute_one_bootstrap(data, r_thresh), ...
        1:n_boot, 'UniformOutput', false) ...
);

% Bootstrap loop — NORMAL
for b = 1:n_boot
    idx = randsample(size(data_norm,1), size(data_norm,1), true);
    Xb = data_norm(idx,:);
     [R, P] = corr(Xb, 'type', 'Pearson', 'rows', 'pairwise');
    R(isnan(R)) = 0; P(isnan(P)) = 1;
    % FDR correction
    upper_idx = find(triu(true(n_rois), 1));
    raw_p = P(upper_idx);
    [~, ~, ~, fdr_adj] = fdr_bh(raw_p, alpha, 'pdep', 'no');
    % Binary mask for significant + strong correlations
    sig = false(n_rois);
    sig(upper_idx) = (fdr_adj < alpha) & (abs(R(upper_idx)) > r_thresh);
    sig = sig | sig';  % Symmetrize

    R_boot_norm(:,:,b) = R;
    A_boot_norm(:,:,b) = sig;
end

% Bootstrap loop — PRE-DIABETIC
for b = 1:n_boot
    idx = randsample(size(data_predm,1), size(data_predm,1), true);
    Xb = data_predm(idx,:);
    [R, P] = corr(Xb, 'type', 'Pearson', 'rows', 'pairwise');
    R(isnan(R)) = 0; P(isnan(P)) = 1;
    % FDR correction
    upper_idx = find(triu(true(n_rois), 1));
    raw_p = P(upper_idx);
    [~, ~, ~, fdr_adj] = fdr_bh(raw_p, alpha, 'pdep', 'no');
    % Binary mask for significant + strong correlations
    sig = false(n_rois);
    sig(upper_idx) = (fdr_adj < alpha) & (abs(R(upper_idx)) > r_thresh);
    sig = sig | sig';  % Symmetrize

    R_boot_predm(:,:,b) = R;
    A_boot_predm(:,:,b) = sig;
end

% Bootstrap loop — DIABETIC
for b = 1:n_boot
    idx = randsample(size(data_dm,1), size(data_dm,1), true);
    Xb = data_dm(idx,:);
    [R, P] = corr(Xb, 'type', 'Pearson', 'rows', 'pairwise');
    R(isnan(R)) = 0; P(isnan(P)) = 1;
    % FDR correction
    upper_idx = find(triu(true(n_rois), 1));
    raw_p = P(upper_idx);
    [~, ~, ~, fdr_adj] = fdr_bh(raw_p, alpha, 'pdep', 'no');
    % Binary mask for significant + strong correlations
    sig = false(n_rois);
    sig(upper_idx) = (fdr_adj < alpha) & (abs(R(upper_idx)) > r_thresh);
    sig = sig | sig';  % Symmetrize

    R_boot_dm(:,:,b) = R;
    A_boot_dm(:,:,b) = sig;
end

% 8b -------------得到每一次boot中，对应网络的average metrics信息， 存在struct里----------------------------------------------------------
% Preallocate
n_boot = size(R_boot_norm, 3);

% Initialize struct array with expected fields
empty_struct = struct( ...
    'Density', [], ...
    'AvgDegree', [], ...
    'AvgStrength', [], ...
    'AvgClustering', []);

metrics_N   = repmat(empty_struct, n_boot, 1);
metrics_Pre = repmat(empty_struct, n_boot, 1);
metrics_DM  = repmat(empty_struct, n_boot, 1);


for b = 1:n_boot
    % === NORMAL group
    R = R_boot_norm(:,:,b);
    A = A_boot_norm(:,:,b);
    [metrics_N(b), ~, ~] = compute_network_metrics(R, A);

    % === PRE-DIABETIC group
    R = R_boot_predm(:,:,b);
    A = A_boot_predm(:,:,b);
    [metrics_Pre(b), ~, ~] = compute_network_metrics(R, A);

    % === DIABETIC group
    R = R_boot_dm(:,:,b);
    A = A_boot_dm(:,:,b);
    [metrics_DM(b), ~, ~] = compute_network_metrics(R, A);
end

density_N     = [metrics_N.Density]';
Degree_N   = [metrics_N.AvgDegree]';
Strength_N = [metrics_N.AvgStrength]';
Cluster_N  = [metrics_N.AvgClustering]';

density_Pre     = [metrics_Pre.Density]';
Degree_Pre   = [metrics_Pre.AvgDegree]';
Strength_Pre = [metrics_Pre.AvgStrength]';
Cluster_Pre  = [metrics_Pre.AvgClustering]';

density_DM     = [metrics_DM.Density]';
Degree_DM   = [metrics_DM.AvgDegree]';
Strength_DM = [metrics_DM.AvgStrength]';
Cluster_DM  = [metrics_DM.AvgClustering]';

% Compute means for each group and each metric
mean_density = [mean(density_N), mean(density_Pre), mean(density_DM)];
mean_degree  = [mean(Degree_N),  mean(Degree_Pre),  mean(Degree_DM)];
mean_strength = [mean(Strength_N), mean(Strength_Pre), mean(Strength_DM)];
mean_cluster  = [mean(Cluster_N),  mean(Cluster_Pre),  mean(Cluster_DM)];

% Create a table
T = table(mean_density', mean_degree', mean_strength', mean_cluster', ...
    'VariableNames', {'Density', 'AvgDegree', 'AvgStrength', 'AvgClustering'}, ...
    'RowNames', {'Normoglycemic', 'Pre-diabetic', 'Diabetic'});

% Display the table
disp(T);


% -------计算pairwise的比较，得到p 值 -----------------------------------------------
% -----Anova 测三组的差别--------------------------------

metrics = {'Density', 'Degree', 'Strength', 'Clustering'};
group_names = {'Normal', 'Pre-diabetic', 'Diabetic'};
pair_names = {'Normal vs Pre-diabetic', 'Normal vs Diabetic', 'Pre-diabetic vs Diabetic'};

% Define metric data
metric_data = {
    [density_N; density_Pre; density_DM], ...
    [Degree_N; Degree_Pre; Degree_DM], ...
    [Strength_N; Strength_Pre; Strength_DM], ...
    [Cluster_N; Cluster_Pre; Cluster_DM]
};

% Prepare group labels for all
n = length(density_N);  % assume same N for each group
group_labels = [ ...
    repmat({'Normal'}, n, 1);
    repmat({'Pre-diabetic'}, n, 1);
    repmat({'Diabetic'}, n, 1)];

% Initialize result table
Metric = {};
Group_Comparison = {};
ANOVA_p = {};
PostHoc_p = {};

% Loop through each metric to Run ANOVA
for i = 1:length(metrics)
    % Run ANOVA
    [p_anova, ~, stats] = anova1(metric_data{i}, group_labels, 'off');

    % Run post hoc test
    c = multcompare(stats, 'CType', 'bonferroni', 'Display', 'off');

    % Store all pairwise comparisons
    for j = 1:3  % 3 pairwise comparisons
        Metric{end+1,1} = metrics{i};
        Group_Comparison{end+1,1} = pair_names{j};

        % Format p-values
        p_val = c(j,6);
        if p_val < 0.001
            p_str = '< 0.001';
        else
            p_str = sprintf('%.3f', p_val);
        end
        PostHoc_p{end+1,1} = p_str;

        % Format ANOVA p
        if p_anova < 0.001
            anova_str = '< 0.001';
        else
            anova_str = sprintf('%.3f', p_anova);
        end
        ANOVA_p{end+1,1} = anova_str;
    end
end

% Create final table
StatsTable = table(Metric, Group_Comparison, ANOVA_p, PostHoc_p);
disp(StatsTable);

%% use previous saved variable to create line plot (用之前的,不能跑第二次）
close all
% R = 0.3
mean_density_03;
mean_degree_03;
mean_strength_03;
mean_cluster_03;
% R =0.5
mean_density_05;
mean_degree_05;
mean_strength_05;
mean_cluster_05;
% R = 0.7
mean_density_07;
mean_strength_07;
mean_cluster_07;

x = 1:3;

% --- Plot Network Density ---
figure;
plot(x, mean_density_03, '-o', 'LineWidth', 2); hold on;
plot(x, mean_density_05, '-s', 'LineWidth', 2);
plot(x, mean_density_07, '-^', 'LineWidth', 2);
% title('Network Density Across Glycemic States');
xticks([1,2,3])
xticklabels({'Normoglycemic', 'Pre-Diabetic', 'Diabetic'})
ylabel('Density');
legend('r > 0.3', 'r > 0.5', 'r > 0.7', 'Location', 'northeast');
grid on;

% --- Plot Avg Degree ---
figure;
plot(x, mean_degree_03, '-o', 'LineWidth', 2); hold on;
plot(x, mean_degree_05, '-s', 'LineWidth', 2);
plot(x, mean_degree_07, '-s', 'LineWidth', 2);
plot(x, NaN(size(x)), 'k-^');  % Placeholder for legend consistency
title('Average Node Degree Across Glycemic States');
xticks([1,2,3]); xticklabels({'Normoglycemic', 'Pre-Diabetic', 'Diabetic'});
ylabel('Avg Degree');
legend('r > 0.3', 'r > 0.5', 'N/A for r > 0.7', 'Location', 'northeast'); % if not available
grid on;

% --- Plot Avg Strength ---
figure;
plot(x, mean_strength_03, '-o', 'LineWidth', 2); hold on;
plot(x, mean_strength_05, '-s', 'LineWidth', 2);
plot(x, mean_strength_07, '-^', 'LineWidth', 2);
title('Average Strength Across Glycemic States');
xticks([1,2,3]); xticklabels({'Normoglycemic', 'Pre-Diabetic', 'Diabetic'});
ylabel('Avg Strength');
legend('r > 0.3', 'r > 0.5', 'r > 0.7', 'Location', 'northeast');
grid on;

% --- Plot Avg Clustering Coefficient ---
figure;
plot(x, mean_cluster_03, '-o', 'LineWidth', 2); hold on;
plot(x, mean_cluster_05, '-s', 'LineWidth', 2);
plot(x, mean_cluster_07, '-^', 'LineWidth', 2);
title('Average Clustering Coefficient Across Glycemic States');
xticks([1,2,3]); xticklabels({'Normoglycemic', 'Pre-Diabetic', 'Diabetic'});
ylabel('Avg Clustering');
legend('r > 0.3', 'r > 0.5', 'r > 0.7', 'Location', 'northeast');
grid on;


%% 8b:---node level analysis (没用）--得到每一次boot中，对应每个ROI的 metrics信息 ----------------------

[strength_N, cluster_N]     = compute_node_metrics_boot(R_boot_norm, A_boot_norm);
[strength_Pre, cluster_Pre] = compute_node_metrics_boot(R_boot_predm, A_boot_predm);
[strength_DM, cluster_DM]   = compute_node_metrics_boot(R_boot_dm, A_boot_dm);


% === Healthy Controls ===
mean_strength_norm = mean(strength_N, 1);
mean_clustering_norm = mean(cluster_N, 1);
%  Organize into table ===
organ_hub_table_norm = table(organ_cols', mean_strength_norm', mean_clustering_norm', ...
    'VariableNames', {'Organ', 'MeanStrength', 'MeanClustering'});
%  Rank by strength or clustering ===
organ_hub_table_norm = sortrows(organ_hub_table_norm, 'MeanStrength', 'descend');  % or 'MeanClustering'
% Display
disp(organ_hub_table_norm);

% === Pre-DM ===
mean_strength_Pre = mean(strength_Pre, 1);
mean_clustering_Pre = mean(cluster_Pre, 1);
%  Organize into table ===
organ_hub_table_Pre = table(organ_cols', mean_strength_Pre', mean_clustering_Pre', ...
    'VariableNames', {'Organ', 'MeanStrength', 'MeanClustering'});
%  Rank by strength or clustering ===
organ_hub_table_Pre = sortrows(organ_hub_table_Pre, 'MeanStrength', 'descend');  % or 'MeanClustering'
% Display
disp(organ_hub_table_Pre);

% === DM ===
mean_strength_DM = mean(strength_DM, 1);
mean_clustering_DM = mean(cluster_DM, 1);
%  Organize into table ===
organ_hub_table_DM = table(organ_cols', mean_strength_DM', mean_clustering_DM', ...
    'VariableNames', {'Organ', 'MeanStrength', 'MeanClustering'});
%  Rank by strength or clustering ===
organ_hub_table_DM = sortrows(organ_hub_table_DM, 'MeanStrength', 'descend');  % or 'MeanClustering'
% Display
disp(organ_hub_table_DM);

%% 9a step:  计算差异网络

close all

thresh= 0;

R_all = {z_norm, z_predm, z_dm};

group_pairs = {[1, 3], [1, 2]};  % Normal vs Diabetic, Normal vs Pre-diabetic
titles = { '\Delta Network: Diabetic − Normal', '\Delta Network: Pre-diabetic − Normal' };

for k = 1:2
    i1 = group_pairs{k}(1);  % Index for Normal
    i2 = group_pairs{k}(2);  % Index for Diabetic or Pre-diabetic
    
    R1 = R_all{i1};  % baseline group (Normal)
    R1(isnan(R1)) = 0; 
    R2 = R_all{i2};  % comparison group
    R2(isnan(R2)) = 0; 

    % Step 1: Compute Δr
    R_delta_filt = R2 - R1;
    inv_z  = @(z) (exp(2*z) - 1) ./ (exp(2*z) + 1);     % Inverse Fisher
    R_diff = inv_z(R_delta_filt); % inverse Fisher transform

    % Set threshold
    R_diff(abs(R_diff) < thresh) = 0;      % Zero out weak correlations with |r| < 0.5
   
    % Step 2:Create masks
    R_pos = R_diff;
    R_pos(R_pos <= 0) = 0;
    
    R_neg = R_diff;
    R_neg(R_neg >= 0) = 0;

    abs_delta_filt_p = abs(R_pos);
    abs_delta_filt_n = abs(R_neg);

    % Step 3: Symmetrize the matrix
    R_delta_sym_p = (abs_delta_filt_p + abs_delta_filt_p') / 2;
    R_delta_sym_n = (abs_delta_filt_n + abs_delta_filt_n') / 2;
 
    % Color
    node_colors_red = repmat([1 0 0], size(R_pos, 1), 1);  % red RGB
    node_colors_blue = repmat([0 0 1], size(R_neg, 1), 1);  % blue RGB

    % Step 5: Plot circular graph
    fig = figure;
    circularGraph(R_delta_sym_p, 'Label', clean_labels, 'ColorMap', node_colors_red);
    hold on;
    circularGraph(R_delta_sym_n, 'Label', clean_labels, 'ColorMap', node_colors_blue);
    fig_names = {'diabetic_norm.fig','pre_norm.fig'};
    % figure_path = fullfile('/Volumes/Extreme Pro/SZNormal/Pre_diab/network_plot/delta_network_plot', fig_names{k});
    % saveas(fig, figure_path);

end

%% 9b step:  差异网络用柱状图表示
close all

% R_all = {R_norm, R_predm, R_dm};

R_all = {z_norm, z_predm, z_dm};
thresh = 0;
group_pairs = {[1, 2], [1, 3]};  

% Define Nature-style colors
red_color = [178, 24, 43]/255;
blue_color = [33, 102, 172]/255;

for k = 1:2
    i1 = group_pairs{k}(1);  % Norm
    i2 = group_pairs{k}(2);  % Pre-Diabetic or Diabetic
    
    R1 = R_all{i1}; R1(isnan(R1)) = 0;
    R2 = R_all{i2}; R2(isnan(R2)) = 0;

    % Δr calculation
    % R_diff = abs(R2) - abs(R1);
    R_delta_filt = abs(R2) - abs(R1);
    inv_z  = @(z) (exp(2*z) - 1) ./ (exp(2*z) + 1);     % Inverse Fisher
    R_diff = inv_z(R_delta_filt); % inverse Fisher transform



    R_diff(abs(R_diff) < thresh) = 0;
    R_delta_sym = (R_diff + R_diff') / 2;

    % Extract upper triangle
    n = length(clean_labels);
    [rows, cols] = find(triu(true(n), 1));
    ROI_1 = clean_labels(rows);
    ROI_2 = clean_labels(cols);
    delta_r_values = R_delta_sym(sub2ind([n, n], rows, cols));
    PairLabels = strcat(ROI_1, "–", ROI_2);

    % Table + sort
    DeltaTable = table(PairLabels', delta_r_values, 'VariableNames', {'Pair', 'Delta_r'});
    [~, sort_idx] = sort(abs(DeltaTable.Delta_r), 'descend');
    TopN = 20;
    TopTable = DeltaTable(sort_idx(1:TopN), :);

    % Bar plot setup
    delta_r = TopTable.Delta_r;
    num_bars = numel(delta_r);
    bar_colors = zeros(num_bars, 3);
    bar_colors(delta_r > 0, :) = repmat(red_color, sum(delta_r > 0), 1);
    bar_colors(delta_r < 0, :) = repmat(blue_color, sum(delta_r < 0), 1);

    % Create figure
    figure('Color', 'w', 'Position', [100 100 1000 400]);
    b = bar(delta_r, 'FaceColor', 'flat');
    b.CData = bar_colors;

    % Styling
    xticks(1:num_bars);
    xticklabels(TopTable.Pair);
    xtickangle(60);
    ylabel('\Deltar (Correlation Change)', 'FontName', 'Helvetica', 'FontSize', 10);
    ylim([-0.6 0.6]);
    line([0 num_bars+1], [0 0], 'Color', 'k', 'LineWidth', 0.75); % baseline

    % Axes formatting
    ax = gca;
    ax.FontName = 'Helvetica';
    ax.FontSize = 8;
    ax.XColor = 'k';
    ax.YColor = 'k';
    ax.TickDir = 'out';
    ax.Box = 'off';
    ax.LineWidth = 1;
    ax.YGrid = 'on';
    ax.GridLineStyle = '--';

    % Optional: title or export label
    % title('Top 20 Inter-Organ \Deltar Changes', 'FontWeight', 'bold');

end
%% 10a step: Individualized network (leave one out test)
% 按照年龄组，性别，BMI组分别建立normative network
% 其中 refNET的 FBG 要都是正常的。 
close all

organ_cols = {'liver','adrenalGlandLeft_norm', 'adrenalGlandRight_norm', ...
              'brain_norm', 'colon_norm','gluteusMaximusLeft_norm',	'gluteusMaximusRight_norm',...
          	'heartMyocardium_norm',	'iliopsoasLeft_norm',	'iliopsoasRight_norm',	'kidneyLeft_norm',...
        	'kidneyRight_norm', 'lungLeft_norm'	,'lungRight_norm',...
         	'pancreas_norm',	'smallBowel_norm','spleen_norm', 'subcutaneousFat_norm',...
        	'thyroidLeft_norm',	'thyroidRight_norm','visceralFat_norm'};

FullData_healthy = readtable('/Volumes/Extreme Pro/SZNormal/LookupTable/Mean_SUL_1000_OW&OB.xlsx');
FullData_healthy_normFBG = FullData_healthy(strcmp(FullData_healthy.FBG_Category,'Normal'),:);

% Groups 
Age_Group = {'20-39', '40-59', '60-79'};
Gender = {'F', 'M'};
BMI_Category = {'Normal', 'OW&Obese'};

% Set network parameter
min_population = 20; 
mag_threshold = 0;



% Loop through combinations
for a = 1:length(Age_Group)
    for g = 1:length(Gender)
        for b = 1:length(BMI_Category)
          
            % Subset the data
            subset = FullData_healthy_normFBG( ...
                strcmp(FullData_healthy_normFBG.Age_Group, Age_Group{a}) & ...
                strcmp(FullData_healthy_normFBG.Gender, Gender{g}) & ...
                strcmp(FullData_healthy_normFBG.BMI_Category, BMI_Category{b}), :);

            % Skip if not enough samples
            if height(subset) < min_population 
                fprintf('Skipped: %s, %s, %s (n=%d)\n', ...
                    Age_Group{a}, Gender{g}, BMI_Category{b}, height(subset));
                continue;
            end

            % Extract the SUL values for the 21 ROIs
            data = subset(:, organ_cols);
            X = table2array(data);
  
            % Full correlation matrix
            R_full = corr(X, 'Rows', 'pairwise', 'Type', 'Pearson');
            upper_idx = find(triu(true(size(R_full)), 1));
            R_full_vec = R_full(upper_idx);

            % Leave-One-Out validation
            n_sub = size(X,1);
            corrs = zeros(n_sub,1);
            for i = 1:n_sub
                X_loo = X([1:i-1, i+1:end], :);
                R_loo = corr(X_loo, 'Rows', 'pairwise', 'Type', 'Pearson');
                R_loo_vec = R_loo(upper_idx);
                corrs(i) = corr(R_full_vec, R_loo_vec, 'Type', 'Spearman', 'Rows','complete');
            end

            % Report result
            fprintf('Group: %s_%s_%s\n', Age_Group{a}, Gender{g}, BMI_Category{b});
            fprintf('Mean Spearman Correlation with Full Network: %.3f (SD = %.3f)\n', mean(corrs), std(corrs));
            fprintf('Group number: %.1f\n',   height(subset));
       
        end
    end
end

%% 10a: normative group network: 一个一个的图分开展示

FullData_healthy = readtable('/Volumes/Extreme Pro/SZNormal/LookupTable/Mean_SUL_1000_OW&OB.xlsx');
FullData_healthy_normFBG = FullData_healthy(strcmp(FullData_healthy.FBG_Category,'Normal'),:);

% Groups 
Age_Group = {'20-39', '40-59', '60-79'};
Gender = {'F', 'M'};
BMI_Category = {'Normal', 'OW&Obese'};

% Parameters
min_population = 20; 
mag_threshold = 0.5;
output_dir = '/Volumes/Extreme Pro/SZNormal/Pre_diab/SubgroupMax_0.7';
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

% Loop through combinations
for a = 1:length(Age_Group)
    for g = 1:length(Gender)
        for b = 1:length(BMI_Category)

            % Subset the data
            subset = FullData_healthy_normFBG( ...
                strcmp(FullData_healthy_normFBG.Age_Group, Age_Group{a}) & ...
                strcmp(FullData_healthy_normFBG.Gender, Gender{g}) & ...
                strcmp(FullData_healthy_normFBG.BMI_Category, BMI_Category{b}), :);

            % Skip if not enough samples
            if height(subset) < min_population 
                fprintf('Skipped: %s, %s, %s (n=%d)\n', ...
                    Age_Group{a}, Gender{g}, BMI_Category{b}, height(subset));
                continue;
            end

            % Extract ROI data
            data = subset(:, organ_cols);
            X = table2array(data);

            % Correlation and FDR
            [R, P] = corr(X, 'Rows', 'pairwise', 'Type', 'Pearson');
            r_thresh = mag_threshold;
            n = size(R, 1);  
            upper_mask = triu(true(n), 1);

            pvals = P(upper_mask);
            [h_norm, ~, ~] = fdr_bh(pvals, 0.05, 'pdep', 'no');
            rvals = R(upper_mask);

            sig_mask = (h_norm == 1) & (abs(rvals) > r_thresh);
            adjacency = false(n);
            adjacency(upper_mask) = sig_mask;
            adjacency = adjacency | adjacency';

            % Store matrix
            group_id = sprintf('%s_%s_%s', Age_Group{a}, Gender{g}, BMI_Category{b});
            fprintf('Group: %s (n=%d)\n', group_id, height(subset));
            save(fullfile(output_dir, sprintf('CorrMatrix_%s.mat', group_id)), ...
                'R', 'adjacency', 'group_id');

            % Plot in separate figure
            A = adjacency;
            R_filt = R .* A;
            R_abs = abs(R_filt);
            R_abs_sym = (R_abs + R_abs') / 2;

            fig = figure('Name', group_id, 'NumberTitle', 'off');
            circularGraph(R_abs_sym, 'Label', clean_labels, 'ColorMap', node_colors);

            % Save figure
            saveas(fig, fullfile(output_dir, [group_id '.fig']));
            

        end
    end
end


%% 10b Step : Get individualized connection values
% set parameter
mag_threshold = 0;
min_population = 20;

FullData_healthy = readtable('/Volumes/Extreme Pro/SZNormal/LookupTable/Mean_SUL_1000_OW&OB.xlsx');
FullData_healthy_normFBG = FullData_healthy(strcmp(FullData_healthy.FBG_Category,'Normal'),:);
FullData_healthy_PreFBG = FullData_healthy(strcmp(FullData_healthy.FBG_Category,'Diabetic'),:);
filename = 'DM_indiv.xlsx';


out_dir = '/Volumes/Extreme Pro/SZNormal/Pre_diab/Individualized';  
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

% Define table
n = length(organ_cols);
[roi_i, roi_j] = find(triu(ones(n), 1));  % upper triangle indices
pair_names = cell(numel(roi_i), 1);
for k = 1:numel(roi_i)
    pair_names{k} = sprintf('%s_%s', organ_cols{roi_i(k)}, organ_cols{roi_j(k)});
end

% Initialize full variable names
var_names = ['PatientID', pair_names'];
empty_data = cell(height(FullData_healthy_PreFBG), numel(var_names));
all_deltas = cell2table(empty_data, 'VariableNames', var_names);

% Fill PatientID column
all_deltas.PatientID = FullData_healthy_PreFBG.Patient;

% Loop through each pre-diabetic subject
% for i = 1:height(FullData_healthy_PreFBG)
for i = 5
    subj = FullData_healthy_PreFBG(i, :);

    % Step 1: Find normative group match
    match_norm = FullData_healthy_normFBG(...
        strcmp(FullData_healthy_normFBG.Age_Group, subj.Age_Group) & ...
        strcmp(FullData_healthy_normFBG.Gender, subj.Gender) & ...
        strcmp(FullData_healthy_normFBG.BMI_Category, subj.BMI_Category), :);

    if height(match_norm)< min_population
        warning(['Not enough normative matchs for subject: ', FullData_healthy_PreFBG.Patient{i}]);
        continue
    end

    % Step 2: Merge subject with matched normative group
    merged_group = [match_norm; subj];

    % Step 3: Compute correlation and adjacency on merged
    merged_vals = table2array(merged_group(:, organ_cols));
    [R_merged, P_merged] = corr(merged_vals, 'Rows', 'pairwise', 'Type', 'Pearson');
    % FDR correction
    r_thresh = mag_threshold;  % Threshold for correlation magnitude
    n = size(R_merged, 1);  
    upper_mask = triu(true(n), 1);
    pvals = P_merged(upper_mask);
    [h_pre, ~, ~] = fdr_bh(pvals, 0.05, 'pdep', 'no'); 
    rvals = R_merged(upper_mask);
    sig_mask = (h_pre == 1) & (abs(rvals) > r_thresh);
    adjacency = false(n);
    adjacency(upper_mask) = sig_mask;
    adj_merged = adjacency | adjacency';

    % Step 4: Compute normative baseline
    norm_vals = table2array(match_norm(:, organ_cols));
    [R_norm, P_norm] = corr(norm_vals, 'Rows', 'pairwise', 'Type', 'Pearson');
    pvals = P_norm(upper_mask);% FDR correction
    [h_norm, ~, ~] = fdr_bh(pvals, 0.05, 'pdep', 'no'); % Benjamini–Hochberg function
    rvals = R_norm(upper_mask);
    sig_mask = (h_norm == 1) & (abs(rvals) > r_thresh);
    adjacency = false(n);
    adjacency(upper_mask) = sig_mask;
    adj_norm = adjacency | adjacency';

    % Step 5: Δ delta
    R_delta = R_merged.* adj_merged - R_norm.*adj_norm;

    % Step 6: Save to table 
    r_vals = zeros(numel(roi_i), 1);
    for k = 1:numel(roi_i)
        r_vals(k) = R_delta(roi_i(k), roi_j(k));
    end

    % Add to table
    all_deltas{i, 2:end} = num2cell(r_vals');
end


all_deltas.Properties.VariableNames = ['Patient', pair_names'];
empty_rows = all(cellfun(@isempty, table2cell(all_deltas(:,2:end))), 2);
all_deltas(empty_rows, :) = [];

% Define filename
% writetable(all_deltas, fullfile(out_dir, filename));

%% 计算个人之间的差异性
pre_dm_deltas = readtable('/Volumes/Extreme Pro/SZNormal/Pre_diab/Individualized/Pre_indiv.xlsx');

delta_values = table2array(pre_dm_deltas(:, 2:end));
roi_pairs = pre_dm_deltas.Properties.VariableNames(2:end);
% Compute variance per ROI pair
variances = var(delta_values, 0, 1, 'omitnan');

% Sort by variance (descending)
[sorted_var, sort_idx] = sort(variances, 'descend');
top_roi_pairs = roi_pairs(sort_idx(1:10));  % top 10
top_values = sorted_var(1:10);

% % Plot histogram of variances
% figure;
% histogram(variances, 30);
% xlabel('Variance of ΔR');
% ylabel('Frequency');
% title('Distribution of ROI-pair Variance Across Individuals');

% Plot boxplot
figure;
boxplot(variances, 'Notch', 'on');
ylabel('Variance');
title('Boxplot of ROI-pair Variability');

% Display top variable ROI pairs
fprintf('Top 10 high-variance ROI pairs:\n');
for i = 1:10
    fprintf('%2d. %-30s  Var = %.4f\n', i, top_roi_pairs{i}, top_values(i));
end

%% Variance 画图 （nature版本）
% Load data
pre_dm_deltas = readtable('/Volumes/Extreme Pro/SZNormal/Pre_diab/Individualized/Pre_indiv.xlsx');
dm_deltas = readtable('/Volumes/Extreme Pro/SZNormal/Pre_diab/Individualized/DM_indiv.xlsx');

% Extract numeric data
delta_pre = table2array(pre_dm_deltas(:, 2:end));
delta_dm = table2array(dm_deltas(:, 2:end));

% Compute variance across individuals for each ROI pair
var_pre = var(delta_pre, 0, 1, 'omitnan');
var_dm = var(delta_dm, 0, 1, 'omitnan');

% Combine data for boxplot
all_var = [var_pre, var_dm]';
group_labels = [repmat({'Pre-diabetic'}, length(var_pre), 1); repmat({'Diabetic'}, length(var_dm), 1)];

% Nature-style colors (colorblind-friendly)
pre_color = [33 102 172] / 255;    % Blue
dm_color  = [178 24 43] / 255;     % Red
box_colors = [pre_color; dm_color];

% Create figure
figure('Color', 'w', 'Units', 'normalized', 'Position', [0.3, 0.3, 0.35, 0.5]);

% Draw boxplot
boxplot(all_var, group_labels, 'Notch', 'on', 'Widths', 0.5, 'Symbol', '');
hold on;

% Overlay scatter with jitter
jitter_strength = 0.08;
x_pre = 1 + (rand(size(var_pre)) - 0.5) * jitter_strength;
x_dm = 2 + (rand(size(var_dm)) - 0.5) * jitter_strength;
% scatter(x_pre, var_pre, 10, 'k', 'filled', 'MarkerFaceAlpha', 0.3);
% scatter(x_dm, var_dm, 10, 'k', 'filled', 'MarkerFaceAlpha', 0.3);


% Color boxes manually
h = findobj(gca, 'Tag', 'Box');
for j = 1:length(h)
    patch(get(h(j), 'XData'), get(h(j), 'YData'), box_colors(j, :), ...
        'FaceAlpha', 0.6, 'EdgeColor', 'k', 'LineWidth', 1.2);
end

% Find and modify outlier handles
outliers = findobj(gca, 'Tag', 'Outliers');
for i = 1:length(outliers)
    outliers(i).Marker = 'o';        
    outliers(i).MarkerEdgeColor = [0 0 0];  % black
    outliers(i).MarkerSize = 5;
    outliers(i).LineWidth = 0.9;
end


% Style axes
set(gca, ...
    'FontName', 'Helvetica', ...
    'FontSize', 9, ...
    'LineWidth', 1, ...
    'TickDir', 'out', ...
    'Box', 'off', ...
    'YGrid', 'on', ...
    'GridLineStyle', '--', ...
    'XTickLabel', {'Pre-diabetic', 'Diabetic'}, ...
    'XTickLabelRotation', 0);

ylabel('Variance of \Deltar', 'FontSize', 10, 'FontName', 'Helvetica');

%% 11a Step: 把pre和DM结合在同一个table里--merge
% Load data (replace with actual filename)
FullData_healthy = readtable('/Volumes/Extreme Pro/SZNormal/LookupTable/Mean_SUL_1000_OW&OB.xlsx');
FullData_healthy_DMFBG = FullData_healthy(strcmp(FullData_healthy.FBG_Category,'Diabetic'),:);
T = readtable('/Volumes/Extreme Pro/SZNormal/Pre_diab/Individualized/DM_indiv.xlsx');  
merged_DM = innerjoin(FullData_healthy_DMFBG, T, 'Keys', 'Patient');

FullData_healthy_PreFBG = FullData_healthy(strcmp(FullData_healthy.FBG_Category,'Pre-diabetic'),:);
T_pre = readtable('/Volumes/Extreme Pro/SZNormal/Pre_diab/Individualized/Pre_indiv.xlsx');  
merged_pre = innerjoin(FullData_healthy_PreFBG, T_pre, 'Keys', 'Patient');

CombinedTable = [merged_DM; merged_pre];
roi_cols = pair_names;  % ROI-pair columns

delta_data = CombinedTable(:,roi_cols);
group_labels = CombinedTable.FBG_Category;

writetable(CombinedTable, '/Volumes/Extreme Pro/SZNormal/Pre_diab/Individualized/CombinedTable.xlsx');

%% 计算北大和河南merge的inter-organ 数据
% Groups 
% set parameter
mag_threshold = 0;
min_population = 20;

organ_cols = {'liver','adrenalGlandLeft_norm', 'adrenalGlandRight_norm','thyroidLeft_norm', ...
              'brain_norm', ...
          	'heartMyocardium_norm',	...
            'gluteusMaximusLeft_norm',	'gluteusMaximusRight_norm','iliopsoasLeft_norm',	'iliopsoasRight_norm',...
        	'kidneyLeft_norm','kidneyRight_norm',...
            'lungLeft_norm'	,'lungRight_norm',...
         	'pancreas_norm',	'smallBowel_norm','colon_norm','spleen_norm',...
            'subcutaneousFat_norm',...
        	'visceralFat_norm'};

FullData_healthy = readtable('/Volumes/Extreme Pro/SZNormal/LookupTable/Mean_SUL_1000_OW&OB.xlsx');
FullData_healthy_normFBG = FullData_healthy(strcmp(FullData_healthy.FBG_Category,'Normal'),:);

% load 北大医院数据
beida = readtable('/Volumes/Extreme Pro/SZNormal/LookupTable/北大和河南的数据merge.xlsx');

FullData_healthy_PreFBG = beida(strcmp(beida.FBG_Category,'Diabetic'),:);
filename = 'dm_indiv.xlsx';


out_dir = '/Volumes/Extreme Pro/SZNormal/Pre_diab/Individualized_beida&henan';  
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

% Define table
n = length(organ_cols);
[roi_i, roi_j] = find(triu(ones(n), 1));  % upper triangle indices
pair_names = cell(numel(roi_i), 1);
for k = 1:numel(roi_i)
    pair_names{k} = sprintf('%s_%s', organ_cols{roi_i(k)}, organ_cols{roi_j(k)});
end

% Initialize full variable names
var_names = ['PatientID', pair_names'];
empty_data = cell(height(FullData_healthy_PreFBG), numel(var_names));
all_deltas = cell2table(empty_data, 'VariableNames', var_names);

% Fill PatientID column
all_deltas.PatientID = FullData_healthy_PreFBG.Patient;

% Loop through each pre-diabetic subject
for i = 1:height(FullData_healthy_PreFBG)

    subj = FullData_healthy_PreFBG(i, :);

    % Step 1: Find normative group match
    match_norm = FullData_healthy_normFBG(...
        strcmp(FullData_healthy_normFBG.Age_Group, subj.Age_Group) & ...
        strcmp(FullData_healthy_normFBG.Gender, subj.Gender) & ...
        strcmp(FullData_healthy_normFBG.BMI_Category, subj.BMI_Category), :);

    if height(match_norm)< min_population
        warning(['Not enough normative matchs for subject: ', FullData_healthy_PreFBG.Patient{i}]);
        continue
    end

    % Step 2: Merge subject with matched normative group
    match_norm = match_norm(:, organ_cols);
    subj = subj(:, organ_cols);
    merged_group = [match_norm; subj];

    % Step 3: Compute correlation and adjacency on merged
    merged_vals = table2array(merged_group(:, organ_cols));
    [R_merged, P_merged] = corr(merged_vals, 'Rows', 'pairwise', 'Type', 'Pearson');
    % FDR correction
    r_thresh = mag_threshold;  % Threshold for correlation magnitude
    n = size(R_merged, 1);  
    upper_mask = triu(true(n), 1);
    pvals = P_merged(upper_mask);
    [h_pre, ~, ~] = fdr_bh(pvals, 0.05, 'pdep', 'no'); 
    rvals = R_merged(upper_mask);
    sig_mask = (h_pre == 1) & (abs(rvals) > r_thresh);
    adjacency = false(n);
    adjacency(upper_mask) = sig_mask;
    adj_merged = adjacency | adjacency';

    % Step 4: Compute normative baseline
    norm_vals = table2array(match_norm(:, organ_cols));
    [R_norm, P_norm] = corr(norm_vals, 'Rows', 'pairwise', 'Type', 'Pearson');
    pvals = P_norm(upper_mask);% FDR correction
    [h_norm, ~, ~] = fdr_bh(pvals, 0.05, 'pdep', 'no'); % Benjamini–Hochberg function
    rvals = R_norm(upper_mask);
    sig_mask = (h_norm == 1) & (abs(rvals) > r_thresh);
    adjacency = false(n);
    adjacency(upper_mask) = sig_mask;
    adj_norm = adjacency | adjacency';

    % Step 5: Δ delta
    R_delta = R_merged.* adj_merged - R_norm.*adj_norm;

    % Step 6: Save to table 
    r_vals = zeros(numel(roi_i), 1);
    for k = 1:numel(roi_i)
        r_vals(k) = R_delta(roi_i(k), roi_j(k));
    end

    % Add to table
    all_deltas{i, 2:end} = num2cell(r_vals');
end


all_deltas.Properties.VariableNames = ['Patient', pair_names'];
empty_rows = all(cellfun(@isempty, table2cell(all_deltas(:,2:end))), 2);
all_deltas(empty_rows, :) = [];

writetable(all_deltas, fullfile(out_dir, filename));

%% Merge beida into one table 
FullData_healthy = readtable('/Volumes/Extreme Pro/SZNormal/LookupTable/北大和河南的数据merge.xlsx');
FullData_healthy_DMFBG = FullData_healthy(strcmp(FullData_healthy.FBG_Category,'Diabetic'),:);
T = readtable('/Volumes/Extreme Pro/SZNormal/Pre_diab/Individualized_beida&henan/dm_indiv.xlsx');  
merged_DM = innerjoin(FullData_healthy_DMFBG, T, 'Keys', 'Patient');

FullData_healthy_PreFBG = FullData_healthy(strcmp(FullData_healthy.FBG_Category,'Pre-diabetic'),:);
T_pre = readtable('/Volumes/Extreme Pro/SZNormal/Pre_diab/Individualized_beida&henan/pre_indiv.xlsx');  
merged_pre = innerjoin(FullData_healthy_PreFBG, T_pre, 'Keys', 'Patient');

CombinedTable = [merged_DM; merged_pre];
roi_cols = pair_names;  % ROI-pair columns

delta_data = CombinedTable(:,roi_cols);
group_labels = CombinedTable.FBG_Category;

writetable(CombinedTable, '/Volumes/Extreme Pro/SZNormal/Pre_diab/Individualized_beida&henan/CombinedTable.xlsx');

% 北大一共22 个人入选，DM 4人，pre-dm 18人
% 北大河南共 41人，  dm 7人， pre-dm 34人。
%% （不用）Logistic regression with external （using inter-organ data)

ext_AUC = [];
ext_ACC = [];
ext_F1 = [];
ext_MCC = [];

n_models = 5;
for m = 1:n_models
    rng(m); 

% Load data
CombinedTable = readtable('/Volumes/Extreme Pro/SZNormal/Pre_diab/Individualized/CombinedTable.xlsx');
roi_cols = CombinedTable.Properties.VariableNames(74:end);%inter only

% roi_cols = CombinedTable.Properties.VariableNames(44:74); % norm

% roi_cols = CombinedTable.Properties.VariableNames(14:43); % original PET values

% roi_cols = CombinedTable.Properties.VariableNames(46:end); % combined

X = table2array(CombinedTable(:, roi_cols));
Y = categorical(CombinedTable.FBG_Category);

 % Convert labels to binary: Diabetic vs Others
Y_bin = strcmp(string(Y), 'Diabetic');

% Outer CV
    outerCV = cvpartition(Y_bin, 'KFold', 5);

    auc_vals = zeros(outerCV.NumTestSets, 1);
    acc_vals = zeros(outerCV.NumTestSets, 1);
    sens_vals = zeros(outerCV.NumTestSets, 1);
    spec_vals = zeros(outerCV.NumTestSets, 1);

    for i = 1:outerCV.NumTestSets
        trainIdx = training(outerCV, i);
        testIdx = test(outerCV, i);

        X_train = X(trainIdx, :);
        Y_train = Y_bin(trainIdx);
        X_test = X(testIdx, :);
        Y_test = Y_bin(testIdx);

        % Standardize
        mu = mean(X_train, 1);
        sigma = std(X_train, [], 1);
        sigma(sigma==0) = 1;
        X_train = (X_train - mu) ./ sigma;
        X_test = (X_test - mu) ./ sigma;

        % Undersample majority class
        idx_pos = find(Y_train == 1);
        idx_neg = find(Y_train == 0);
        n_min = min(length(idx_pos), length(idx_neg));
        sel_pos = datasample(idx_pos, n_min);
        sel_neg = datasample(idx_neg, n_min);
        sel = [sel_pos; sel_neg];
        X_bal = X_train(sel, :);
        Y_bal = Y_train(sel);

        % remove low variance features
        tol = 0.01;
        std_dev = std(X_bal, 0, 1);
        keep_features = std_dev > tol;
        X_bal = X_bal(:, keep_features);
        X_test = X_test(:, keep_features);

        % Inner CV for feature selection and lambda
        [B, FitInfo] = lassoglm(X_bal, Y_bal, 'binomial', 'CV', 5,'Options',...
            statset('UseParallel', false, 'MaxIter', 1e8));
        idxLambda1SE = FitInfo.Index1SE;
         selected = find(B(:, idxLambda1SE) ~= 0);
        coef_sel = [FitInfo.Intercept(idxLambda1SE); B(selected, idxLambda1SE)];
        X_test_sel = X_test(:, selected);
        probs = glmval(coef_sel, X_test_sel, 'logit');
       
        % Threshold
        preds = probs > 0.5;

        % Confusion matrix
        TP = sum((preds==1)&(Y_test==1));
        TN = sum((preds==0)&(Y_test==0));
        FP = sum((preds==1)&(Y_test==0));
        FN = sum((preds==0)&(Y_test==1));

        acc_vals(i) = (TP + TN) / length(Y_test);
        sens_vals(i) = TP / (TP + FN + eps);
        spec_vals(i) = TN / (TN + FP + eps);
        [~, ~, ~, auc_vals(i)] = perfcurve(Y_test, probs, 1);
    end

    fprintf('AUC: %.3f (95%% CI: %.3f–%.3f)\n', mean(auc_vals), quantile(auc_vals, [0.025 0.975]));
    fprintf('Accuracy: %.2f%%\n', mean(acc_vals)*100);


% === External validation ===
    External = readtable('/Volumes/Extreme Pro/SZNormal/Pre_diab/Individualized_beida&henan/CombinedTable.xlsx');
    X_ext = table2array(External(:, roi_cols));
    Y_ext = strcmp(External.FBG_Category, 'Diabetic');

    % Use final model on all training data
    mu_all = mean(X);
    sigma_all = std(X); sigma_all(sigma_all==0) = 1;
    X_std_all = (X - mu_all) ./ sigma_all;

    % Undersample
    idx_pos = find(Y_bin == 1);
    idx_neg = find(Y_bin == 0);
    n_min = min(length(idx_pos), length(idx_neg));
    sel = [datasample(idx_pos, n_min); datasample(idx_neg, n_min)];
    X_bal_all = X_std_all(sel, :);
    Y_bal_all = Y_bin(sel);

    % Train final logistic model
    [B_all, FitInfo_all] = lassoglm(X_bal_all, Y_bal_all, 'binomial', 'CV', 5);
    idx_final = FitInfo_all.Index1SE;
    coef_final = [FitInfo_all.Intercept(idx_final); B_all(:, idx_final)];
    selected_final = find(B_all(:, idx_final) ~= 0);

    % Predict external
    X_ext_std = (X_ext - mu_all) ./ sigma_all;
    X_ext_matched = X_ext_std(:, keep_features);
    coef_sel_final = [FitInfo_all.Intercept(idx_final); B_all(selected_final, idx_final)];
    X_ext_sel = X_ext_matched(:, selected_final);
    probs_ext = glmval(coef_sel_final, X_ext_sel, 'logit');
    preds_ext = probs_ext > 0.5;

    TP = sum((preds_ext==1)&(Y_ext==1));
    TN = sum((preds_ext==0)&(Y_ext==0));
    FP = sum((preds_ext==1)&(Y_ext==0));
    FN = sum((preds_ext==0)&(Y_ext==1));

    acc = (TP+TN)/length(Y_ext);
    precision = TP / (TP + FP + eps);
    recall = TP / (TP + FN + eps);
    f1 = 2 * precision * recall / (precision + recall + eps);
    mcc = ((TP * TN) - (FP * FN)) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) + eps);
    [~, ~, ~, auc_ext] = perfcurve(Y_ext, probs_ext, 1);

    ext_AUC = [ext_AUC; auc_ext];
    ext_ACC = [ext_ACC; acc];
    ext_F1 = [ext_F1; f1];
    ext_MCC = [ext_MCC; mcc];

    fprintf('External AUC: %.3f | ACC: %.2f%% | F1: %.3f | MCC: %.3f\n', auc_ext, acc*100, f1, mcc);
end

% Summary
fprintf('\nMean External AUC: %.3f\n', mean(ext_AUC));
fprintf('AUC 95%% CI: %.3f–%.3f\n', quantile(ext_AUC, [0.025 0.975]));

%% （不用）random forest with external （using inter-organ data)
% Initialize external evaluation metrics
ext_AUC = []; ext_ACC = []; ext_F1 = []; ext_MCC = [];

n_models = 10;  % number of seeds / runs
for m = 1:n_models
    rng(m);

    % === Load Data ===
    CombinedTable = readtable('/Volumes/Extreme Pro/SZNormal/Pre_diab/Individualized/CombinedTable.xlsx');
    % inter_organ = CombinedTable.Properties.VariableNames(75:end); % inter organ

   
    original= CombinedTable.Properties.VariableNames(14:43);
    roi_cols = [original,'Age','BMI'];
    X = table2array(CombinedTable(:, roi_cols));
    Y = CombinedTable.FBG_Category;
    Y_bin = strcmp(string(Y), 'Diabetic');  % Binary target

    % === Feature Selection: Variance + Correlation ===
    std_thresh = 0.1;
    std_vals = std(X, 0, 1);
    keep_std = std_vals > std_thresh;
    X = X(:, keep_std);  % remove low variance

    corr_matrix = corr(X);
    corr_matrix(logical(eye(size(corr_matrix)))) = 0;
    corr_thresh = 0.95;
    to_remove = any(abs(corr_matrix) > corr_thresh, 2);
    X = X(:, ~to_remove);  % remove highly correlated

    keep_features = keep_std; keep_features(keep_std) = ~to_remove;  % for external match

    % === Outer Cross-validation ===
    outerCV = cvpartition(Y_bin, 'KFold', 5);
    auc_vals = zeros(outerCV.NumTestSets,1);
    acc_vals = zeros(outerCV.NumTestSets,1);
    sens_vals = zeros(outerCV.NumTestSets,1);
    spec_vals = zeros(outerCV.NumTestSets,1);

    for i = 1:outerCV.NumTestSets
        trainIdx = training(outerCV, i);
        testIdx = test(outerCV, i);

        X_train = X(trainIdx, :);
        Y_train = Y_bin(trainIdx);
        X_test = X(testIdx, :);
        Y_test = Y_bin(testIdx);

        % Standardize
        mu = mean(X_train);
        sigma = std(X_train);
        sigma(sigma == 0) = 1;
        X_train = (X_train - mu) ./ sigma;
        X_test = (X_test - mu) ./ sigma;

        % Undersample majority class
        idx_pos = find(Y_train == 1);
        idx_neg = find(Y_train == 0);
        n_min = min(length(idx_pos), length(idx_neg));
        sel_pos = datasample(idx_pos, n_min);
        sel_neg = datasample(idx_neg, n_min);
        sel = [sel_pos; sel_neg];
        X_bal = X_train(sel, :);
        Y_bal = Y_train(sel);

        % Optional: apply SMOTE (if needed)
        [X_bal, Y_bal] = smote(X_bal, [], 5, 'Class', double(Y_bal)+1);

        % Train Random Forest
        RF = TreeBagger(100, X_bal, categorical(Y_bal), ...
            'Method', 'classification', 'OOBPrediction', 'on');

        % Predict
        [~, scores] = predict(RF, X_test);
        preds = scores(:,2) > 0.5;

        TP = sum(preds==1 & Y_test==1);
        TN = sum(preds==0 & Y_test==0);
        FP = sum(preds==1 & Y_test==0);
        FN = sum(preds==0 & Y_test==1);

        acc_vals(i) = (TP+TN)/length(Y_test);
        sens_vals(i) = TP / (TP + FN + eps);
        spec_vals(i) = TN / (TN + FP + eps);
        [~,~,~,auc_vals(i)] = perfcurve(Y_test, scores(:,2), 1);
    end

    % === External Validation ===
    External = readtable('/Volumes/Extreme Pro/SZNormal/Pre_diab/Individualized_beida&henan/CombinedTable_Copy.xlsx');
    X_ext = table2array(External(:, roi_cols));
    X_ext = X_ext(:, keep_features);  % apply same feature selection
    Y_ext = External.FBG_Category;
    Y_ext_bin = strcmp(string(Y_ext), 'Diabetic');

    % Standardize external set using training stats
    mu_all = mean(X);
    sigma_all = std(X); sigma_all(sigma_all==0) = 1;
    X_ext_std = (X_ext - mu_all) ./ sigma_all;

    % Re-train final RF on all data (undersampled)
    idx_pos = find(Y_bin == 1);
    idx_neg = find(Y_bin == 0);
    n_min = min(length(idx_pos), length(idx_neg));
    sel = [datasample(idx_pos, n_min); datasample(idx_neg, n_min)];
    X_train_final = (X(sel,:) - mu_all) ./ sigma_all;
    Y_train_final = Y_bin(sel);

    RF_final = TreeBagger(200, X_train_final, categorical(Y_train_final), ...
        'Method', 'classification');

    % Predict external
    [~, score_ext] = predict(RF_final, X_ext_std);
    pred_ext = score_ext(:,2) > 0.5;

    TP = sum(pred_ext==1 & Y_ext_bin==1);
    TN = sum(pred_ext==0 & Y_ext_bin==0);
    FP = sum(pred_ext==1 & Y_ext_bin==0);
    FN = sum(pred_ext==0 & Y_ext_bin==1);

    ext_acc = (TP+TN)/length(Y_ext_bin);
    precision = TP / (TP + FP + eps);
    recall = TP / (TP + FN + eps);
    f1 = 2 * (precision * recall) / (precision + recall + eps);
    mcc = ((TP*TN)-(FP*FN)) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) + eps);
    [~, ~, ~, ext_auc] = perfcurve(Y_ext_bin, score_ext(:,2), 1);

    ext_AUC = [ext_AUC; ext_auc];
    ext_ACC = [ext_ACC; ext_acc];
    ext_F1 = [ext_F1; f1];
    ext_MCC = [ext_MCC; mcc];

    fprintf('External AUC: %.3f | ACC: %.2f%% | F1: %.3f | MCC: %.3f\n', ext_auc, ext_acc*100, f1, mcc);
end
% === Summary Report ===
fprintf('\n--- External Summary ---\n');
fprintf('Mean External AUC: %.3f (95%% CI: %.3f–%.3f)\n', mean(ext_AUC), quantile(ext_AUC,[0.025 0.975]));
fprintf('Mean Accuracy: %.2f%% (95%% CI: %.2f–%.2f)\n', mean(ext_ACC)*100, quantile(ext_ACC,[0.025 0.975])*100);
fprintf('Mean F1: %.3f | Mean MCC: %.3f\n', mean(ext_F1), mean(ext_MCC));

%% RF without feature selection (使用）
ext_AUC = [];
ext_ACC = [];
ext_SEN = [];
ext_SPE = [];

n_models = 1;
for m = 1:n_models
    rng(m); 

% Load data
CombinedTable = readtable('/Volumes/Extreme Pro/SZNormal/Pre_diab/Individualized/CombinedTable.xlsx');
inter = CombinedTable.Properties.VariableNames(75:end); %inter only

% inter = CombinedTable.Properties.VariableNames(14:43); % original 
% roi_cols = CombinedTable.Properties.VariableNames(44:74); % norm

% roi_cols = CombinedTable.Properties.VariableNames(46:end); % combined

% roi_cols = ['Age','BMI', inter]; 
roi_cols = inter;

X = table2array(CombinedTable(:, roi_cols));
Y = categorical(CombinedTable.FBG_Category);

% Create stratified k-folds
cv = cvpartition(Y, 'KFold', 5);

% Initialize performance metrics
auc_vals = zeros(cv.NumTestSets, 1);
acc_vals = zeros(cv.NumTestSets, 1);
sens_vals = zeros(cv.NumTestSets, 1);
spec_vals = zeros(cv.NumTestSets, 1);

for i = 1:cv.NumTestSets
    % Split training and test data
    trainIdx = training(cv, i);
    testIdx = test(cv, i);

    X_train = X(trainIdx, :);
    Y_train = Y(trainIdx);
    X_test = X(testIdx, :);
    Y_test = Y(testIdx);

    mu = mean(X_train, 1);
    sigma = std(X_train, 0, 1);
    sigma(sigma == 0) = 1;  % avoid division by zero
    X_train_std = (X_train - mu) ./ sigma;
    X_test_std  = (X_test - mu) ./ sigma; 

    % Convert Y to numeric for SMOTE if needed
    [Y_cat, ~, Y_numeric] = unique(Y_train);  

    % === Apply SMOTE on training data only ===
    [X_smoted, Y_smoted_numeric] = smote(X_train_std, [], 5, 'Class', Y_numeric);  % balance all classes
    Y_smoted = Y_cat(Y_smoted_numeric);  % convert back to categorical

    % Train model
    RF = TreeBagger(500, X_smoted, Y_smoted, 'Method', 'classification', ...
        'OOBPrediction', 'On', ...
    'OOBPredictorImportance', 'on',...
    'MinLeafSize', 5, ...
    'NumPredictorsToSample', floor(sqrt(size(X_smoted,2))));

    % Predict
    [preds, scores] = predict(RF, X_test_std);
    preds = categorical(preds);

    % Confusion matrix
    C = confusionmat(Y_test, preds, 'Order', categories(Y));
    if size(C,1) < 2  % handle rare folds with 1 class
        continue;
    end
    TP = C(1,1); FN = C(1,2); FP = C(2,1); TN = C(2,2);
    acc_vals(i) = (TP + TN) / sum(C(:));
    sens_vals(i) = TP / (TP + FN);
    spec_vals(i) = TN / (TN + FP);


    % Compute AUC
    true_binary = double(Y_test == 'Diabetic');
    [~, ~, ~, auc] = perfcurve(true_binary, scores(:,1), 1);
    auc_vals(i) = auc;
end

% === Final Performance Summary ===
fprintf('AUC: %.3f (95%% CI: %.3f–%.3f)\n', mean(auc_vals), quantile(auc_vals, [0.025 0.975]));
fprintf('AUC: %.3f (std: %.3f)\n', mean(auc_vals), std(auc_vals));

fprintf('Accuracy: %.3f%% (95%% CI: %.3f–%.3f)\n', mean(acc_vals)*100, quantile(acc_vals, [0.025 0.975])*100);
fprintf('Sensitivity: %.3f%%(95%% CI: %.3f–%.3f)\n', mean(sens_vals)*100,quantile(sens_vals, [0.025 0.975]));
fprintf('Specificity: %.3f%%(95%% CI: %.3f–%.3f)\n', mean(spec_vals)*100,quantile(spec_vals, [0.025 0.975]));


%  ========= External Validation =============================================
External = readtable('/Volumes/Extreme Pro/SZNormal/Pre_diab/Individualized_beida&henan/CombinedTable_Copy.xlsx');

mu_train = mean(X);   % or mean(X_os) if training on all oversampled data
sigma_train = std(X); % or std(X_os)
sigma_train(sigma_train == 0) = 1;

X_ext = table2array(External(:, roi_cols));
Y_ext = categorical(External.FBG_Category);
Y_ext_bin = Y_ext == "Diabetic";
X_ext_std = (X_ext - mu_train) ./ sigma_train;

% Predict
[pred_ext, score_ext] = predict(RF, X_ext);
pred_ext = categorical(pred_ext);

% Confusion matrix
C_ext = confusionmat(Y_ext, pred_ext);
if size(C_ext,1) < 2, C_ext(2,2) = 0; end
TP = C_ext(1,1); FN = C_ext(1,2); FP = C_ext(2,1); TN = C_ext(2,2);
ext_acc = (TP + TN) / sum(C_ext(:));
ext_sens = TP / (TP + FN);
ext_spec = TN / (TN + FP);

precision = TP / (TP + FP + eps);
recall = TP / (TP + FN + eps);  % same as sensitivity
f1 = 2 * (precision * recall) / (precision + recall + eps);
numerator = (TP * TN) - (FP * FN);
denominator = sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN) + eps);
mcc = numerator / denominator;

[~, ~, ~, ext_auc] = perfcurve(Y_ext_bin, score_ext(:,1), 1);

ext_AUC = [ext_AUC;ext_auc];
ext_ACC = [ext_ACC;ext_acc];
ext_SEN = [ext_SEN;ext_sens];
ext_SPE = [ext_SPE;ext_spec];
% Report
fprintf('\n--- External Validation ---\n');
fprintf('AUC: %.3f\n', ext_auc);
fprintf('Accuracy: %.3f%%\n', ext_acc*100);
fprintf('Sensitivity: %.3f%%\n', ext_sens*100);
fprintf('Specificity: %.3f%%\n', ext_spec*100);
fprintf('MCC: %.3f\n', mcc);
fprintf('F1 Score: %.3f\n', f1);

end

fprintf('\n--- External Summary ---\n');
fprintf('Mean External AUC: %.3f (95%% CI: %.3f–%.3f)\n', mean(ext_AUC), quantile(ext_AUC,[0.025 0.975]));
fprintf('Mean External ACC: %.3f (95%% CI: %.3f–%.3f)\n', mean(ext_ACC), quantile(ext_ACC,[0.025 0.975]));

fprintf('Mean External SEN: %.3f (95%% CI: %.3f–%.3f)\n', mean(ext_SEN), quantile(ext_SEN,[0.025 0.975]));
fprintf('Mean External SPE: %.3f (95%% CI: %.3f–%.3f)\n', mean(ext_SPE), quantile(ext_SPE,[0.025 0.975]));


%% 11b Step: show the difference in each group 
% Step 1: In the pre-diabetic group
Pre_sum = mean(table2array(merged_pre(:,roi_cols)),1);

[~, idx] = maxk(abs(Pre_sum), 10); % Get the top 10 ROI-pairs by absolute sum
top_values = Pre_sum(idx);
top_labels = roi_cols(idx);

% Capitalize after underscore for readability
formatted_names = erase(top_labels, "_norm");  % optional: remove '_norm'
formatted_names = strrep(formatted_names, '_', '--');  % replace underscore
formatted_names = regexprep(formatted_names, '(^|-)(\w)', '${upper($2)}');  % capitalize
top_labels_clean = regexprep(formatted_names, '(?<=^|-)(\w)', '${upper($1)}');

% Set bar colors: red for positive, blue for negative
bar_colors = repmat([0 0 1], 10, 1);  % default blue
bar_colors(top_values > 0, :) = repmat([1 0 0], sum(top_values > 0), 1);  % red for positive

% Plot
fig = figure;
b = bar(top_values, 'FaceColor', 'flat');
b.CData = bar_colors;

xticks(1:10);
xticklabels(top_labels_clean);
xtickangle(45);
ylabel('Mean Δr');
title('Top 10 ROI-Pair Δr Means in Pre-diabetic Group');
grid on;

% Define output directory and filename
out_dir = '/Volumes/Extreme Pro/SZNormal/Pre_diab/Individualized';  
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

% Save figure
saveas(fig, fullfile(out_dir, 'Top10_PreDiabetic_DeltaR_BarPlot.png'));
saveas(fig, fullfile(out_dir, 'Top10_PreDiabetic_DeltaR_BarPlot.fig'));


%% Step 2: In the diabetic group
DM_sum = mean(table2array(merged_DM(:,roi_cols)),1);

[~, idx] = maxk(abs(DM_sum), 10); % Get the top 10 ROI-pairs by absolute sum
top_values = DM_sum(idx);
top_labels = roi_cols(idx);

% Capitalize after underscore for readability
formatted_names = erase(top_labels, "_norm");  % optional: remove '_norm'
formatted_names = strrep(formatted_names, '_', '--');  % replace underscore
formatted_names = regexprep(formatted_names, '(^|-)(\w)', '${upper($2)}');  % capitalize
top_labels_clean = regexprep(formatted_names, '(?<=^|-)(\w)', '${upper($1)}');

% Set bar colors: red for positive, blue for negative
bar_colors = repmat([0 0 1], 10, 1);  % default blue
bar_colors(top_values > 0, :) = repmat([1 0 0], sum(top_values > 0), 1);  % red for positive

% Plot
fig = figure;
b = bar(top_values, 'FaceColor', 'flat');
b.CData = bar_colors;

xticks(1:10);
xticklabels(top_labels_clean);
xtickangle(45);
ylabel('Mean Δr');
title('Top 10 ROI-Pair Δr Means in Diabetic Group');
grid on;

% Define output directory and filename
out_dir = '/Volumes/Extreme Pro/SZNormal/Pre_diab/Individualized';  
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

% Save figure
saveas(fig, fullfile(out_dir, 'Top10_Diabetic_DeltaR_BarPlot.png'));
saveas(fig, fullfile(out_dir, 'Top10_Diabetic_DeltaR_BarPlot.fig'));


%% 12a step: RF important features (解释性）
% Prepare feature matrix and labels

% === Step 5: Feature Importance (Top 10 only) ===
importance = RF.OOBPermutedPredictorDeltaError;
[sorted_importance, sort_idx] = sort(importance, 'descend');

top_k = 10;
top_idx = sort_idx(1:top_k);
top_features = roi_cols(top_idx);
top_scores = sorted_importance(1:top_k);

% Format feature names (replace '_' with '–' and capitalize first letter)
formatted_names = erase(top_features, "_norm");  % optional: remove '_norm'
formatted_names = strrep(formatted_names, '_', '--');  % replace underscore
formatted_names = regexprep(formatted_names, '(^|-)(\w)', '${upper($2)}');  % capitalize
formatted_names = regexprep(formatted_names, '(?<=^|-)(\w)', '${upper($1)}');

% Colors
bar_color = [112, 128, 144]/255;  % subtle slate gray (neutral)
% bar_color = [178 24 43] / 255;  % Nature red
% bar_color = [33 102 172] / 255;  % Nature blue

text_font = 'Helvetica';

fig = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.2 0.3 0.5 0.5]);

% Plot
bar(top_scores, 'FaceColor', bar_color, 'EdgeColor', 'k', 'LineWidth', 0.8);

% Axes
xticks(1:top_k);
xticklabels(formatted_names);
xtickangle(45);
ylabel('Mean Absolute SHAP Contribution (A.U.)', ...
       'FontName', text_font, 'FontSize', 10);

% Axis styling
ax = gca;
ax.FontName = text_font;
ax.FontSize = 9;
ax.LineWidth = 1.2;
ax.XColor = 'k';
ax.YColor = 'k';
ax.TickDir = 'out';
ax.Box = 'off';
ax.YGrid = 'on';
ax.GridLineStyle = '--';


% Save figure
saveas(fig, '/Volumes/Extreme Pro/SZNormal/Pre_diab/Individualized/Top10_Features_PrevsDM.png');   % Save as PNG
savefig(fig, '/Volumes/Extreme Pro/SZNormal/Pre_diab/Individualized/Top10_Features_PrevsDM.fig');  % Save as MATLAB .fig

%% 12b: correlation of feature to FBG
[r, p] = corr(X_mat, CombinedTable.FBG, 'type', 'Spearman');
% Identify significant indices
sig_idx = find(p < 0.05);

% Extract values
sig_rois = roi_cols(sig_idx);
sig_r = r(sig_idx);
sig_p = p(sig_idx);

% Create table
SignificantCorrelationTable = table(sig_rois, sig_r, sig_p, ...
    'VariableNames', {'Organ_Pair', 'Spearman_r', 'p_value'});

% Display the table
disp(SignificantCorrelationTable);

%% Function in Use
function q = bh_fdr(p)
    % Input: p - vector of raw p-values
    % Output: q - vector of FDR-corrected p-values (q-values)
    
    p = p(:);
    [sorted_p, sort_idx] = sort(p);
    m = length(p);
    q = zeros(size(p));
    
    for i = 1:m
        q(sort_idx(i)) = sorted_p(i) * m / i;
    end
    
    % Ensure q-values are monotonic
    for i = m-1:-1:1
        q(sort_idx(i)) = min(q(sort_idx(i)), q(sort_idx(i+1)));
    end
    
    % Cap at 1
    q = min(q, 1);
end


% Clustering function for binary undirected graphs
function C = clustering_coef_bu(A)
    K = sum(A);  % degree
    C = zeros(size(A,1),1);
    for u = 1:length(K)
        neighbors = find(A(u,:));
        k = length(neighbors);
        if k >= 2
            subgraph = A(neighbors, neighbors);
            C(u) = sum(subgraph(:)) / (k*(k-1));
        end
    end
end


% Capitalize first letter 
function out = capitalize_first(str)
    if isempty(str)
        out = str;
    else
        out = [upper(str(1)), str(2:end)];
    end
end


% 在图片里将p<0.05的用* 表示
function s = get_star(p)
    if p < 0.001
        s = '***';
    elseif p < 0.01
        s = '**';
    elseif p < 0.05
        s = '*';
    else
        s = 'n.s.';
    end
end

% Define bootstrap function (packaged) ---------------
function [R_mean,z_mean] = bootstrap_network(data, n_boot, alpha)
    n_rois = size(data, 2);
    z_boot = @(r) 0.5 * log((1 + r) ./ (1 - r));       % Fisher z-transform
    inv_z  = @(z) (exp(2*z) - 1) ./ (exp(2*z) + 1);     % Inverse Fisher

    z_stack = zeros(n_rois, n_rois, n_boot);           % Fisher z-matrices
   

    for b = 1:n_boot
        idx = randi(size(data, 1), size(data, 1), 1);   % Bootstrap resampling
        sample = data(idx, :);

        [R, P] = corr(sample, 'Rows', 'pairwise', 'Type', 'Pearson');
        R(isnan(R)) = 0; P(isnan(P)) = 1;

        % FDR correction
        upper_idx = find(triu(true(n_rois), 1));
        raw_p = P(upper_idx);
        [~, ~, ~, fdr_adj] = fdr_bh(raw_p, alpha, 'pdep', 'yes');

        % Binary mask for significant + strong correlations
        sig = false(n_rois);
        sig(upper_idx) = (fdr_adj < alpha);
        sig = sig | sig';  % Symmetrize

        % Keep significant values only
        R_sig = R .* sig;

        % Apply Fisher transform to filtered matrix
        R_sig(R_sig == 0) = NaN;  % remove zeros (non-significant) from averaging
       
        % Store z-transformed matrix
        z_stack(:, :, b) = z_boot(R_sig);
    end

    % Average correlation across bootstraps
    z_mean = mean(z_stack, 3, 'omitnan');
    R_mean = inv_z(z_mean); % inverse Fisher transform


end


% --- Define a function to compute network metrics in one single network---
function [metrics, node_strength, node_clustering] = compute_network_metrics(R, A)
    n = size(R,1);  % number of nodes

    % Network Density
    density = sum(A(:)) / (n*(n-1));  % undirected, no self-loops

    % Degree: number of connections per node
    degree = sum(A, 2);  % node degree vector
    avg_degree = mean(degree);

    % Strength: sum of absolute correlation weights of connections
    R(isnan(R)) = 0;         % remove NaNs
    R(1:size(R,1)+1:end) = 1;  % set diagonal to 1
    strength = sum(abs(R .* A), 2);  % weighted node strength
    avg_strength = mean(strength);

    % Clustering Coefficient
    clustering = zeros(n,1);
    for i = 1:n
        neighbors = find(A(i,:));
        k = length(neighbors);
        if k < 2
            clustering(i) = 0;
        else
            subgraph = A(neighbors, neighbors);
            actual_edges = sum(subgraph(:)) / 2;
            possible_edges = k * (k - 1) / 2;
            clustering(i) = actual_edges / possible_edges;
        end
    end
    avg_clustering = mean(clustering);

    % Store summary metrics
    metrics = struct( ...
        'Density', density, ...
        'AvgDegree', avg_degree, ...
        'AvgStrength', avg_strength, ...
        'AvgClustering', avg_clustering);

    % Return node-level metrics
    node_strength = strength;
    node_clustering = clustering;
end


% ---------function to compute node metrics across boots------
function [strength_mat, cluster_mat] = compute_node_metrics_boot(R_stack, A_stack)
    n_boot = size(R_stack, 3);
    n_rois = size(R_stack, 1);

    strength_mat = zeros(n_boot, n_rois);
    cluster_mat  = zeros(n_boot, n_rois);

    for b = 1:n_boot
        R = R_stack(:,:,b);
        A = A_stack(:,:,b);

        R(isnan(R)) = 0;
        R(1:n_rois+1:end) = 1;  % set diagonal to 1

        % Strength
        strength_mat(b,:) = sum(abs(R .* A), 2)';

        % Clustering
        cluster_vals = zeros(n_rois,1);
        for i = 1:n_rois
            neighbors = find(A(i,:));
            k = length(neighbors);
            if k >= 2
                subgraph = A(neighbors, neighbors);
                actual_edges = sum(subgraph(:))/2;
                possible_edges = k*(k-1)/2;
                cluster_vals(i) = actual_edges / possible_edges;
            end
        end
        cluster_mat(b,:) = cluster_vals';
    end
end

% Mantel test
function [r_mantel, p_value] = mantel_test(A, B, n_perm)
    % A, B: input square symmetric matrices (e.g., correlation or distance)
    % n_perm: number of permutations
    % Output: r_mantel - observed correlation, p_value - permutation p

    % Make sure matrices are the same size
    if ~isequal(size(A), size(B))
        error('Matrices must be the same size');
    end

    % Get upper triangle (excluding diagonal)
    mask = triu(true(size(A)), 1);
    a_vals = A(mask);
    b_vals = B(mask);

    % Remove NaNs
    valid = ~isnan(a_vals) & ~isnan(b_vals);
    a_vals = a_vals(valid);
    b_vals = b_vals(valid);

    % Observed Pearson correlation
    r_mantel = corr(a_vals, b_vals, 'Type', 'Pearson');

    % Permutation test
    n = size(A, 1);
    r_perm = zeros(n_perm, 1);
    for i = 1:n_perm
        perm_idx = randperm(n);
        B_perm = B(perm_idx, perm_idx);
        b_perm_vals = B_perm(mask);
        b_perm_vals = b_perm_vals(valid);
        r_perm(i) = corr(a_vals, b_perm_vals, 'Type', 'Pearson');
    end

    % Calculate empirical p-value
    p_value = mean(abs(r_perm) >= abs(r_mantel));
end

function cmap = bluewhitered(m)
% Returns a diverging colormap: blue–white–red, centered at zero
if nargin < 1
    m = 256;
end
bottom = [0 0 1];
middle = [1 1 1];
top = [1 0 0];
n = fix(m/2);
cmap = [interp1([0 1], [bottom; middle], linspace(0,1,n)); ...
        interp1([0 1], [middle; top], linspace(0,1,m-n))];
end
