%%% and create lookup table
%%% by Age, Gender, BMI category, and Fasting Blood Glucose (FBG) level

addpath(genpath('/Volumes/Extreme Pro/SZNormal/Code'));

%% Step 1a: Load PatientInfo data （武警医院）
PatientInfo = readtable("/Volumes/Extreme Pro/SZNormal/Info/武警医院表.xlsx");
PatientInfo.Properties.VariableNames = {'Patient', 'ExaminationSeries', 'Name', 'Gender', 'Age', 'Weight', 'Height', 'FBG', 'Dose','BMI'};
PatientInfo = rmmissing(PatientInfo, 'DataVariables', 'BMI'); % Remove Nan
PatientInfo.Gender = categorical(PatientInfo.Gender);
PatientInfo = PatientInfo(~isnan(PatientInfo.FBG), :); % Clean EmptyGlucose (remove NaN)

if ~iscategorical(PatientInfo.Gender) % Convert Gender to categorical
    PatientInfo.Gender = categorical(PatientInfo.Gender);
end


PatientInfo.BMI_Category = strings(height(PatientInfo), 1); % BMI Category
PatientInfo.BMI_Category(PatientInfo.BMI < 18.5) = "Underweight";
PatientInfo.BMI_Category(PatientInfo.BMI >= 18.5 & PatientInfo.BMI < 23.0) = "Normal";
PatientInfo.BMI_Category(PatientInfo.BMI >= 23.0 & PatientInfo.BMI < 25.0) = "Overweight";
PatientInfo.BMI_Category(PatientInfo.BMI >= 25.0 ) = "Obese";
% PatientInfo.BMI_Category(PatientInfo.BMI >= 30.0) = "Obese_II";
PatientInfo.BMI_Category = categorical(PatientInfo.BMI_Category);

% Age Category
% edges = [20, 30, 40, 50, 60, 70, 100]; % Age ranges from <20, 20-29, ..., 80-89
% labels = {'20-29', '30-39', '40-49', '50-59', '60-69', '>70'}; % Age group labels

% edges = [10, 30, 40, 50, 60,70, 100]; % Age ranges from <20, 20-29, ..., 80-89
% labels = {'<30', '30-39', '40-49','50-59','60-69', '≥70'}; % Age group labels

edges = [10, 40, 60, 80, 100]; % Age ranges from <20, 20-29, ..., 80-89
labels = {'20-39', '40-59', '60-79','≥80'}; % Age group labels
% 20-39 Peak organ function， 40-59 Early metabolic changes， 60–79
% Aging-related decline，80+ Frailty, comorbidities 
% Ideal when aligning with clinical cutoffs, e.g., for dementia, cancer screening, etc.
PatientInfo.Age_Group = discretize(PatientInfo.Age, edges, 'categorical', labels);

% Fasting Blood Glucose (FBG) Category
% Normal: <5.6, Pre-diabetes: 5.6–6.9, Diabetes: ≥7.0 (mmol/L)
glucose = PatientInfo.FBG;
PatientInfo.FBG_Category = strings(height(PatientInfo),1);
PatientInfo.FBG_Category(glucose < 5.6) = "Normal";
PatientInfo.FBG_Category(glucose >= 5.6 & glucose < 7.0) = "Pre-diabetic";
PatientInfo.FBG_Category(glucose >= 7.0) = "Diabetic";
PatientInfo.FBG_Category = categorical(PatientInfo.FBG_Category);

% Join PatientInfo with SUVTable on 'Patient'
SUVTable = readtable('/Volumes/Extreme Pro/SZNormal/Extracted_SUV_meanmedianmaxmin/SUL_statistics_1000.xlsx', ...
    'Sheet','Median');

roi_cols = SUVTable.Properties.VariableNames(2:end);  % skip SubjectID
reference = SUVTable.liver;

for i = 1:length(roi_cols)
    col = roi_cols{i};
    if ~strcmp(col, 'liver')
        SUVTable.([col '_norm']) = SUVTable.(col) ./ reference;
    end
end

FullData_healthy = innerjoin(PatientInfo, SUVTable, 'Keys', 'Patient'); % Inner join by 'Patient' column

output_filename = '/Volumes/Extreme Pro/SZNormal/LookupTable/Median_SUL_1000.xlsx';
writetable(FullData_healthy, output_filename);


metadata_cols = {'Patient','ExaminationSeries', 'Name','Age', 'Gender', 'Weight', 'Height', 'Dose','FBG', ...
                 'BMI', 'BMI_Category', 'Age_Group', 'FBG_Category'};


%% Step 1b: Load PatientInfo data (北大深圳） 
PatientInfo = readtable("/Volumes/Extreme Pro/SZNormal/Info/54综合表.xlsx");
PatientInfo.Properties.VariableNames = {'Patient', 'ExaminationSeries', 'Name', 'Gender', 'Age', 'Weight', 'Height', 'FBG', 'Dose','BMI'};
% Remove Nan
PatientInfo = rmmissing(PatientInfo, 'DataVariables', 'BMI');
PatientInfo.Gender = categorical(PatientInfo.Gender);
% Clean EmptyGlucose (optional: remove NaN)
PatientInfo = PatientInfo(~isnan(PatientInfo.FBG), :);

% Convert Gender to categorical if needed
if ~iscategorical(PatientInfo.Gender)
    PatientInfo.Gender = categorical(PatientInfo.Gender);
end

% Step 2: Convert to Categories
% BMI Category
PatientInfo.BMI_Category = strings(height(PatientInfo), 1);

PatientInfo.BMI_Category(PatientInfo.BMI < 18.5) = "Underweight";
PatientInfo.BMI_Category(PatientInfo.BMI >= 18.5 & PatientInfo.BMI < 23.0) = "Normal";
PatientInfo.BMI_Category(PatientInfo.BMI >= 23.0 & PatientInfo.BMI < 25.0) = "Overweight";
PatientInfo.BMI_Category(PatientInfo.BMI >= 25.0 ) = "Obese";
% PatientInfo.BMI_Category(PatientInfo.BMI >= 30.0) = "Obese_II";
PatientInfo.BMI_Category = categorical(PatientInfo.BMI_Category);


edges = [10, 40, 60, 80, 100]; % Age ranges from <20, 20-29, ..., 80-89
labels = {'20-39', '40-59', '60-79','≥80'}; % Age group labels

PatientInfo.Age_Group = discretize(PatientInfo.Age, edges, 'categorical', labels);

% Fasting Blood Glucose (FBG) Category
% Normal: <5.6, Pre-diabetes: 5.6–6.9, Diabetes: ≥7.0 (mmol/L)
glucose = PatientInfo.FBG;
PatientInfo.FBG_Category = strings(height(PatientInfo),1);
PatientInfo.FBG_Category(glucose < 5.6) = "Normal";
PatientInfo.FBG_Category(glucose >= 5.6 & glucose < 7.0) = "Pre-diabetic";
PatientInfo.FBG_Category(glucose >= 7.0) = "Diabetic";
PatientInfo.FBG_Category = categorical(PatientInfo.FBG_Category);

% Join PatientInfo with SUVTable on 'Patient'
SUVTable = readtable('/Volumes/Extreme Pro/SZNormal/Extracted_SUV_meanmedianmaxmin/SUL_statistics_54.xlsx', ...
    'Sheet','Mean');
% Inner join by 'Patient' column
FullData_healthy_54 = innerjoin(PatientInfo, SUVTable, 'Keys', 'Patient');

output_filename = '/Volumes/Extreme Pro/SZNormal/LookupTable/Mean_SUL_54.xlsx';
writetable(FullData_healthy_54, output_filename);

%% 河南省医院
PatientInfo = readtable("/Volumes/Extreme Pro/SZNormal/DICOM_Metadata.xlsx");
PatientInfo.Properties.VariableNames = {'Patient', 'Series', 'Gender',  'Age', 'Height', 'Weight', 'BMI', 'Dose (mCi)', 'Dose', 'FBG'};
% Remove Nan
PatientInfo = rmmissing(PatientInfo, 'DataVariables', 'BMI');
PatientInfo = rmmissing(PatientInfo, 'DataVariables', 'Patient');
PatientInfo.Gender = categorical(PatientInfo.Gender);
% Clean EmptyGlucose (optional: remove NaN)
PatientInfo = PatientInfo(~isnan(PatientInfo.FBG), :);

% Convert Gender to categorical if needed
if ~iscategorical(PatientInfo.Gender)
    PatientInfo.Gender = categorical(PatientInfo.Gender);
end

% Step 2: Convert to Categories
% BMI Category
PatientInfo.BMI_Category = strings(height(PatientInfo), 1);

PatientInfo.BMI_Category(PatientInfo.BMI < 18.5) = "Underweight";
PatientInfo.BMI_Category(PatientInfo.BMI >= 18.5 & PatientInfo.BMI < 23.0) = "Normal";
% PatientInfo.BMI_Category(PatientInfo.BMI >= 23.0 & PatientInfo.BMI < 25.0) = "Overweight";
% PatientInfo.BMI_Category(PatientInfo.BMI >= 25.0 ) = "Obese";
% PatientInfo.BMI_Category(PatientInfo.BMI >= 30.0) = "Obese_II";
PatientInfo.BMI_Category(PatientInfo.BMI >= 23.0 ) = "OW&Obese";

PatientInfo.BMI_Category = categorical(PatientInfo.BMI_Category);


edges = [10, 40, 60, 80, 100]; % Age ranges from <20, 20-29, ..., 80-89
labels = {'20-39', '40-59', '60-79','≥80'}; % Age group labels

PatientInfo.Age_Group = discretize(PatientInfo.Age, edges, 'categorical', labels);

% Fasting Blood Glucose (FBG) Category
% Normal: <5.6, Pre-diabetes: 5.6–6.9, Diabetes: ≥7.0 (mmol/L)
glucose = PatientInfo.FBG;
PatientInfo.FBG_Category = strings(height(PatientInfo),1);
PatientInfo.FBG_Category(glucose < 5.6) = "Normal";
PatientInfo.FBG_Category(glucose >= 5.6 & glucose < 7.0) = "Pre-diabetic";
PatientInfo.FBG_Category(glucose >= 7.0) = "Diabetic";
PatientInfo.FBG_Category = categorical(PatientInfo.FBG_Category);

% Join PatientInfo with SUVTable on 'Patient'
SUVTable = readtable('/Volumes/Extreme Pro/SZNormal/Extracted_SUV_meanmedianmaxmin/SUL_statistics_pathology.xlsx', ...
    'Sheet','Mean');
% Inner join by 'Patient' column
FullData_healthy = innerjoin(PatientInfo, SUVTable, 'Keys', 'Patient');

output_filename = '/Volumes/Extreme Pro/SZNormal/LookupTable/Mean_SUL_pathology_OW&OB.xlsx';
writetable(FullData_healthy, output_filename);

%% Step 2: Statistical Group comparison per Subgroup
% Do statistical analysis showing whether SUV change between subgroups
% Grouping variables
group_vars = {'Age_Group', 'Gender', 'BMI_Category', 'FBG_Category'};
anova_results = struct();

for g = 1:numel(group_vars)
    group_name = group_vars{g};
    fprintf('\n\n--- Analyzing SUV Differences by %s ---\n', group_name);

    for r = 1:numel(roi_cols)
        roi = roi_cols{r};
        y = FullData_healthy.(roi);
        if iscell(y)
            y = cellfun(@double, y);  % Convert cells to numbers
        end
        g_labels = FullData_healthy.(group_name);
        
        valid_idx = ~isnan(y) & ~ismissing(g_labels);
        y = y(valid_idx);
        g_labels = g_labels(valid_idx);


        % Check normality per group
        unique_groups = categories(categorical(g_labels));
        normal_flags = true(numel(unique_groups),1);

        for k = 1:numel(unique_groups)
            group_data = y(g_labels == unique_groups{k});
            if numel(group_data) >= 4  % Lilliefors requires at least 4 points
                normal_flags(k) = lillietest(group_data) == 0;  % 0 = normal
            else
                normal_flags(k) = false;  % not enough data
            end
        end

        all_groups_normal = all(normal_flags);

        if all_groups_normal && numel(unique_groups) > 1
            p = anova1(y, g_labels, 'off');
            test_used = 'ANOVA';
        else
            p = kruskalwallis(y, g_labels, 'off');
            test_used = 'Kruskal-Wallis';
        end

        % Save results
        anova_results.(group_name).(roi).p_value = p;
        anova_results.(group_name).(roi).test = test_used;

        fprintf('%s (ROI: %s): p = %.4f [%s]\n', group_name, roi, p, test_used);
    end
end

% Step 4: Adjust for Multiple Comparisons (FDR) 
all_p = [];
roi_labels = {};
group_labels = {};

for g = 1:numel(group_vars)
    for r = 1:numel(roi_cols)
        p = anova_results.(group_vars{g}).(roi_cols{r}).p_value;
        all_p(end+1) = p;
        roi_labels{end+1} = roi_cols{r};
        group_labels{end+1} = group_vars{g};
    end
end

% Apply FDR correction
[~, ~, adj_p] = fdr_bh(all_p);

% Report significant findings
fprintf('\n--- Significant Differences after FDR Correction (q < 0.05) ---\n');
for i = 1:length(all_p)
    pvals = all_p(i);
    [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals);
    if adj_p < 0.05
        fprintf('ROI: %-15s | Group: %-15s | adj. p = %.4f\n', ...
            roi_labels{i}, group_labels{i}, adj_p);
    end
end

%% Visualization
roi_to_plot = 'AdrenalGlandLeft';
group_by = 'Age_Group';  % or Age_Group, Gender, FBG_Category

figure;
boxplot(FullData.(roi_to_plot), FullData.(group_by));
title(sprintf('%s SUV by %s', roi_to_plot, group_by), 'Interpreter', 'none');
ylabel('SUV');
xlabel(group_by);

%% Step 5a:  LookUp Table Descriptive
% group_vars = {'Age_Group', 'Gender', 'BMI_Category', 'FBG_Category'};
group_vars = {'Age_Group', 'Gender', 'BMI_Category'};
[G, group_labels] = findgroups(FullData(:, group_vars));

% Prepare cell array to hold rows
summary_data = {};

for g = 1:max(G)
    idx = G == g;
    group_data = FullData(idx, :);
    label_values = table2cell(group_labels(g,:));  % e.g., {'30–39', 'Male', 'Normal', 'Diabetic'}

    % Store sample size
    n_participants = sum(idx);

    % Check normality of each ROI
    normality_flags = cell(1, numel(roi_cols));
    for r = 1:numel(roi_cols)
        x = group_data.(roi_cols{r});
        if iscell(x)
            x = cellfun(@double, x);
        end
        x = x(~isnan(x));
        
        if numel(x) < 4
            normality_flags{r} = "Too few";
        else
            try
                [h, ~] = swtest(x);  % Shapiro-Wilk: h = 0 means normal
                if h == 0
                    normality_flags{r} = "Normal";
                else
                    normality_flags{r} = "Not normal";
                end
            catch
                 % Fallback to Lilliefors if swtest is unavailable or fails
                h = lillietest(x);
                if h == 0
                    normality_flags{r} = "Normal*";  % fallback method
                else
                    normality_flags{r} = "Not normal*";
                end
            end
        end
    end

    % Combine this row
    summary_data(end+1, :) = [label_values, {n_participants}, normality_flags];
end

% Construct table
summary_table = cell2table(summary_data, ...
    'VariableNames', [group_vars, {'N'}, roi_cols]);


output_filename = '/Volumes/Extreme Pro/SZNormal/Extracted_SUV_meanmedianmaxmin/Data_exploration.xlsx';
writetable(summary_table, output_filename);

%% Step 5b: Build LookUp Table
% group_vars = {'Age_Group', 'Gender', 'BMI_Category', 'FBG_Category'};
group_vars = {'Age_Group', 'Gender', 'BMI_Category'};
[G, group_labels] = findgroups(FullData(:, group_vars));  % Create group indices

lookup_rows = [];  % Will collect each row of the lookup table
lookup_headers = [group_vars, roi_cols];  % Final column names
lookup_data = {};  % Cell array for row-wise results

for g = 1:max(G)
    idx = (G == g);
    group_data = FullData(idx, :);

    % Group label values (e.g., Age_Group = "30-39", etc.)
    label_values = table2cell(group_labels(g,:));

    % Compute summary stat for each ROI
    roi_summaries = cell(1, numel(roi_cols));
    for r = 1:numel(roi_cols)
        roi_vals = group_data.(roi_cols{r});
        if iscell(roi_vals)
            roi_vals = cellfun(@double, roi_vals);  % Convert from cell if needed
        end
        roi_summaries{r} = summaryStat(roi_vals);
    end

    % Append to result table
    lookup_data(end+1, :) = [label_values, roi_summaries];
end

% Convert to table
lookup_table = cell2table(lookup_data, 'VariableNames', lookup_headers);

output_filename = '/Volumes/Extreme Pro/SZNormal/Result/LookupTable_mean_withoutFBG.xlsx';
writetable(lookup_table, output_filename);

%% Step 6a: Build LookUpTable only for normal individual
normal_idx = (FullData.BMI_Category == "Normal") & (FullData.FBG_Category == "Normal");
FullData_normal = FullData(normal_idx, :);

group_vars = {'Age_Group', 'Gender'};
[G, group_labels] = findgroups(FullData_normal(:, group_vars));  % Create group indices

lookup_rows = [];  % Will collect each row of the lookup table
lookup_headers = [group_vars, roi_cols];  % Final column names
lookup_data_normal = {};  % Cell array for row-wise results

for g = 1:max(G)
    idx = (G == g);
    group_data = FullData_normal(idx, :);

    % Group label values (e.g., Age_Group = "30-39", etc.)
    label_values = table2cell(group_labels(g,:));

    % Compute summary stat for each ROI
    roi_summaries = cell(1, numel(roi_cols));
    for r = 1:numel(roi_cols)
        roi_vals = group_data.(roi_cols{r});
        if iscell(roi_vals)
            roi_vals = cellfun(@double, roi_vals);  % Convert from cell if needed
        end
        roi_summaries{r} = summaryStat(roi_vals);
    end

    % Append to result table
    lookup_data_normal(end+1, :) = [label_values, roi_summaries];
end

% Convert to table
lookup_table_normal = cell2table(lookup_data_normal, 'VariableNames', lookup_headers);

output_filename = '/Volumes/Extreme Pro/SZNormal/Result/LookupTable_normal.xlsx';
writetable(lookup_table_normal, output_filename,'Sheet', 'Mean');


%% Define function
function s = summaryStat(x)
    % Computes either mean ± 3 std (if normal), or median (IQR) (if not)
    % Uses Shapiro-Wilk test (swtest.m must be on path)

    x = x(~isnan(x));  % Remove NaNs

    if numel(x) < 4
        s = "Insufficient data";
        return;
    end

    mu = mean(x);
    sigma = std(x);
    s = sprintf('%.2f ± %.2f', mu, 3*sigma);
   
end
