%%% Extract SUV mean, max, min, median

%% Extract SUV min, max, mean, median 
% Initialize arrays to store statistics for each ROI
patients_path = '/Volumes/Extreme Pro/SZNormal/SUV_Distribution/SUL_distribution_pathology';
imagedir = dir(patients_path);
patnames_ = {imagedir(1:end).name}';
patnames = patnames_(~startsWith(patnames_, '.')); % All patients 
numeric_parts = cellfun(@(x) str2double(x(2:end-4)), patnames);
[sorted_values, sort_order] = sort(numeric_parts);
patnames = patnames(sort_order);

load('/Volumes/Extreme Pro/SZNormal/SUV_Distribution/SUL_distribution_pathology/001_WANG_YOU_FANG.mat');
ROI_names = featcell(1,:);

num_patients = length(patnames);
% Initialize four tables
mean_SUV = zeros(num_patients, length(ROI_names));
median_SUV = zeros(num_patients, length(ROI_names));
max_SUV = zeros(num_patients, length(ROI_names));


% Loop through each patient
for i = 1:num_patients
    files_path = fullfile(patients_path,patnames{i});
    load(files_path);
    disp(['Loaded data for: ', patnames{i}]);

    % Loop through each ROI
    for j = 1:length(featcell)
        SUV_values = featcell{2,j};
        if isempty(SUV_values)
            SUV_values = featcell{2,j+1};
        end
        mean_SUV(i,j) = mean(SUV_values);
        median_SUV(i,j) = median(SUV_values);
        max_SUV(i,j) = max(SUV_values);
        
    end
end

%% write in table
output_filename = '/Volumes/Extreme Pro/SZNormal/Extracted_SUV_meanmedianmaxmin/SUL_statistics_pathology.xlsx';

patnames_clean = erase(patnames, '.mat'); 

% Convert numeric matrices to tables and assign column names
mean_SUV_table = array2table(mean_SUV, 'VariableNames', matlab.lang.makeValidName(ROI_names));
median_SUV_table = array2table(median_SUV, 'VariableNames', matlab.lang.makeValidName(ROI_names));
max_SUV_table = array2table(max_SUV, 'VariableNames', matlab.lang.makeValidName(ROI_names));

% Add patient names as the first column
mean_SUV_table = addvars(mean_SUV_table, patnames_clean, 'Before', 1, 'NewVariableNames', 'Patient');
median_SUV_table = addvars(median_SUV_table, patnames_clean, 'Before', 1, 'NewVariableNames', 'Patient');
max_SUV_table = addvars(max_SUV_table, patnames_clean, 'Before', 1, 'NewVariableNames', 'Patient');

writetable(mean_SUV_table, output_filename, 'Sheet', 'Mean');
writetable(median_SUV_table, output_filename, 'Sheet', 'Median');
writetable(max_SUV_table, output_filename, 'Sheet', 'Max');
disp('Tables have been saved to Excel files.');