%% c1 extract static SUV second batch
addpath(genpath('/Volumes/Extreme Pro/糖代谢/code'));
addpath(genpath('/Volumes/Extreme Pro/糖代谢/code姜然'));

%% extract static SUV distribution within given ROIs （54）
PET_dir_name = '/Volumes/Extreme Pro/SZNormal/BDSZ_54_PET_CT_Brain_Organ_Cardiac/BDSZ_54_PET_CT_Brain_Organ_Cardiac';
imagedir = dir(PET_dir_name);
patnames_ = {imagedir(1:end).name}';
patnames = patnames_(~startsWith(patnames_, '.')); % All patients 

label_table = readtable('/Volumes/Extreme Pro/SZNormal/Info/BDSZ_54例_生理信息及对应序号.xlsx');
label_table.Properties.VariableNames = {'Label', 'ExaminationSeries', 'Series','Name', 'Gender', 'Age', 'Weight', 'Height', 'EmptyGlucose', 'mCi','Dose','BMI'};

out_dir_name = '/Volumes/Extreme Pro/SZNormal/SUV_Distribution/SUV_distribution_54';
mask_dir_name = '/Volumes/Extreme Pro/SZNormal/masks/masks_54';

% Loop through each patient:
for i = 1:numel(patnames) 
    disp(['Working on patient ' patnames{i} '...']);
    mask_dir_ = dir(fullfile(mask_dir_name,patnames{i}));
    mask_dir = mask_dir_(~startsWith({mask_dir_.name}, '.'));

    pet_folder = fullfile(PET_dir_name, patnames{i});
    pet_files = dir(pet_folder);
    pet_file = pet_files(~startsWith({pet_files.name}, '.')&endsWith({pet_files.name},'.nii')&startsWith({pet_files.name},'PT'));
    pet_path = fullfile(pet_file.folder, pet_file.name);

    % retrieve injected dose and body weight
    injected_dose_Bq = label_table.Dose(strcmp(label_table.Label,patnames{i}))* 1e6;
    body_weight_g = label_table.Weight(strcmp(label_table.Label,patnames{i}))*1000;
    gender = label_table.Gender(strcmp(label_table.Label,patnames{i}));
    body_BMI = label_table.BMI(strcmp(label_table.Label,patnames{i}));
    
    % Load the raw PET image (NIfTI file)
    pet_image_data = niftiread(pet_path);

    % Calculate SUV: SUV = (PET voxel value) / (Injected dose / Body weight)
     SUV_image = double(pet_image_data)* body_weight_g /injected_dose_Bq;
     vol_SUV = double(SUV_image);

    % Calculate SUV distribution per Body ROI
    featcell = cell([2 numel(mask_dir)]); 
    featcell(1,1:end) = strrep(erase({mask_dir(1:end).name}','.nii'),'_',' ');

    for j = 1:numel(mask_dir)
        disp(['Extracting from ROI: ' mask_dir(j).name '...']);
        mask = niftiread(fullfile(mask_dir(j).folder,mask_dir(j).name));

        SUV_distr = vol_SUV(mask==1);
        if size(mask,3) == size(vol_SUV,3) % check vol and mask length
            featcell{2,j} = SUV_distr;
        else
            featcell{2,j} = [];
            disp('Vol and mask not aligned!');
        end
    end
    save(fullfile(out_dir_name,patnames{i}),"featcell");
end

%% extract static SUL distribution within given ROIs （54）
PET_dir_name = '/Volumes/Extreme Pro/SZNormal/BDSZ_54_PET_CT_Brain_Organ_Cardiac/BDSZ_54_PET_CT_Brain_Organ_Cardiac';
imagedir = dir(PET_dir_name);
patnames_ = {imagedir(1:end).name}';
patnames = patnames_(~startsWith(patnames_, '.')); % All patients 

label_table = readtable('/Volumes/Extreme Pro/SZNormal/Info/BDSZ_54例_生理信息及对应序号.xlsx');
label_table.Properties.VariableNames = {'Label', 'ExaminationSeries', 'Series','Name', 'Gender', 'Age', 'Weight', 'Height', 'EmptyGlucose', 'mCi','Dose','BMI'};

out_dir_name = '/Volumes/Extreme Pro/SZNormal/SUV_Distribution/SUL_distribution_54';
mask_dir_name = '/Volumes/Extreme Pro/SZNormal/masks/masks_54';

% Loop through each patient:
for i = 1:numel(patnames) 
    disp(['Working on patient ' patnames{i} '...']);
    mask_dir_ = dir(fullfile(mask_dir_name,patnames{i}));
    mask_dir = mask_dir_(~startsWith({mask_dir_.name}, '.'));

    pet_folder = fullfile(PET_dir_name, patnames{i});
    pet_files = dir(pet_folder);
    pet_file = pet_files(~startsWith({pet_files.name}, '.')&endsWith({pet_files.name},'.nii')&startsWith({pet_files.name},'PT'));
    pet_path = fullfile(pet_file.folder, pet_file.name);

    % retrieve injected dose and body weight
    injected_dose_Bq = label_table.Dose(strcmp(label_table.Label,patnames{i}))* 1e6;
    body_weight_kg = label_table.Weight(strcmp(label_table.Label,patnames{i}));
    gender = label_table.Gender(strcmp(label_table.Label,patnames{i}));
    body_BMI = label_table.BMI(strcmp(label_table.Label,patnames{i}));
    
    % Load the raw PET image (NIfTI file)
    pet_image_data = double(niftiread(pet_path));

    if strcmp(gender,'M') % for male
        LBM = 9270 * body_weight_kg / (6680 + 216 * body_BMI);
    else
        LBM = 9270 * body_weight_kg / (8780 + 244 * body_BMI);
    end

    LBM_g = LBM * 1000;

    SUL = pet_image_data * LBM_g / injected_dose_Bq;
    
    vol_SUL = double(SUL);

    % Calculate SUV distribution per Body ROI
    featcell = cell([2 numel(mask_dir)]); 
    featcell(1,1:end) = strrep(erase({mask_dir(1:end).name}','.nii'),'_',' ');

    for j = 1:numel(mask_dir)
        disp(['Extracting from ROI: ' mask_dir(j).name '...']);
        mask = niftiread(fullfile(mask_dir(j).folder,mask_dir(j).name));

        SUV_distr = vol_SUL(mask==1);
        if size(mask,3) == size(vol_SUL,3) % check vol and mask length
            featcell{2,j} = SUV_distr;
        else
            featcell{2,j} = [];
            disp('Vol and mask not aligned!');
        end
    end
    save(fullfile(out_dir_name,patnames{i}),"featcell");
end


%% extract static SUV distribution within given ROIs (1000)
PET_dir_name = '/Volumes/Extreme Pro/SZNormal/SZWJ_1185_PET_CT_Brain_Organ_Cardiac/SZWJ_1185_PET_CT_Brain_Organ_Cardiac';
imagedir = dir(PET_dir_name);
patnames_ = {imagedir(1:end).name}';
patnames = patnames_(~startsWith(patnames_, '.')); % All patients 

label_table = readtable('/Volumes/Extreme Pro/SZNormal/Info/szwj1189.xlsx','Range','A2:M1190');
label_table.Properties.VariableNames = {'Label', 'ExaminationSeries','Name', ...
    'Age', 'Gender','Weight', 'measured_mm', 'EmptyGlucose', 'Dose1','Dose',...
    'Height', 'BMI','categories'};

out_dir_name = '/Volumes/Extreme Pro/SZNormal/SUV_Distribution/SUV_ditribution_1000';
mask_dir_name = '/Volumes/Extreme Pro/SZNormal/masks/masks_1000';

% Loop through each patient:
for i = 1:numel(patnames) 
    disp(['Working on patient ' patnames{i} '...']);
    mask_dir_ = dir(fullfile(mask_dir_name,patnames{i}));
    mask_dir = mask_dir_(~startsWith({mask_dir_.name}, '.'));

    pet_folder = fullfile(PET_dir_name, patnames{i});
    pet_files = dir(pet_folder);
    pet_file = pet_files(~startsWith({pet_files.name}, '.')&endsWith({pet_files.name},'.nii')&startsWith({pet_files.name},'converted'));
    pet_path = fullfile(pet_file.folder, pet_file.name);

    % retrieve injected dose and body weight
    injected_dose_Bq = label_table.Dose(strcmp(label_table.Label,patnames{i}))* 1e6;
    body_weight_g = label_table.Weight(strcmp(label_table.Label,patnames{i}))*1000;
    gender = label_table.Gender(strcmp(label_table.Label,patnames{i}));
    body_BMI = label_table.BMI(strcmp(label_table.Label,patnames{i}));
    
   
    % Load the raw PET image (NIfTI file)
    pet_image_data = niftiread(pet_path);

    % Calculate SUV: SUV = (PET voxel value) / (Injected dose / Body weight)
     SUV_image = double(pet_image_data)* body_weight_g /injected_dose_Bq;
     vol_SUV = double(SUV_image);

    % Calculate SUV distribution per Body ROI
    featcell = cell([2 numel(mask_dir)]); 
    featcell(1,1:end) = strrep(erase({mask_dir(1:end).name}','.nii'),'_',' ');

    for j = 1:numel(mask_dir)
        disp(['Extracting from ROI: ' mask_dir(j).name '...']);
        mask = niftiread(fullfile(mask_dir(j).folder,mask_dir(j).name));

        SUV_distr = vol_SUV(mask==1);
        if size(mask,3) == size(vol_SUV,3) % check vol and mask length
            featcell{2,j} = SUV_distr;
        else
            featcell{2,j} = [];
            disp('Vol and mask not aligned!');
        end
    end
    save(fullfile(out_dir_name,patnames{i}),"featcell");
end


%% extract static SUL distribution within given ROIs (1000)
PET_dir_name = '/Volumes/Extreme Pro/SZNormal/SZWJ_1185_PET_CT_Brain_Organ_Cardiac/SZWJ_1185_PET_CT_Brain_Organ_Cardiac';
imagedir = dir(PET_dir_name);
patnames_ = {imagedir(1:end).name}';
patnames = patnames_(~startsWith(patnames_, '.')); % All patients 

label_table = readtable('/Volumes/Extreme Pro/SZNormal/Info/szwj1189.xlsx','Range','A2:M1190');
label_table.Properties.VariableNames = {'Label', 'ExaminationSeries','Name', 'Age', 'Gender','Weight', 'measured_mm', 'EmptyGlucose', 'Dose1','Dose','Height', 'BMI','categories'};

out_dir_name = '/Volumes/Extreme Pro/SZNormal/SUV_Distribution/SUL_distribution_1000';
mask_dir_name = '/Volumes/Extreme Pro/SZNormal/masks/masks_1000';

% Loop through each patient:
for i = 1:numel(patnames) 

    disp(['Working on patient ' patnames{i} '...']);
    mask_dir_ = dir(fullfile(mask_dir_name,patnames{i}));
    mask_dir = mask_dir_(~startsWith({mask_dir_.name}, '.'));

    pet_folder = fullfile(PET_dir_name, patnames{i});
    pet_files = dir(pet_folder);
    pet_file = pet_files(~startsWith({pet_files.name}, '.')&endsWith({pet_files.name},'.nii')&startsWith({pet_files.name},'converted'));
    pet_path = fullfile(pet_file.folder, pet_file.name);

      % retrieve injected dose and body weight
    injected_dose_Bq = label_table.Dose(strcmp(label_table.Label,patnames{i}))* 1e6;
    body_weight_kg = label_table.Weight(strcmp(label_table.Label,patnames{i}));
    gender = label_table.Gender(strcmp(label_table.Label,patnames{i}));
    body_BMI = label_table.BMI(strcmp(label_table.Label,patnames{i}));
    
    % Load the raw PET image (NIfTI file)
    pet_image_data = double(niftiread(pet_path));

    if strcmp(gender,'M') % for male
        LBM = 9270 * body_weight_kg / (6680 + 216 * body_BMI);
    else
        LBM = 9270 * body_weight_kg / (8780 + 244 * body_BMI);
    end

    LBM_g = LBM * 1000;

    SUL = pet_image_data * LBM_g / injected_dose_Bq;
    
    vol_SUL = double(SUL);

    % Calculate SUV distribution per Body ROI
    featcell = cell([2 numel(mask_dir)]); 
    featcell(1,1:end) = strrep(erase({mask_dir(1:end).name}','.nii'),'_',' ');

    for j = 1:numel(mask_dir)
        disp(['Extracting from ROI: ' mask_dir(j).name '...']);
        mask = niftiread(fullfile(mask_dir(j).folder,mask_dir(j).name));

        SUV_distr = vol_SUL(mask==1);
        if size(mask,3) == size(vol_SUL,3) % check vol and mask length
            featcell{2,j} = SUV_distr;
        else
            featcell{2,j} = [];
            disp('Vol and mask not aligned!');
        end
    end
    save(fullfile(out_dir_name,patnames{i}),"featcell");
end


%% extract static SUV distribution within given ROIs (pathology)
PET_dir_name = '/Volumes/Extreme Pro/糖代谢/Images/PET';

label_table = readtable('/Volumes/Extreme Pro/SZNormal/DICOM_Metadata.xlsx');
label_table = rmmissing(label_table);

out_dir_name = '/Volumes/Extreme Pro/SZNormal/SUV_Distribution/SUV_distribution_pathology';
mask_dir_name = '/Volumes/Extreme Pro/糖代谢/Images/organized_masks';

% Loop through each patient:
for i = 1:length(label_table.Label) 

    labels = label_table.Label;
    disp(['Working on patient ' labels{i} '...']);
    mask_dir_ = dir(fullfile(mask_dir_name,labels{i}));
    mask_dir = mask_dir_(~startsWith({mask_dir_.name}, '.'));

    pet_folder = fullfile(PET_dir_name, labels{i});
    pet_files = dir(pet_folder);
    pet_file = pet_files(~startsWith({pet_files.name}, '.')&endsWith({pet_files.name},'.nii')&startsWith({pet_files.name},'converted'));
    pet_path = fullfile(pet_file.folder, pet_file.name);

    % retrieve injected dose and body weight
    injected_dose_Bq = label_table.Dose(i);
    body_weight_g = label_table.Weight(i)*1000;
    gender = label_table.Sex(i);
    body_BMI = label_table.BMI(i);
    

    pet_image_data = niftiread(pet_path);   % Load the raw PET image (NIfTI file)

    % Calculate SUV: SUV = (PET voxel value) / (Injected dose / Body weight)
     SUV_image = double(pet_image_data)* body_weight_g /injected_dose_Bq;
     vol_SUV = double(SUV_image);

    % Calculate SUV distribution per Body ROI
    featcell = cell([2 numel(mask_dir)]); 
    featcell(1,1:end) = strrep(erase({mask_dir(1:end).name}','.nii'),'_',' ');

    for j = 1:numel(mask_dir)
        disp(['Extracting from ROI: ' mask_dir(j).name '...']);
        mask = niftiread(fullfile(mask_dir(j).folder,mask_dir(j).name));

        SUV_distr = vol_SUV(mask==1);
        if size(mask,3) == size(vol_SUV,3) % check vol and mask length
            featcell{2,j} = SUV_distr;
        else
            featcell{2,j} = [];
            disp('Vol and mask not aligned!');
        end
    end
    save(fullfile(out_dir_name,label_table.Label{i}),"featcell");
end

%%  extract static SUL distribution within given ROIs (pathology)
PET_dir_name = '/Volumes/Extreme Pro/糖代谢/Images/PET';

label_table = readtable('/Volumes/Extreme Pro/SZNormal/DICOM_Metadata.xlsx');
label_table = rmmissing(label_table);

out_dir_name = '/Volumes/Extreme Pro/SZNormal/SUV_Distribution/SUL_distribution_pathology';
mask_dir_name = '/Volumes/Extreme Pro/糖代谢/Images/organized_masks';

% Loop through each patient:
for i = 1:length(label_table.Label) 

    labels = label_table.Label;
    disp(['Working on patient ' labels{i} '...']);
    mask_dir_ = dir(fullfile(mask_dir_name,labels{i}));
    mask_dir = mask_dir_(~startsWith({mask_dir_.name}, '.'));

    pet_folder = fullfile(PET_dir_name, labels{i});
    pet_files = dir(pet_folder);
    pet_file = pet_files(~startsWith({pet_files.name}, '.')&endsWith({pet_files.name},'.nii')&startsWith({pet_files.name},'converted'));
    pet_path = fullfile(pet_file.folder, pet_file.name);

    % retrieve injected dose and body weight
    injected_dose_Bq = label_table.Dose(i);
    body_weight_kg = label_table.Weight(i);
    gender = label_table.Sex(i);
    body_BMI = label_table.BMI(i);
   
    % Load the raw PET image (NIfTI file)
    pet_image_data = double(niftiread(pet_path));

    if strcmp(gender,'M') % for male
        LBM = 9270 * body_weight_kg / (6680 + 216 * body_BMI);
    else
        LBM = 9270 * body_weight_kg / (8780 + 244 * body_BMI);
    end

    LBM_g = LBM * 1000;

    SUL = pet_image_data * LBM_g / injected_dose_Bq;
    
    vol_SUL = double(SUL);


    % Calculate SUV distribution per Body ROI
    featcell = cell([2 numel(mask_dir)]); 
    featcell(1,1:end) = strrep(erase({mask_dir(1:end).name}','.nii'),'_',' ');

    for j = 1:numel(mask_dir)
        disp(['Extracting from ROI: ' mask_dir(j).name '...']);
        mask = niftiread(fullfile(mask_dir(j).folder,mask_dir(j).name));

        SUL_distr = vol_SUL(mask==1);
        if size(mask,3) == size(vol_SUL,3) % check vol and mask length
            featcell{2,j} = SUL_distr;
        else
            featcell{2,j} = [];
            disp('Vol and mask not aligned!');
        end
    end
    save(fullfile(out_dir_name,label_table.Label{i}),"featcell");
end