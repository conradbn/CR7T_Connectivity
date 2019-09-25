%% CR7T Connectivity Project - Create stimulus onset files script
% This script runs routine to create stimulus onset files for use in AFNI
% programs (i.e. onset time in seconds from start of run, with separate
% runs corresponding to each line of file)

cd(dir_nii)

%% Create All Stimulus Onset File
% This "allstims" file is to be used for beta-series estimation, i.e. not 
% used by afni_proc.py script
stimvars = load(stim_file(1).name);
onsets_filename = [onset_dir 'onsets_allstims_' subject_label '.txt'];
runs = {'identify_A','identify_B','compare_A','compare_B'}; %fieldnames(stimvars);
for kk = 1:size(runs,2)
    onsets_sec(kk,:) = (stimvars.(runs{kk}).Stimuli0x2EOnsetTime ...
                        - stimvars.(runs{kk}).ISIfixationSTART0x2EOnsetTime)...
                        ./1000;
end
dlmwrite(onsets_filename,onsets_sec,'delimiter',' ');    

%% Condition-wise Onsets - Get indices/onsets
% Create condition-wise onset times structure from eprime data
stim_inds = CR7T_get_regressor_indices(dir_nii); % **FUNCTION**

% Load onsets from all stimuli (previously saved)  
onsets_allstims = load(onsets_filename);

% Catenate all runs for following multiplication
onsets_allstims_cat = reshape(onsets_allstims',1,[]);

onsets.identify =            stim_inds.task(1,:)   .* onsets_allstims_cat; 
onsets.compare =             stim_inds.task(2,:)   .* onsets_allstims_cat; 
onsets.symbolic =            stim_inds.format(1,:) .* onsets_allstims_cat; 
onsets.nonsymbolic =         stim_inds.format(2,:) .* onsets_allstims_cat;
onsets.IS =          stim_inds.task(1,:) .* stim_inds.format(1,:)  .* onsets_allstims_cat; 
onsets.IN =          stim_inds.task(1,:) .* stim_inds.format(2,:)  .* onsets_allstims_cat; 
onsets.CS =          stim_inds.task(2,:) .* stim_inds.format(1,:)  .* onsets_allstims_cat; 
onsets.CN =          stim_inds.task(2,:) .* stim_inds.format(2,:)  .* onsets_allstims_cat; 
% onsets.identify_correct =    stim_inds.task(1,:)   .* onsets_allstims_cat; 
% onsets.compare_correct =     stim_inds.task(2,:)   .* onsets_allstims_cat; 
% onsets.symbolic_correct =    stim_inds.format(1,:) .* onsets_allstims_cat; 
% onsets.nonsymbolic_correct = stim_inds.format(2,:) .* onsets_allstims_cat; 
% onsets.errors =              stim_inds.err(1,:)    .* onsets_allstims_cat; 
% onsets.omissions =           stim_inds.omi(2,:)    .* onsets_allstims_cat; 

%% Condition-wise onsets - Write files
% Get all condition names
conditions_all = fieldnames(onsets)';

% Separate onset time vectors back into runs
for kk = 1:size(conditions_all,2)
    temp_cond = onsets.(conditions_all{kk});
    temp_cond_rows = vertcat(temp_cond(1,1:80),...
                             temp_cond(1,81:160),...
                             temp_cond(1,161:240),...
                             temp_cond(1,241:320)); 

    % If stimulus does not occur during run, put a -1 placeholder for 3dDeconvolve
    for jj = 1:4
        if sum(temp_cond_rows(jj,:)) == 0
            temp_cond_rows(jj,1) = -1; 
        end
    end

    % Make zeros NaN
    temp_cond_rows(temp_cond_rows == 0) = NaN;

    % Save temporary stimulus onset file
    dlmwrite('temp.txt',temp_cond_rows,'delimiter',' '); 

    % Remove NaN's from files using unix commands, send to final file
    stim_filename = [onset_dir 'onsets_' conditions_all{kk} '_' subject_label '.txt'];
    unix('sed -e "s/-1/*/g" temp.txt >temp2.txt');
    unix(['tr -s "NaN" " " <temp2.txt >' stim_filename]);
    unix('rm -f temp*');
end
