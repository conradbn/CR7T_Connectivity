%% Price 7T Cortical Representation Project - Preprocessing Pipeline 
%% Preprocessing Script
% by Benjamin Neal Conrad, January 2019
% Vanderbilt University

% Now with correct slice-timing information
% afni_proc now runs the 3dREML activation analysis, instead of
% beta-series estimation (since I do beta-series estimation now via 3dLSS)

%% Start process timer
tic

%% Clear workspace, add matfiles and Python, AFNI binary files to path
%     clear variables; close all;
%     addpath(genpath('/Users/benconrad/Desktop/Price7T/matfiles/'));
add_AFNI_to_path;

%% Subject label and directories
subject_dir_raw = pwd;
[upperPath, deepestFolder, ~] = fileparts(subject_dir_raw);
subject_label = deepestFolder;

% Check to be sure you are in a subject's raw directory
if contains(subject_dir_raw,'/MRI_raw/CR_0') ~= 1
    disp('You are not in a subject''s raw directory! This script cannot proceed.');
    return
end 

% Create subject nifti (nii) directory
% dir_nii = strrep(subject_dir_raw,'raw','nii'); 
temp_dir = '/media/sf_afni_tmp/';
dir_nii = [temp_dir subject_label]; 
unix(['rm -rf ' dir_nii '_AProc']); % remove if already exists
mkdir([dir_nii '_AProc']);
dir_nii = [dir_nii '_AProc/'];

% Save matlab command window ouput
diary([dir_nii 'command_window_output_' subject_label '.txt']); diary on;

% Make directory for "original" nifti files - orig_dir
orig_dir = [dir_nii 'orig/'];
mkdir(orig_dir);

% Make directory for onset files - onset_dir
onset_dir = [dir_nii 'onsets/'];
mkdir(onset_dir); 

%% Convert PARREC files to nifti format
% Anatomical(s) (MP2RAGE)
anatomical_parrecs = dir('*WIP_MP2RAGE_1x1x1_SENSE*.PAR');
for ii = 1:size(anatomical_parrecs,1)
    tmp_name = anatomical_parrecs(ii).name;
    tmp_nii_name = strrep(anatomical_parrecs(ii).name,'PAR','nii');
    % Don't do it if file already exists
    if exist([orig_dir tmp_nii_name],'file') ~= 2
        dicm2nii(tmp_name,orig_dir,'nii');
    end
end

% Functional (fMRI) Runs
functional_parrecs = dir('*WIP_fMRI_2o5mm_iso_SENSE*.PAR');
for ii = 1:size(functional_parrecs,1)
    tmp_name = functional_parrecs(ii).name;
    tmp_nii_name = strrep(functional_parrecs(ii).name,'PAR','nii');
    % Don't do it if file already exists
    if exist([orig_dir tmp_nii_name],'file') ~= 2
        dicm2nii(tmp_name,orig_dir,'nii');
    end
end

% Stimulus matlab file
stim_file = dir('*_stimvars+rt+acc.mat');
copyfile(stim_file(1).name,dir_nii);

% Clean up workspace
clear ans ii tmp_nii_name tmp_name deepestFolder ...
        anatomical_parrecs B0_parrecs functional_parrecs; 

%% MP2RAGE Correction
% This step implements the correction procedure described in Marques et
% al., 2010 NeuroImage, combining the two volumes within the MP2RAGE to
% create a corrected version which is more purely T1 weighted. The first is
% a standard MPRAGE acquisition (TI1) and the second is exactly the same but
% is acquired after waiting a longer time (TI2) after the initial adiabatic
% inversion pulse. This process eliminates the non-T1w effects from B1
% inhomogeneity, proton density (M0), and T2* which are present in the TI1
% image.

cd(orig_dir)
% Load 4D volume (contains TI1 & TI2)
anat_nii_name = dir('*_WIP_MP2RAGE_1x1x1_SENSE*1_magnitude.nii');
anat_nii_name = anat_nii_name(1).name;
anat_nii = load_untouch_nii(anat_nii_name);

% Get separate 3D volumes
anat_TI1 = double(anat_nii.img(:,:,:,1)) ;%./ 4.2310e-04;% scale factor
anat_TI2 = double(anat_nii.img(:,:,:,2)) ;%./ 4.2310e-04;% scale factor

% Combination equation
beta = 0; %.5e04;
anat_TIfinal = ((anat_TI1.*anat_TI2) - beta)./ ...
               ((anat_TI1.^2 + anat_TI2.^2) + beta*2) ;

% Plot results
plot_mp2rage(anat_TI1,anat_TI2,anat_TIfinal,130,100,130);
print([orig_dir 'MP2RAGE_correction_' subject_label],'-dpng');

% Create individual anatomical volumes
unix(['3dTcat -prefix ' orig_dir 'anat_orig_' anat_nii_name ' ' anat_nii_name '[0]']);
unix(['3dTcat -prefix ' orig_dir 'anat_pd_' anat_nii_name ' ' anat_nii_name '[1]']);

% Create corrected anatomical volume
anat_pd_nii_name = dir('anat_pd_*'); anat_pd_nii_name = anat_pd_nii_name(1).name;
anat_pd_nii = load_untouch_nii(anat_pd_nii_name);
anat_corr_nii = anat_pd_nii; anat_corr_nii.img = anat_TIfinal;

save_untouch_nii(anat_corr_nii,[orig_dir 'anat_corr_' anat_nii_name]);

%% Skullstripping of Anatomical
% Will use TI2 MP2RAGE (mostly PDw) image as it provides the nicest outline of brain due
% to the proton density weighting (i.e. contrast of skull/tissue)

% First normalize intensities in PD volume
unix(['3dUnifize -T2 -input ' orig_dir 'anat_pd_' anat_nii_name ...
      ' -prefix anat_pdnorm_' anat_nii_name]);

% Run skullstripping on normalized (bias corrected) PD volume
unix(['3dSkullStrip -shrink_fac_bot_lim 0.7 -mask_vol -touchup -touchup'...
      ' -input ' orig_dir 'anat_pdnorm_' anat_nii_name ...
      ' -prefix anat_bmask_alledges_' anat_nii_name]); % -push_to_edge -shrink_fac 0.8
% Run skullstripping on original PD volume (necessary for some subjects)
% unix(['3dSkullStrip -shrink_fac_bot_lim 0.7 -mask_vol -touchup -touchup'...
%       ' -input ' orig_dir 'anat_pd_' anat_nii_name ...
%       ' -prefix anat_bmask_alledges_' anat_nii_name]);

% Create mask with only inner voxels 
unix(['3dcalc -a anat_bmask_alledges_' anat_nii_name ' -expr "ispositive(a-5)"'...
      ' -prefix anat_bmask_' anat_nii_name]); 

% Apply brain mask to anatomicals
unix(['3dcalc -a ' orig_dir 'anat_orig_' anat_nii_name ...
            ' -b ' orig_dir 'anat_bmask_' anat_nii_name ...
            ' -prefix ' orig_dir 'k_anat_orig_' anat_nii_name ...
            ' -expr a*b']);

unix(['3dcalc -a ' orig_dir 'anat_corr_' anat_nii_name ...
            ' -b ' orig_dir 'anat_bmask_' anat_nii_name ...
            ' -prefix ' orig_dir 'k_anat_corr_' anat_nii_name ...
            ' -expr a*b']);

unix(['3dcalc -a ' orig_dir 'anat_pd_' anat_nii_name ...
            ' -b ' orig_dir 'anat_bmask_' anat_nii_name ...
            ' -prefix ' orig_dir 'k_anat_pd_' anat_nii_name ...
            ' -expr a*b']);

% Load brainmasked (k), corrected MP2RAGE
k_anat_corr_nii_name = dir('k_anat_corr_*'); k_anat_corr_nii_name = k_anat_corr_nii_name(1).name;
k_anat_corr_nii = load_untouch_nii(k_anat_corr_nii_name);

% Plot results
plot_mp2rage(anat_TI1,anat_TI2,k_anat_corr_nii.img,130,100,130);
print([orig_dir 'MP2RAGE_brainmask_correction_' subject_label],'-dpng');

%% Create All Stimulus Onset File
CR7T_create_stimulus_onset_files;

% Clean up
clear kk jj temp_cond temp_cond_rows temp_inds

%% Save workspace before AFNI Proc Script
save(['workspace_variables_pre_proc_' subject_label]);  

%% AFNI Proc Script
% Create process script
unix(['afni_proc.py -subj_id ' subject_label ...
        ' -dsets ' orig_dir '*_WIP_fMRI_2o5mm_iso_SENSE*.nii'...
        ' -blocks despike tshift align tlrc volreg blur mask scale regress'...
        ' -copy_anat ' orig_dir k_anat_corr_nii_name ...
        ' -tshift_opts_ts'...
        '    -verbose'...
        '    -tpattern alt+z'...
        ' -align_opts_aea'...
        '    -cost lpc+ZZ'...
        '    -big_move'...
        ' -anat_has_skull no'...
        ' -align_unifize_epi yes'...
        ' -tlrc_base ~/abin/MNI152_T1_2009c+tlrc'...
        ' -volreg_warp_dxyz 2.5'...
        ' -volreg_align_to MIN_OUTLIER'...
        ' -volreg_tlrc_warp'...
        ' -volreg_compute_tsnr yes'...
        ' -blur_size 4'...
        ' -regress_stim_times '...
            onset_dir '*IS*.txt '...
            onset_dir '*IN*.txt '...  
            onset_dir '*CS*.txt '...  
            onset_dir '*CN*.txt '...
        ' -regress_stim_labels IS IN CS CN'...
        ' -regress_basis "GAM(8.6,.547,.5)"'...
        ' -regress_censor_motion 0.3'...  
        ' -regress_censor_outliers 0.05'...
        ' -regress_apply_mot_types demean deriv'...
        ' -regress_opts_3dD'...
        '    -bout'...
        '    -jobs 4'...
        '    -noFDR'...
        '    -allzero_OK'...
        '    -GOFORIT 3'...
        '    -gltsym "SYM: +CS -CN"                             -glt_label 1 diff_CS_CN'...
        '    -gltsym "SYM: +0.5*CS +0.5*CN"                     -glt_label 2 mean_CS_CN'...
        '    -gltsym "SYM: +0.5*CS +0.5*CN -0.5*IS -0.5*IN"     -glt_label 3 diff_C_I'...
        '    -gltsym "SYM: +0.5*CS -0.5*CN +0.5*IS -0.5*IN"     -glt_label 4 diff_S_N'...
        '    -gltsym "SYM: +IS -IN"                             -glt_label 5 diff_IS_IN'...
        '    -gltsym "SYM: +0.5*IS +0.5*IN"                     -glt_label 6 mean_IS_IN'...
        '    -gltsym "SYM: +0.25*IS +0.25*IN +0.25*CS +0.25*CN" -glt_label 7 mean_Allstims'...
        ' -regress_3dD_stop'...
        ' -regress_reml_exec'...
        ' -regress_opts_reml'...
        '    -GOFORIT 3'...
        '    -Rwherr errts.' subject_label '_REMLwh'...
        ' -regress_run_clustsim yes'...
        ' -regress_make_ideal_sum sum_ideal.1D'...
        ' -regress_est_blur_epits'...
        ' -regress_est_blur_errts']);

% Run process script
unix(['tcsh -xef proc.' subject_label]); %|& tee output.proc.CR_040_315023

%% Compress all BRIK and NII files to save space
cd(dir_nii);
unix('find . -name \*.BRIK -exec pigz -v {} \;');
unix('find . -name \*.nii -exec pigz -v {} \;');

%% Set (afni) results directory
results_dir = [dir_nii subject_label '.results/'];

%% Fix QC html document
CR7T_fix_qc_html

%% Set final directory location on server
serv_dir = '/mnt/CR7T_Connectivity/MRI_proc/';
dir_nii_orig = dir_nii;
dir_nii = strrep(dir_nii,temp_dir,serv_dir);
orig_dir = strrep(orig_dir,temp_dir,serv_dir);
onset_dir = strrep(onset_dir,temp_dir,serv_dir);
results_dir = strrep(results_dir,temp_dir,serv_dir);

%% Save workspace variables
cd(dir_nii_orig)
% First clean up a bit
clear anat_nii anat_corr_nii k_anat_corr_nii anat_pd_nii func3_nii fnii_1 anat_TI1 anat_TI2 anat_TIfinal ans beta
save(['workspace_variables_post_proc_' subject_label]);

%% Turn off diary
toc % Report total processing time
diary off;

%% Move folder to server
mkdir(dir_nii);
disp('...Moving files to server...');
unix(['mv ' dir_nii_orig '* ' dir_nii]);



