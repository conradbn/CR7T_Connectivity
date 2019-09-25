function output = CR7T_automask_with_clfrac(runs,subject_label,results_dir,clfrac)
%% Create EPI mask's using 3dAutomask with default settings, no dilation
% The traditional behavior of afni_proc was to dilate the EPI mask by 1
% voxel, but now has changed. Using no dilation should represent a tighter,
% more accurate mask of the EPI data. This script gets an intersection
% across all runs.

% Set subject ID and go to results (afni_proc) directory
subj = subject_label;
cd(results_dir)

%% ================================== mask ==================================
% create 'full_mask' dataset (union mask)
for run = 1:numel(runs)
    unix(['3dAutomask -clfrac ' num2str(clfrac,2)...
          ' -prefix rm.mask_r' runs{run} ' pb04.' subj '.r' runs{run} '.blur+tlrc']);
end

% Make output directory if it doesn't exist
if exist('Automask','dir') == 0
    mkdir('Automask');
end

% Remove new mask if it exists
unix(['rm -f Automask/full_mask.nodilate.clfrac' num2str(clfrac,'%.2f') '.' subj '.nii.gz']);

% create union of inputs, output type is byte
unix(['3dmask_tool -inputs rm.mask_r*+tlrc.HEAD -union'...
      ' -prefix Automask/full_mask.nodilate.clfrac' num2str(clfrac,'%.2f') '.' subj '.nii.gz']);

% remove temp files
unix('rm -f rm*');