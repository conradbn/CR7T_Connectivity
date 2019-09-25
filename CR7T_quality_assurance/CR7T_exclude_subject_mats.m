%% Exclude subjects from original sample
purge
load('Z:\CR7T_Connectivity\matfiles\workspace_zscore_mats.mat');

% Exclude subjects
subj_dirs_orig = subj_dirs;
subs_exclude = [5,6,8,11,12,28,32];
subj_dirs(subs_exclude) = [];
zscore_mats.BN202_Compare_All(:,:,subs_exclude) = [];
zscore_mats.BN202_Compare_Symbolic(:,:,subs_exclude) = [];
zscore_mats.BN202_Compare_Nonsymbolic(:,:,subs_exclude) = [];
zscore_mats_prenorm.BN202_Compare_All(:,:,subs_exclude) = [];
zscore_mats_prenorm.BN202_Compare_Symbolic(:,:,subs_exclude) = [];
zscore_mats_prenorm.BN202_Compare_Nonsymbolic(:,:,subs_exclude) = [];

% Save
save('Z:\CR7T_Connectivity\matfiles\workspace_zscore_mats_33subs.mat');