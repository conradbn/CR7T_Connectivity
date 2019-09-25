%% Check number of stimuli used per subject/condition (after censoring)
purge
cd Z:\CR7T_Connectivity\MRI_Proc

%% Symbolic
fnames = subdir('*BN202_Compare_Symbolic_num_TRs_used.mat');
for ii = 1:numel(fnames)
    load(fnames(ii).name);
    num_stims_used(ii,1) = num_Stims_for_df;
end

%% Nonsymbolic
fnames = subdir('*BN202_Compare_Nonsymbolic_num_TRs_used.mat');
for ii = 1:numel(fnames)
    load(fnames(ii).name);
    num_stims_used(ii,2) = num_Stims_for_df;
end
   
%% Get MAD below the threshold
mad_sym = mad(num_stims_used(:,1));
mad_nonsym = mad(num_stims_used(:,2));