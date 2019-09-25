%% Group-level consensus clustering, permuations
% 10/29/18
% Run group-level consensus clustering on pair-permuted subject-level partitions
% Can use output for testing differences between Symbolic and Nonsymbolic conditions
% in true wholebrain/module-level allegiance values, as well as other metrics
% at the group level like consistency and diversity

tic
purge;
% Shuffle the random number generator to get a unique set of permutations
rng('shuffle');

% Start up the parallel pool of CPU "workers"
%parpool

% Set number of permutations to run
num_perms = 50000; 

% Load subject agreement matrices 
load('workspace_consensus_clustering_gamma_2pt45','C_all');
cnames = fieldnames(C_all);

% Set community index mats for permuations
C_perms(:,:,1) = squeeze(C_all.BN202_Compare_Symbolic);
C_perms(:,:,2) = squeeze(C_all.BN202_Compare_Nonsymbolic);
num_subs = size(C_perms,2); % Number of subjects

% Convert to single to save space
C_perms = single(C_perms);

%% First get random relabeling of each community-index 
for p = 1:num_perms
    t = randi([1 2],num_subs,1);
    perm_inds(:,p) = single(t);
end

% Permutations with forced balancing, i.e. 50% swapped conditions 
half_perms = num_perms/2;
for p = 1:num_perms/2
    t = randi([1 2],num_subs,1);
    perm_inds(:,p) = single(t);
end
for p = 1:(num_perms/2) + 1
    t = randi([1 2],num_subs,1);
    perm_inds(:,p) = single(t);
end

% Create pair-permuted community partition (C) vectors
for p = 1:size(perm_inds,2)
    % Get subject index/opposite index for this permutation
    inds_cond1 = perm_inds(:,p);
    inds_cond2 = inds_cond1-1;
    inds_cond2(inds_cond2 == 0)=2; 
    % Use index vectors to extract subject community partitions for this
    % permutation
    for s = 1:num_subs
        C_all_perm(:,s,1,p) = C_perms(:,s,inds_cond1(s));
        C_all_perm(:,s,2,p) = C_perms(:,s,inds_cond2(s));
    end
end
   
%% Run consensus clustering across permutations
% pag = zeros(size(C_perms,1),size(C_perms,1),num_perms,2); % Actually
% slower if preallocated..
% For each "condition" - this helps reduce overhead to cpu nodes
for ii = 1:2
    % Progress bar
    parfor_progress(num_perms);
    ca = squeeze(C_all_perm(:,:,ii,:));
    % Loop through subjects using parallel computations
    parfor p = 1:num_perms
        % Subject-level consensus partitions (pair permuted)
        c = ca(:,:,p);
        % Percent agreement matrix, across subject
        pa = agreement(c)./num_subs;
        % Consensus partition, across group
        cg(:,p,ii) = single(consensus_und(pa,0.5,50));
        % Percent agreement matrices across group
        pag(:,:,p,ii) = pa;
        % Print current progress
        parfor_progress;
    end
end

% Get the final group consensus vectors and percent agree mats from pair permutation
C_group_perm = cg;
Pa_group_perm = pag;

%% Save data
clear ans c ca cg pa pag p s t ii C_perms inds_cond1 inds_cond2
save(['Z:\CR7T_Connectivity\matfiles\CR7T_network_analysis\permutations_group_level\'...
      'workspace_group_consensus_pair_permutations_' num2str(num_perms) 'iters_' datestr(now,'mm_DD_HH_MM_SS')],'-v7.3');
% Report total compute time
toc








