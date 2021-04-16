function [C_all,M_all,Q_all,Pa_all,Pa_group,C_group] = CR7T_louvain_with_consensus_clustering_group_partition(mats_struct,gammas,tau,niters)
%% MODULARITY ANALYSIS 
% 10/20/18
% Calculate community partiion for each subject/condition using Louvain
% modularity maximization with consensus clustering across niters for each subject.
% Then group level consensus clustering of group agreement matrix

% Check if structure, get/set condition names
if isstruct(mats_struct)
    conds = fieldnames(mats_struct);
else
    disp('Warning - Matrices not in a structured variable, using "mats" as a prefix')
    conds = {'mats'};
    ms = mats_struct;
    clear mats_struct;
    mats_struct.mats = ms;
end

%% Subject Level - Louvain partitioning plus consensus clustering 
% For each condition
for ii = 1:size(conds)
    mats = mats_struct.(conds{ii});
    % Replace Inf and NaN values
    mats(isnan(mats)) = 0;
    mats(isinf(mats)) = 0;
    % For all gamma values
    for jj = 1:numel(gammas)
        % For all subjects
        for kk = 1:size(mats,3)
            disp(['Subject ' num2str(kk) ' of ' num2str(size(mats,3))...
                  ' - ' conds{ii} ' - Gamma ' num2str(gammas(jj))]);
            % Run Louvain algorithm, with subject-level consensus clustering 
            % over n iterations
            [M,Mi,Qi,Pa] = CR7T_community_louvain_with_iteration(mats(:,:,kk),gammas(jj),tau,[],'negative_asym',niters);
            % Communisty (module) index
            C_all.(conds{ii})(:,jj,kk) = M;
            M_all.(conds{ii})(:,:,jj,kk) = Mi;
            Q_all.(conds{ii})(:,jj,kk) = Qi;
            Pa_all.(conds{ii})(:,:,jj,kk) = Pa;
        end
    end
end

%% Group level - consensus clustering on percent agreement matrix
% Get number of subjects
ns = size(mats,3);
% For each condition 
for ii = 1:size(conds)
    disp(conds{ii});
    % For each gamma value
    for jj = 1:numel(gammas)
        Ci = squeeze(C_all.(conds{ii})(:,jj,:));
        % Percent agreement matrices (have to divide by number of subjects)
        P = agreement(Ci)./ns;
        Pa_group.(conds{ii})(:,:,jj) = P;
        % Consensus partition across group
        C_group.(conds{ii})(:,jj) = consensus_und(P,tau,niters);
    end
end

