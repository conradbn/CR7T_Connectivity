function [M,Mi,Qi,Magree] = CR7T_community_louvain_with_iteration(W,gamma,tau,M0,B,iters)

% Run Louvain community detection algorithm many times (n = iters)
for ii = 1:iters
    [Mi(:,ii),Qi(ii)] = community_louvain(W,gamma,M0,B);
end

% Calculate agreement matrix of all partitioning
Magree = agreement_weighted(Mi, Qi);
%Magree = agreement(Mi);

% Convert to probability *NOT NEEDED IF USING agreement_weighted
% Mprob = Magree./iters;

% Use agreement matrix to provide consensus partition for the given matrix
M = consensus_und(Magree,tau,iters);

