function out = CR7T_convert_vec_BN202_to_BN246(vec)
% Load labels
labels = [];
load('workspace_labels.mat','T');
labels = T.Index_new;
ind_good = find(labels~=0);
for ii = 1:length(labels)
    if ismember(ii,ind_good)
        ind = labels(ii);
        out(ii,1) = vec(ind);
    else
        out(ii,1) = 0;
    end
end

