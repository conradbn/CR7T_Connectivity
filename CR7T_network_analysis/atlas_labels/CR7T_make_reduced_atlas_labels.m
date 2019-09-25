purge;
T = readtable('workspace_labels.csv');

excluded = find(T.Excluded == 1);
ind = 1;
for ii = 1:size(T,1)
    if ismember(ii,excluded)
        T.Index_new(ii,1) = 0;
        T.Labels_reduced(ii,1) = {''};
    else
        T.Index_new(ii,1) = ind;
        T.Labels_reduced(ii,1) = T.Region(ii);
        ind = ind+1;
    end
end

clear ii ind
save('workspace_labels.mat')