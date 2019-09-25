%% Plot all subject matrices
purge
load('Z:/CR7T_Connectivity/matfiles/CR7T_WORKSPACES/workspace_zscore_mats_33subs.mat');

%zscore_mats = zscore_mats_prenorm;
condition_names = fieldnames(zscore_mats);

% Load RdBu color map
RdBu = cbrewer('div', 'RdBu', 256,'PCHIP');
RdBu = flip(RdBu);

%% Display the Zscore matrices 
for cond = 1:numel(condition_names)
    % Load condition mats, remove inf/nan vals
    zmats = zscore_mats.(condition_names{cond});
    zmats(isnan(zmats)) = 0;
    zmats(isinf(zmats)) = 0;
    % Set figure window size for better looking display
    figure('Position', [100, 100, 1700, 1000]);  % 1700,1000
    % Get condition name
    condsplit = strsplit(condition_names{cond},'_');
    condsplit = condsplit{1,end};
    % Plot matrices
    for jj = 1:size(zmats,3)
        %subaxis(5,8,jj,'Spacing', 0.02, 'Padding', 0, 'Margin', 0.05); 
        subaxis(5,7,jj,'Spacing', 0.02, 'Padding', 0, 'Margin', 0.05); 
        imagesc(zmats(:,:,jj)); colormap(RdBu);
        caxis([-6 6])
        axis square;
        axis off; 
        subj_name = subj_dirs(jj).name;
        title(['Zscores CR-' subj_name(4:6) ' ' condsplit ]);
%         if jj == size(zmats,3)
%             colorbar
%         end
    end
end

%% Plot the subject level difference matrices
% Load condition mats, remove inf/nan vals
zmats = zscore_mats.BN202_Compare_Symbolic - zscore_mats.BN202_Compare_Nonsymbolic;
zmats(isnan(zmats)) = 0;
zmats(isinf(zmats)) = 0;
% Set figure window size for better looking display
figure('Position', [100, 100, 1700, 1000]);  % 1700,1000
% Get condition name
condsplit = strsplit(condition_names{cond},'_');
condsplit = condsplit{1,end};
% Plot matrices
for jj = 1:size(zmats,3)
    %subaxis(5,8,jj,'Spacing', 0.02, 'Padding', 0, 'Margin', 0.05); 
    subaxis(5,7,jj,'Spacing', 0.02, 'Padding', 0, 'Margin', 0.05); 
    imagesc(zmats(:,:,jj)); colormap(RdBu);
    caxis([-3 3])
    axis square;
    axis off; 
    subj_name = subj_dirs(jj).name;
    title(['Z-Diff (Sym-Nonsym) CR-' subj_name(4:6)]);
%         if jj == size(zmats,3)
%             colorbar
%         end
end

%% Show the mean matrices
for cond = 1:numel(condition_names)
    zmats = zscore_mats.(condition_names{cond});
    zmats(isnan(zmats)) = 0;
    zmats(isinf(zmats)) = 0;
    zscore_mats_mean(:,:,cond) = mean(zmats,3);
end
figure('Position', [100, 100, 1200, 400]);
for cond = 1:numel(condition_names)
    subaxis(1,3,cond,'Spacing', 0.02, 'Padding', 0, 'Margin', 0.05);
    imagesc(zscore_mats_mean(:,:,cond)); colormap(RdBu); 
    caxis([-6 6])
    axis square;
    axis off; 
    title(['Mean Zscores ' condition_names{cond}],'interpreter','none');
%     if cond == numel(condition_names)
%         colorbar
%     end
end

%% Diff matrix
diff = zscore_mats_mean(:,:,2)-zscore_mats_mean(:,:,3);
figure
imagesc(diff(:,:)); colormap(RdBu); colorbar
caxis([-1 1])
axis square;
axis off; 
title('Zscores Difference: Symbolic - Nonsymbolic','interpreter','none');



