%% Calculate histograms of voxelwise mean signal intensity
% across 3dAutomask clip fractions. Then get group mean distributions and plot

purge
cd /mnt/CR7T_Connectivity/MRI_proc
fnames = subdir('workspace_variables_post*.mat');
clfracs = [0.1,0.2,0.25,0.3,0.31,0.32,0.33,0.34,0.35,0.4,0.45,0.5];

%% Loop through all clipfractions and subjects
for cf = 1:numel(clfracs)
    clfrac = clfracs(cf);
    for ii = 1:numel(fnames)
        load(fnames(ii).name,'subject_label','results_dir');
        subj = subject_label;
        cd(results_dir)
        disp(results_dir);
        %% Create overall mean image for each run, if needed
        if ~exist(['mean.pb04.' subj '.r04.blur.nii.gz'],'file')
            runs = {'01','02','03','04'};
            for run = 1:numel(runs)
                unix(['3dTstat -prefix mean.pb04.' subj '.r' runs{run} '.blur.nii.gz pb04.' subj '.r' runs{run} '.blur+tlrc']);
            end
        end
        %% Create global mean image, if needed
        if ~exist(['mean.pb04.' subj '.all.blur.nii.gz'],'file')
            unix(['3dMean -prefix mean.pb04.' subj '.all.blur.nii.gz'...
                    ' mean.pb04.' subj '.r01.blur.nii.gz'...
                    ' mean.pb04.' subj '.r02.blur.nii.gz'...
                    ' mean.pb04.' subj '.r03.blur.nii.gz'...
                    ' mean.pb04.' subj '.r04.blur.nii.gz']);
        end
        % Now calculate the histogram
        data = ['mean.pb04.' subj '.all.blur.nii.gz'];
        mask = subdir(['full_mask.nodilate.clfrac' num2str(clfrac,'%.2f') '.*.nii.gz']);
        mask = mask.name;
        [v,c] = CR7T_automask_calculate_histogram(data,mask,100);
        %% Normalize since some subjects have order of magnitude larger values in epi
        vals(:,ii,cf) = normalize(v); 
        counts(:,ii,cf) = normalize(c);
    end
end
%% Save workspace
save('/mnt/CR7T_Connectivity/matfiles/workspace_automask_histograms');

%% Plotting
close all

% Get mean distributions
for cf = 1:numel(clfracs)
    mcounts(:,cf) = mean(counts(:,:,cf),2);
    mvals(:,cf) = mean(vals(:,:,cf),2);
end

% Find local minima
local_min = islocalmin(mcounts);

for cf = 1:numel(clfracs)
    figure('Position',[100,100,800,400]);
    plot(vals(:,:,cf),counts(:,:,cf))
    hold on
    plot(mvals(:,cf),mcounts(:,cf),'-k','LineWidth',3)
    if any(local_min(:,cf) == 1)
        loc = find(local_min(:,cf) == 1);
        line([mvals(loc) mvals(loc)], [min(mcounts(:,cf)) max(mcounts(:,cf))],'linewidth',2,'color','r');
    end
    grid on
    xlabel('Voxelwise signal intensity (arbitrary units)')
    ylabel('Voxel count (arbitrary units)')
    title(['Voxelwise signal intensity distribution - 3dAutomask clip fraction = ' num2str(clfracs(cf))]);
    % Print
    print('-dpng','-r300',['/mnt/CR7T_Connectivity/figures/signal_distribution_3dAutomask_clfrac_' num2str(clfracs(cf),'%.2f') '.png']);
end
