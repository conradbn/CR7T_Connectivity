%% Check EPI coverage across Brainnetome regions
purge
cd /mnt/CR7T_Connectivity/MRI_proc
% Load Brainnetome atlas
atlas = '/mnt/CR7T_Connectivity/atlas/BN_Atlas_246_1mm_resample.nii.gz';
atlas_nii = load_untouch_nii(atlas);
atlas_img = atlas_nii.img;
rois = unique(atlas_img);

% Get subject mask names
fnames = subdir('full_mask.nodilate.clfrac0.32*.nii.gz'); % *HARD CODED string

%% Load subject masks
for ii = 1:numel(fnames)
    disp(fnames(ii).name);
    mask_nii = load_untouch_nii(fnames(ii).name);
    mask_img = mask_nii.img;
    for jj = 2:numel(rois)
        index = rois(jj);
        roi_img = atlas_img;
        roi_img(roi_img ~= index) = 0;
        roi_img(roi_img == index) = 1;
        roi_size = sum(roi_img(:));
        mask_roi = double(mask_img) .* double(roi_img);
        mask_included_size = sum(mask_roi(:));
        %percentage_coverage(ii,jj) = mask_included_size/roi_size *100;
        percentage_coverage(ii,jj) = round(mask_included_size/roi_size *100);
    end
end

%% Plot percentage coverage
figure
imagesc(percentage_coverage(:,2:end)); colorbar;

%% Get good and bad regions
cutoff = 50;
mcutoff = 50;
%ind_good = find(min(percentage_coverage)>=cutoff & median(percentage_coverage)>mcutoff);
ind_good = find(min(percentage_coverage)>=cutoff);
disp(['Good ROIs = ' num2str(numel(ind_good))]);
%ind_bad = find(min(percentage_coverage)<cutoff | median(percentage_coverage)<=mcutoff);
ind_bad = find(min(percentage_coverage)<cutoff);
disp(['Bad ROIs = ' num2str(numel(ind_bad))]);

close all
figure;
boxplot(percentage_coverage(:,ind_good));
rois_good = rois(ind_good);
figure;
boxplot(percentage_coverage(:,ind_bad));
rois_bad = rois(ind_bad);

%% Create new atlas to visualize the dropout
output_prefix = '/mnt/CR7T_Connectivity/atlas/BN_Atlas_1mm_resample_reduced_clfrac0.32'; % *HARD CODED string
unix(['rm -f ' output_prefix '*']);
% First write indices to text files
rois_final = rois;
rois_final(ind_bad) = 0;
rois_final = rois_final(2:end);

dlmwrite([output_prefix '.txt'],rois_final,'delimiter',' ');
  
% Create new map of roi stats
unix(['3dUndump -datum float -ROImask ' atlas...
      ' -prefix ' output_prefix '.nii.gz '...
      output_prefix '.txt']); 
  
%% %%%%%%%%%%%%%%%%% REMOVE SUBJECT %%%%%%%%%%%%%%%%%%%%%%%%%%%
percentage_coverage2 = percentage_coverage;
percentage_coverage2(6,:) = [];
% Get good and bad regions
cutoff = 50;
mcutoff = 50;
%ind_good = find(min(percentage_coverage2)>=cutoff & median(percentage_coverage2)>mcutoff);
ind_good = find(min(percentage_coverage2)>=cutoff);
disp(['Good ROIs = ' num2str(numel(ind_good))]);
%ind_bad = find(min(percentage_coverage2)<cutoff | median(percentage_coverage2)<=mcutoff);
ind_bad = find(min(percentage_coverage2)<cutoff);
disp(['Bad ROIs = ' num2str(numel(ind_bad))]);

close all
figure;
boxplot(percentage_coverage2(:,ind_good));
rois_good = rois(ind_good);
figure;
boxplot(percentage_coverage2(:,ind_bad));
rois_bad = rois(ind_bad);

% Create new atlas to visualize the dropout
output_prefix = '/mnt/CR7T_Connectivity/atlas/BN_Atlas_1mm_resample_reduced_clfrac0.32_rm_CR_006'; % *HARD CODED string
unix(['rm -f ' output_prefix '*']);
% First write indices to text files
rois_final = rois;
rois_final(ind_bad) = 0;
rois_final = rois_final(2:end);

dlmwrite([output_prefix '.txt'],rois_final,'delimiter',' ');
  
% Create new map of roi stats
unix(['3dUndump -datum float -ROImask ' atlas...
      ' -prefix ' output_prefix '.nii.gz '...
      output_prefix '.txt']);   


% Create new atlas to visualize the dropout - FOR REGIONS EXCLUDED
output_prefix = '/mnt/CR7T_Connectivity/atlas/BN_Atlas_1mm_resample_reduced_clfrac0.32_rm_CR_006_excluded'; % *HARD CODED string
unix(['rm -f ' output_prefix '*']);
% First write indices to text files
rois_final_exclude = rois;
rois_final_exclude(ind_good) = 0;
rois_final_exclude = rois_final_exclude(2:end);

dlmwrite([output_prefix '.txt'],rois_final_exclude,'delimiter',' ');
  
% Create new map of roi stats
unix(['3dUndump -datum float -ROImask ' atlas...
      ' -prefix ' output_prefix '.nii.gz '...
      output_prefix '.txt']);   
  
%% Save workspace
save('/mnt/CR7T_Connectivity/matfiles/workspace_region_coverage');

