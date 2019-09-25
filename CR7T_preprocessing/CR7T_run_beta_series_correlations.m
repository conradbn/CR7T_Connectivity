function [zscore_mat,zscore_mat_prenorm] = CR7T_run_beta_series_correlations(subj_dir, beta_stim_index_function, beta_str_tag, atlas, atlas_labels,mask_prefix)
%% Run Beta Series Correlations
% Process to first extract average beta series from regions of interest for the
% specified condition trials, then compute the cross-correlation matrix and
% convert to z-scores, as well as printing out of matrices.
%
% 10/17/18 - Benjamin Conrad
%
% -Uses the output from 3dLSS beta estimation on final blurred/scaled data, 
%  which included WM and CSF principle components with other nuisance regressors. 
%
% -This script also implements a normalization procedure (within condition)
%  to reduce the potential for mean activity level differences to influence
%  connectivity results; But both original and normalized versions are
%  output

cd(subj_dir);

% Get subject stimulus regressors and specific indices for requested stimuli
regressors = CR7T_get_regressor_indices(pwd); % **Custom function** 
beta_stim_index = eval(beta_stim_index_function);

% Get necessary subject variables
variables = dir('workspace_variables_post*.mat');
variables = variables.name; 
load(variables,'subject_label','results_dir');

% Check if beta directory exists
beta_dir = strrep(results_dir,'.results','.beta_series_correlation');
if exist(beta_dir,'dir') == 0
    mkdir(beta_dir);
end

% Remove results sub-directory if it exists
unix(['rm -rf ' beta_dir beta_str_tag]);

% Create results sub-directory
mkdir([beta_dir beta_str_tag]);
beta_dir_sub = [beta_dir beta_str_tag '/'];
beta_label = [subject_label '_' beta_str_tag];

%% Create new out.ss_review text file - corresponding to beta-series estimation
% cd(results_dir)
% % Only do if needed (i.e. the first time)
% if exist([subject_label '_beta_series.out.ss_review.txt'],'file') == 0
%     % Generate scripts to review single subject results
%     unix(['gen_ss_review_scripts.py '...
%           ' -prefix ' results_dir subject_label '_beta_series.'...
%           ' -xmat_regress ' results_dir subject_label '_beta_series.X.xmat.1D'... 
%           ' -mot_limit 0.3'...
%           ' -out_limit 0.05']);
%             %' -xmat_uncensored ' prefix '_beta_series.uncensored.X.xmat.1D'...
% 
%     % Run basic subject review script
%     unix(['./' subject_label '_beta_series.@ss_review_basic '...
%           '|& tee ' results_dir subject_label '_beta_series.out.ss_review.txt']);
% end

%% Open afni_proc "review" txt output file to get number of censored TRs per stimulus trials
fid  = fopen([results_dir subject_label '_beta_series.out.ss_review.txt'],'r');
text = textscan(fid,'%s','Delimiter','');
text = text{1};
fid  = fclose(fid);

% Parse file for string of interest
text = string(text(contains(text,'num TRs censored per stim :')));
% Remove variable name string
text = regexprep(text,'num TRs censored per stim :','');
% Get values for each TR
num_TRs_censored = sscanf(text,'%f');
% Get index of trials to be censored
num_TRs_censored_idx = find(num_TRs_censored > 1);

%% Create Beta series 4D volume
% Go to subject beta sub-directory
cd(beta_dir_sub)
% First get only the indexes to be kept, removing censored trials
keep_idx = setdiff(beta_stim_index,num_TRs_censored_idx);
% Make it a string, to pass to AFNI command
keep_idx_str = regexprep(num2str(keep_idx-1),'\s+',',');
% Set the beta-series volume
betas_nii = [results_dir subject_label '_beta_series_LSS.nii.gz'];
% Create censored beta-series volume
unix(['3dcalc -prefix ' beta_label '_betamaps_censored.nii.gz -a ' betas_nii '[' keep_idx_str '] -expr "a"']);
% Perform normalization (voxelwise subtraction of mean and division by std)
unix(['3dTstat -mean -prefix tmp_mean.nii ' beta_label '_betamaps_censored.nii.gz']);
unix(['3dTstat -stdev -prefix tmp_std.nii ' beta_label '_betamaps_censored.nii.gz']);
unix(['3dcalc -prefix ' beta_label '_betamaps_censored_condnorm.nii.gz'...
      ' -a ' beta_label '_betamaps_censored.nii.gz'...
      ' -b tmp_mean.nii -c tmp_std.nii '...
      ' -expr "(a-b)/c"']);
unix('rm -f tmp*.nii');

%% ****** Beta correlations AFTER normalization *******************************************************************
% *************************************************************************
%% Mask the beta volume before extracting series from ROIs
unix(['3dcalc -prefix ' beta_label '_betamaps_censored_condnorm_masked.nii.gz '...
      ' -a ' beta_label '_betamaps_censored_condnorm.nii.gz'...
      ' -b ' results_dir '/' mask_prefix subject_label '.nii '...
      ' -expr "a*b"']);
  
%% Load ROI labels
labels =  importdata(atlas_labels);
labels(:,2) = labels(:,1);
labels(:,1) = num2cell(1:length(labels));
num_rois = num2str(length(labels));
      
%% Extract mean beta series from ROIs, save to txt file
unix(['3dROIstats -mask ' atlas ' -numROI ' num_rois ' -1DRformat '...
       beta_label '_betamaps_censored_condnorm_masked.nii.gz > '...
       beta_label '_betamaps_censored_condnorm_masked_roivecs.txt']);
   
% %% Move data into condition directory
% unix(['mv *' beta_str_tag '* ' beta_dir_sub]);

%% Load beta series txt file and get vectors
beta_series = importdata([beta_dir_sub beta_label '_betamaps_censored_condnorm_masked_roivecs.txt']);
beta_series = beta_series.data(:,3:end);

%% Reduce down to restricted regions of interest
cols_zeros = find(all(beta_series==0)); % all zeros
beta_series(:,cols_zeros) = [];
labels(cols_zeros,:) = [];

%% Run correlations
beta_series_corr = corr(beta_series);

% Transform to Fisher z values
beta_series_corrz = atanh(beta_series_corr);

% Save Zmatrix
save([beta_dir_sub beta_label '_Zmatrix.mat'], 'beta_series_corrz');

% Calculate Z-score matrix
num_Stims_for_df = size(beta_series,1);
temp_sigma = 1/sqrt(num_Stims_for_df-3); % Standard deviation
beta_series_corrzscore = beta_series_corrz./temp_sigma;

% Save Matrices to subject directory
save([beta_dir_sub beta_label '_RMatrix.mat'], 'beta_series_corr');
save([beta_dir_sub beta_label '_ZMatrix.mat'], 'beta_series_corrz');
save([beta_dir_sub beta_label '_ZscoreMatrix.mat'], 'beta_series_corrzscore');
save([beta_dir_sub beta_label '_num_TRs_used.mat'], 'num_Stims_for_df');
save([beta_dir_sub beta_label '_stim_indices.mat'], 'beta_stim_index');

%% Plotting matrices
RdBu = cbrewer('div', 'RdBu', 256,'PCHIP');
RdBu = flip(RdBu);
% Plot Pearson correlation matrix
fig = figure('Position',[100,100,1400,1200]);
    set(gcf,'visible','off'); % Make figure not visible
    imagesc(beta_series_corr);
    colorbar; colormap(RdBu); caxis([-1 1]);
    Tick=1:size(labels,1);
    % Nicefy plot
    ax = gca;
    ax.Title.String = ['CR-' subject_label(4:6) ' Beta Series Correlation Matrix (Pearsons r) '...
                       '(num trials used = ' num2str(num_Stims_for_df) ') - ' beta_str_tag];
    ax.Title.Interpreter = 'none';  
    ax.FontSize = 10;
    ax.FontWeight = 'bold';
    ax.XTick = Tick;
    ax.YTick = Tick;
    ax.XTickLabelRotation = 271;
    ax.XTickLabel = labels(:,2);
    ax.YTickLabel = labels(:,2);
    ax.XAxis.FontSize = 5;
    ax.YAxis.FontSize = 5;

% Prepare fig for printing and save as png
set(fig,'PaperPositionMode','auto');
print('-dpng','-r300',[beta_dir_sub beta_label '_RMatrix.png']);

% Plot Z-score matrix
fig = figure('Position',[100,100,1400,1200]);
    set(gcf,'visible','off'); % Make figure not visible
    imagesc(beta_series_corrzscore);
    colorbar; colormap(RdBu); caxis([-5 10]); %[-10 25]
    Tick=1:size(labels,1);
    % Nicefy plot
    ax = gca;
    ax.Title.String = ['CR-' subject_label(4:6) ' Beta Series Correlation Matrix (Z-scores) '...
                       '(num trials used = ' num2str(num_Stims_for_df) ') - ' beta_str_tag];
    ax.Title.Interpreter = 'none';  
    ax.FontSize = 10;
    ax.FontWeight = 'bold';
    ax.XTick = Tick;
    ax.YTick = Tick;
    ax.XTickLabelRotation = 271;
    ax.XTickLabel = labels(:,2);
    ax.YTickLabel = labels(:,2);
    ax.XAxis.FontSize = 5;
    ax.YAxis.FontSize = 5;

% Prepare fig for printing and save as png
set(fig,'PaperPositionMode','auto');
print('-dpng','-r300',[beta_dir_sub beta_label '_ZscoreMatrix.png']);

% Save zscore matrix
zscore_mat = beta_series_corrzscore;

%% ****** Beta correlations BEFORE normalization *******************************************************************
% *************************************************************************
%% Mask the beta volume before extracting series from ROIs
unix(['3dcalc -prefix ' beta_label '_betamaps_censored_masked.nii.gz '...
      ' -a ' beta_label '_betamaps_censored.nii.gz'...
      ' -b ' results_dir '/' mask_prefix subject_label '.nii '...
      ' -expr "a*b"']);
  
%% Load ROI labels
labels =  importdata(atlas_labels);
labels(:,2) = labels(:,1);
labels(:,1) = num2cell(1:length(labels));
num_rois = num2str(length(labels));
      
%% Extract mean beta series from ROIs, save to txt file
unix(['3dROIstats -mask ' atlas ' -numROI ' num_rois ' -1DRformat '...
       beta_label '_betamaps_censored_masked.nii.gz > '...
       beta_label '_betamaps_censored_masked_roivecs.txt']);
   
% %% Move data into condition directory
% unix(['mv *' beta_str_tag '* ' beta_dir_sub]);

%% Load beta series txt file and get vectors
beta_series = importdata([beta_dir_sub beta_label '_betamaps_censored_masked_roivecs.txt']);
beta_series = beta_series.data(:,3:end);

%% Reduce down to restricted regions of interest
cols_zeros = find(all(beta_series==0)); % all zeros
beta_series(:,cols_zeros) = [];
labels(cols_zeros,:) = [];

%% Run correlations
beta_series_corr = corr(beta_series);

% Transform to Fisher z values
beta_series_corrz = atanh(beta_series_corr);

% Save Zmatrix
save([beta_dir_sub beta_label '_Zmatrix_prenorm.mat'], 'beta_series_corrz');

% Calculate Z-score matrix
num_Stims_for_df = size(beta_series,1);
temp_sigma = 1/sqrt(num_Stims_for_df-3); % Standard deviation
beta_series_corrzscore = beta_series_corrz./temp_sigma;

% Save Matrices to subject directory
save([beta_dir_sub beta_label '_RMatrix_prenorm.mat'], 'beta_series_corr');
save([beta_dir_sub beta_label '_ZMatrix_prenorm.mat'], 'beta_series_corrz');
save([beta_dir_sub beta_label '_ZscoreMatrix_prenorm.mat'], 'beta_series_corrzscore');
save([beta_dir_sub beta_label '_num_TRs_used_prenorm.mat'], 'num_Stims_for_df');
save([beta_dir_sub beta_label '_stim_indices_prenorm.mat'], 'beta_stim_index');

%% Plotting matrices
RdBu = cbrewer('div', 'RdBu', 256,'PCHIP');
RdBu = flip(RdBu);
% Plot Pearson correlation matrix
fig = figure('Position',[100,100,1400,1200]);
    set(gcf,'visible','off'); % Make figure not visible
    imagesc(beta_series_corr);
    colorbar; colormap(RdBu); caxis([-1 1]);
    Tick=1:size(labels,1);
    % Nicefy plot
    ax = gca;
    ax.Title.String = ['CR-' subject_label(4:6) ' Beta Series Correlation Matrix (Pearsons r) BEFORE Normalization '...
                       '(num trials used = ' num2str(num_Stims_for_df) ') - ' beta_str_tag];
    ax.Title.Interpreter = 'none';  
    ax.FontSize = 10;
    ax.FontWeight = 'bold';
    ax.XTick = Tick;
    ax.YTick = Tick;
    ax.XTickLabelRotation = 271;
    ax.XTickLabel = labels(:,2);
    ax.YTickLabel = labels(:,2);
    ax.XAxis.FontSize = 5;
    ax.YAxis.FontSize = 5;

% Prepare fig for printing and save as png
set(fig,'PaperPositionMode','auto');
print('-dpng','-r300',[beta_dir_sub beta_label '_RMatrix_prenorm.png']);

% Plot Z-score matrix
fig = figure('Position',[100,100,1400,1200]);
    set(gcf,'visible','off'); % Make figure not visible
    imagesc(beta_series_corrzscore);
    colorbar; colormap(RdBu); caxis([-5 10]); %[-10 25]
    Tick=1:size(labels,1);
    % Nicefy plot
    ax = gca;
    ax.Title.String = ['CR-' subject_label(4:6) ' Beta Series Correlation Matrix (Z-scores) BEFORE Normalization '...
                       '(num trials used = ' num2str(num_Stims_for_df) ') - ' beta_str_tag];
    ax.Title.Interpreter = 'none';  
    ax.FontSize = 10;
    ax.FontWeight = 'bold';
    ax.XTick = Tick;
    ax.YTick = Tick;
    ax.XTickLabelRotation = 271;
    ax.XTickLabel = labels(:,2);
    ax.YTickLabel = labels(:,2);
    ax.XAxis.FontSize = 5;
    ax.YAxis.FontSize = 5;

% Prepare fig for printing and save as png
set(fig,'PaperPositionMode','auto');
print('-dpng','-r300',[beta_dir_sub beta_label '_ZscoreMatrix_prenorm.png']);

% Save zscore matrix
zscore_mat_prenorm = beta_series_corrzscore;

