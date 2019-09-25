%% Catenate QC images from all subjects 
purge
process_dir = '/mnt/CR7T_Connectivity/MRI_proc';

% Set file names/output
% image_prefix = 'IMG_02_vol_check_stats_anat.axi.jpg';
% output_name = '/mnt/CR7T_Connectivity/figures/qc_all_subs_Fstat_axial.png';
% image_prefix = 'IMG_02_vol_check_stats_anat.sag.jpg';
% output_name = '/mnt/CR7T_Connectivity/figures/qc_all_subs_Fstat_sagittal.png';
% image_prefix = 'IMG_02_vol_check_stats_anat.cor.jpg';
% output_name = '/mnt/CR7T_Connectivity/figures/qc_all_subs_Fstat_coronal.png';
% image_prefix = 'IMG_05_1D_enorm_mot.jpg';
% output_name = '/mnt/CR7T_Connectivity/figures/qc_all_subs_enorm_motion.png';

%% Vertical catenate
% image_prefix = 'IMG_04_1D_cen_out.jpg';
% output_name = 'Z:/CR7T_Connectivity/figures/qc_all_subs_censored_outlier.png';
% % Go to process directory and get all subject file paths
% cd(process_dir)
% fnames = subdir(image_prefix);

% % Make string of all file paths/names, separated by a space
% fstring = [];
% for ii = 1:numel(fnames)
%     fstring = [fstring ' ' fnames(ii).name];
% end
% % Run catenate command
% unix(['convert ' fstring ' -append ' output_name]);


%% Square catenate
% image_prefix = 'IMG_04_1D_cen_out.jpg';
% output_name = '/mnt/CR7T_Connectivity/figures/qc_all_subs_censored_outlier_square.png';
image_prefix = 'IMG_05_1D_enorm_mot.jpg';
output_name = '/mnt/CR7T_Connectivity/figures/qc_all_subs_enorm_motion_square.png';
image_prefix = 'IMG_02_vol_check_stats_anat.sag.jpg';
output_name = '/mnt/CR7T_Connectivity/figures/qc_all_subs_Fstat_sagittal_square.png';
image_prefix = 'IMG_02_vol_check_stats_anat.cor.jpg';
output_name = '/mnt/CR7T_Connectivity/figures/qc_all_subs_Fstat_coronal_square.png';
% image_prefix = 'IMG_05_1D_enorm_mot.jpg';
% output_name = '/mnt/CR7T_Connectivity/figures/qc_all_subs_enorm_motion_square.png';
cd(process_dir)
fnames = subdir(image_prefix);
% Get as close as possible to square layout
num_tot = size(fnames,1);
num_rows = round(sqrt(num_tot));
num_cols = round(num_tot/num_rows);
% If more rows than columns, reverse it 
if num_rows > num_cols
    n = num_rows;
    num_rows = num_cols;
    num_cols = n;
end

% Make string of all file paths/names, separated by a space
ind = 0;
for ii = 1:num_rows
    fstring = [];
    for jj = 1:num_cols
        if ind+jj < num_tot
            fstring = [fstring ' ' fnames(ind+jj).name];
        end
    end
    ind = ind + num_cols;
    % Run catenate command
    unix(['convert ' fstring ' +append rm' num2str(ii) '.png']);
end

% Run catenate command
unix(['convert rm*.png -append ' output_name]);
% Remove temp files
unix('rm -f rm*.png');


    

