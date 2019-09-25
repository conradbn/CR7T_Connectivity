function CR7T_segment_anat_final(results_dir)
%% CR7T_segment_anat_final
% This function runs a custom segmentation pipeline on the normalized 
% anatomical volumes, to create tissue masks for anatomical principal
% component correction

cd(results_dir)
unix('rm -rf Seg'); % Remove directory if exists
unix('mkdir Seg');
anat = dir('anat_final.*+tlrc.HEAD');
anat = anat(1).name;
master = dir('mask_epi_anat.*+tlrc.HEAD');
master = master(1).name;
%% Remove hyperintense voxels
unix(['3dcalc -prefix Seg/anat_final_rm_hyperintensity.nii'...
      ' -a ' anat...
      ' -expr "a * isnegative(a - 0.45)"']);
      
%% Run AFNI segmentation (no bias correction)

% ** This (including extra class called "Hyper") didnt work for subjects 
% ** which had relatively less hyperintense areas...
% unix(['3dSeg -prefix Seg -anat ' anat ' -mask ' anat...
%       ' -classes "CSF; GM ; WM; Hyper;"'...
%       ' -bias_fwhm 0.0 -mixfrac UNI -Bmrf 1.0 -main_N 5 -blur_meth BFT'],'-echo');

disp('Segmenting....');

unix(['3dSeg -prefix Seg -anat Seg/anat_final_rm_hyperintensity.nii'...
      ' -mask Seg/anat_final_rm_hyperintensity.nii'...
      ' -classes "CSF; GM ; WM;"'...
      ' -bias_fwhm 0.0 -Bmrf 1.0 -main_N 5 -blur_meth BFT'],'-echo');
  
unix('3dcalc -a Seg/Classes+tlrc"<WM>" -expr "step(a)" -prefix Seg/mask_WM.nii');
unix('3dcalc -a Seg/Classes+tlrc"<CSF>" -expr "step(a)" -prefix Seg/mask_CSF.nii');

%% WM 
% Generate eroded/resampled masks for WM
unix(['3dmask_tool -input Seg/Classes+tlrc"<WM>" -dilate_input -1'...
      ' -prefix Seg/mask_WMe.nii']);
unix(['3dresample -master ' master ' -rmode NN'...
      ' -input Seg/mask_WMe.nii -prefix Seg/mask_WMe_resam.nii']);

% Eliminate potential cerebellum voxels (due to hyperintensities in
% inferior regions) + apply epi mask
cmask = '/mnt/CR7T_Connectivity/matfiles/ventricle_cerebellum_masks/cerebellum_exterior_Neuromorphometrics_dilate_resam.nii';
unix(['3dcalc -prefix Seg/mask_WMe_resam_final.nii'...
      ' -a Seg/mask_WMe_resam.nii'...
      ' -b ' cmask...
      ' -c ' master...
      ' -expr "a * not(b) * c"']);
      
%% Ventricles      
% Generate resampled masks for CSF (no erosion) 
unix(['3dresample -master ' master ' -rmode NN'...
      ' -input Seg/mask_CSF.nii -prefix Seg/mask_CSF_resam.nii']);
  
% Restrict to lateral ventricles + apply epi mask
vmask = '/mnt/CR7T_Connectivity/matfiles/ventricle_cerebellum_masks/lat_ventricles_Neuromorphometrics_dilate_resam.nii';
unix(['3dcalc -prefix Seg/mask_Vent_resam.nii'...
      ' -a Seg/mask_CSF_resam.nii'...
      ' -b ' vmask...
      ' -c ' master...
      ' -expr "a * b * c"']);