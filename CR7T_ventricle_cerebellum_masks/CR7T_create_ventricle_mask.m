purge; 

nii_orig = load_untouch_nii('labels_Neuromorphometrics.nii');
new_nii = nii_orig;
new_nii.img = single(new_nii.img);
new_nii.img(new_nii.img == 51) = -1;
new_nii.img(new_nii.img == 52) = -1;
new_nii.img(new_nii.img >= 0) = 0;
new_nii.img(new_nii.img == -1) = 1;

save_untouch_nii(new_nii,'lat_ventricles_Neuromorphometrics.nii');


unix(['3dmask_tool -input lat_ventricles_Neuromorphometrics.nii -dilate_input 1'...
      ' -prefix lat_ventricles_Neuromorphometrics_dilate.nii']);
unix(['3dresample -master mask_epi_anat.CR_001_313912+tlrc.BRIK -rmode NN'...
      ' -input lat_ventricles_Neuromorphometrics_dilate.nii -prefix lat_ventricles_Neuromorphometrics_dilate_resam.nii']);
