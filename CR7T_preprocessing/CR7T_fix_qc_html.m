%% Fix QC html output for the CR7T project
% The corrected MP2RAGE volumes have a value range of ~0-1, which causes
% 3dedge3 (called within the QC html creation routine) to fail and produce no edges
% for the anatomical-template overlay. This script makes a small
% modification, calling a custom file @djunct_edgy_align_check_scale, which
% includes a flag to scale the voxel values before edge detection

cd(results_dir);

% Remove original QC folder
unix(['rm -rf QC_' subject_label]);

% Replace string
fid  = fopen('@ss_review_html','r');
f=fread(fid,'*char')';
fclose(fid);
f = strrep(f,'@djunct_edgy_align_check','@djunct_edgy_align_check_scale');
fid  = fopen('@ss_review_html_3dedge3_scale','w');
fprintf(fid,'%s',f);
fclose(fid);

% Rerun the process
unix('./@ss_review_html_3dedge3_scale');
unix(['apqc_make_html.py -qc_dir QC_' subject_label]);
