function [out,vars] = CR7T_extract_QA_vals(subj_proc_folder)
%% CR7T_extract_QA_vals
% This function uses afni's gen_ss_review_table.py function to extract the
% values from out.ss_review text files. Also get the TSNR from the volreg file,
% which isn't included in out.ss_review process.

%% Setup and run gen_ss_review_tabl.py
infiles = ['/mnt/CR7T_Connectivity/MRI_proc/' subj_proc_folder '/*.results/out.ss_review.*.txt'];
out_name = ['/mnt/CR7T_Connectivity/quality_assurance/review_stats_' subj_proc_folder '.xls'];


unix(['gen_ss_review_table.py -overwrite -verb 2 -separator :'...
      ' -tablefile ' out_name ...
      ' -infiles ' infiles]);

% Read the file
fileID = fopen(out_name);
T  = textscan(fileID, '%s','delimiter', '\t');
T = T{1,1};
% Get rows containing value (middle row)
inds = find(contains(T,'value'));
% Set up reduced tables
T1 = T(1:inds(1)-1); 
T1{end+1} = ''; % Have to add a cell, due to bug in gen_ss_review
T2 = T(inds);
T3 = T(inds(end)+1:numel(T));

vars = {'subject ID',1;...
        'average motion (per TR)',1;...
        'max motion displacement',1;...
        'average outlier frac (TR)',1;...
        'fraction censored per run',4;...
        'censor fraction',1;...
        'fraction TRs censored',4;...
        'TSNR average',1;...
        'global correlation (GCOR)',1;...
        'anat/EPI mask Dice coef',1;...
        'maximum F-stat (masked)',1;...
        'AFNI version',1};

%% Loop through vars and extract vals
ind = 1;
for ii = 1:size(vars,1)
    row = find(contains(T1,vars{ii,1}));
    % Skip if doesn't exist
    if isempty(row)
        out{1,ind} = 'NaN';
        ind = ind + 1;
        continue
    end
    % Check if only one value is expected for this variable
    if vars{ii,2} == 1
        out{1,ind} = T3{row};
        ind = ind + 1;
    else
    % Since more than one value expected, loop through and get each value of this variable
        for jj = 1:vars{ii,2}
            if ~strcmp(T2{row-1+jj},['value_' num2str(jj)])
                disp(['** WARNING - While parsing file for ' subj '... could not find value #' num2str(jj) ' for *' vars{ii,1} '*... check manually!']);
                out{1,ind} = 'NaN';
            else
                out{1,ind} = T3{row-1+jj};
            end
            ind = ind + 1;
        end
    end
end
% 
% %% Get TSNR from volreg file
% tsnr = dir(['/mnt/CR7T_Connectivity/MRI_proc/' subj_proc_folder '/*.results/TSNR.vreg.mean.txt']);
% tsnr = load([tsnr.folder '/' tsnr.name]);
% out{1,ind} = tsnr;




