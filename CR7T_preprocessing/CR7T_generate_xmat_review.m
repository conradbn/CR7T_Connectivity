function CR7T_generate_xmat_review(subject_label)
%% This script creates the output review info necessary for beta-series scrubbing
%% Generate and run review scripts, in order to get censored TR's per stim a la afni_proc.py method
%mkdir([results_dir 'beta_series_review']);


% Generate scripts to review single subject results
unix(['gen_ss_review_scripts.py -prefix  ' prefix '.'...
      ' -xmat_regress ' subject_label '_beta_series.X.xmat.1D'... 
      ' -xmat_uncensored ' subject_label '_beta_series.uncensored.X.xmat.1D'...
      ' -mot_limit 0.3'...
      ' -out_limit 0.05']);
  
%% Run basic subject review script
unix(['./' prefix '.@ss_review_basic |& tee ' prefix '.out.ss_review.txt']);