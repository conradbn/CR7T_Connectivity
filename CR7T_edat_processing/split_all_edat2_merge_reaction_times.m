%% Script to get data from tab delimited, merged edat file for Price 7T project
% Merged file was created in E-merge Eprime program, contains .edat2 data
% for all subjects/runs

clear all; close all;
struct_edat = tdfread('all_edat2_merge.txt');
% clearvars -except struct_edat


subjs = unique(struct_edat.Subject);
runs = unique(cellstr(struct_edat.ExperimentName));
vars = {'ExperimentName','Subject','format','numerosity','ISI'...
        'Stimuli0x2EOnsetTime',...
        'Stimuli0x2EOnsetDelay',...
        'ISIfixationSTART0x2EOnsetTime',...
        'ISIfixationSTART0x2EOnsetDelay'...
        'ISIfixationSTART0x2EOffsetTime'...
        'Stimuli0x2ERT'...
        'Stimuli0x2ERTTime'};

% Loop through each subject
for ii = 1:numel(subjs)
    % Get data from rows unique to subject
    rows_subj = find([struct_edat.Subject] == ii);
    full_subj_vars = structfun(@(x) x(rows_subj,:), struct_edat, 'UniformOutput',false);
    % Loop through experimental runs
    for jj = 1:numel(runs)
        % Get data from rows unique to run
        rows_exp = find(strcmp(cellstr(full_subj_vars.ExperimentName),runs{jj}));
        subj_exp_vars_full.(runs{jj}) = structfun(@(x) x(rows_exp,:), full_subj_vars, 'UniformOutput',false);
        % Loop through variables of interest and create new structure
        for kk = 1:numel(vars)
            stimvars.(runs{jj}).(vars{kk}) = subj_exp_vars_full.(runs{jj}).(vars{kk});          
        end
    end
    if ii == 39 % Skip subject 39, no data
        continue 
    end
    if ii < 10
        save(['CR_00' num2str(stimvars.(runs{3}).Subject(1)) '_stimvars+rt.mat'],'-struct','stimvars');
    else
        save(['CR_0' num2str(stimvars.(runs{3}).Subject(1)) '_stimvars+rt.mat'],'-struct','stimvars');
    end
end

% I am calculating the onset time the following way, which accounts for the
% first stimulus onset delay
x = stimvars.compare_A.Stimuli0x2EOnsetTime - stimvars.compare_A.ISIfixationSTART0x2EOnsetTime;

% Eric's method, assumes 2000ms since the beginning of acquisition (i.e.
% does not account for stimulus onset delay which extends the time to first
% onset +5-40ms. Furthermore neither of these account for initial ISI onset delay... 
y = stimvars.compare_A.Stimuli0x2EOnsetTime(1) - 2000;
z = stimvars.compare_A.Stimuli0x2EOnsetTime - y;

