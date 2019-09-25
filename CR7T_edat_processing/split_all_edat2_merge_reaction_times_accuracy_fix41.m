clear all; close all;
struct_edat = tdfread('all_edat2_merge.txt');
% clearvars -except struct_edat


subjs = unique(struct_edat.Subject);
runs = unique(cellstr(struct_edat.ExperimentName));
vars = {'ExperimentName','Subject','format','numerosity','ISI'...
        'Stimuli0x2EOnsetTime',...
        'Stimuli0x2EOnsetDelay',...
        'ISIfixationSTART0x2EOnsetTime',...
        'ISIfixationSTART0x2EOnsetDelay',...
        'ISIfixationSTART0x2EOffsetTime',...
        'Stimuli0x2ERT',...
        'Stimuli0x2ERTTime',...
        'Stimuli0x2EACC',...
        'Stimuli0x2ERESP'};


for ii = 41
% Get data from rows unique to subject
    rows_subj = find([struct_edat.Subject] == ii);
    full_subj_vars = structfun(@(x) x(rows_subj,:), struct_edat, 'UniformOutput',false);
    
                rows_subj_KEEP = find([full_subj_vars.Session] == 2);
                full_subj_vars = structfun(@(x) x(rows_subj_KEEP,:), full_subj_vars, 'UniformOutput',false);
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
    if ii < 10
        save(['CR_00' num2str(stimvars.(runs{3}).Subject(1)) '_stimvars_fix+rt+acc.mat'],'-struct','stimvars');
    else
        save(['CR_0' num2str(stimvars.(runs{3}).Subject(1)) '_stimvars_fix+rt+acc.mat'],'-struct','stimvars');
    end
end