function output = CR7T_get_regressor_indices(subject_dir)

%% Load experiment variables
cd(subject_dir);
stimvars_file = dir('*_stimvars+rt+acc.mat');
stimvars = load(stimvars_file.name);
runs = {'identify_A','identify_B','compare_A','compare_B'};

task = [];
format = [];
numerosity = [];
acc = [];
rt = [];

for ii = 1:size(runs,2)
    if ii < 3
        task = [task ones(1,size(stimvars.(runs{ii}).Stimuli0x2EOnsetTime,1))];
    else
        task = [task 2*ones(1,size(stimvars.(runs{ii}).Stimuli0x2EOnsetTime,1))];
    end
    format = [format stimvars.(runs{ii}).format'];
    numerosity = [numerosity stimvars.(runs{ii}).numerosity'];
    acc = [acc stimvars.(runs{ii}).Stimuli0x2EACC'];
    rt = [rt stimvars.(runs{ii}).Stimuli0x2ERT'];
end

%num_timepoints = size(task,2);

regressors.task(1,:) = double(task == 1);
regressors.task(2,:) = double(task == 2);

regressors.format(1,:) = double(format == 1);
regressors.format(2,:) = double(format == 2);

regressors.numerosity(1,:) = double(numerosity == 2);
regressors.numerosity(2,:) = double(numerosity == 4);
regressors.numerosity(3,:) = double(numerosity == 6);
regressors.numerosity(4,:) = double(numerosity == 8);

regressors.acc = acc;
regressors.rt = rt; 

output = regressors;