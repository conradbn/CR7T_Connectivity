%% Make module color map for plotting to surface
purge

cmap = ones(210,3);
% Set modules rgb close to Power 2011, Neuron colors
module_id_rgb = [2,1,0,0;... % DMN, red
                 6,1,1,0;... % FPN, yellow
                 10,0,1,1;... % Sensorimotore, cyan
                 12,0,0,0;... % Salience, black
                 16,0.3,0.3,0.3;... % MTG, dark gray
                 17,1,0.5,1;... % Auditory, pink
                 20,0,1,0;... % ITG??, green
                 21,0,0,1;... % Visual, blue
                 22,0,0.5,0.5;... % Hippocampus, teal
                 24,1,0.5,0;... % Ventral Attention (ITG/MTG), orange
                 27,0.5,0.5,1;... % Precuneus, pale blue
                 34,0.3,0,0;... % Basal ganglia, dark brown
                 35,1,0,1;... % Caudate, purple
                 36,0.6,0,0]; % Thalumus, light brown
             
for ii = 1:length(module_id_rgb)
    ind = module_id_rgb(ii);
    cmap(ind,:) = module_id_rgb(ii,2:end);
end