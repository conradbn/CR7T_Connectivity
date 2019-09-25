%% CR7T Batch generation of BrainNet Viewer Images
purge
cd('Z:\CR7T_Connectivity\matfiles\CR7T_network_analysis\condition_comparison\BrainNetViewer_Plotting')

mesh_name = 'BrainMesh_ICBM152_smoothed.nv';
edge_names = dir('Pa_Diff_Sym-Nonsym_Gamma_2.45.*combine_0.12.edge');
node_names = dir('BN202_*0.12.node');
settings_names = {'brainnetview_node_profile_settings_Left.mat';...
                 'brainnetview_node_profile_settings_Right.mat';...
                 'brainnetview_node_profile_settings_Back.mat';...
                 'brainnetview_node_profile_settings_Top.mat'};
dir_label = {'_l';'_r';'_b';'_t'};

for ii = 1:numel(edge_names)
    % Set file names
    edge_name = edge_names(ii).name;
    node_name = node_names(ii).name;
    for jj = 1:numel(dir_label)
        % Set direction/settings
        settings_name = settings_names{jj};
        png_name = [node_name dir_label{jj} '.png'];
        % Run the function
        BrainNet_MapCfg(mesh_name,edge_name,node_name,settings_name,png_name);
    end
end