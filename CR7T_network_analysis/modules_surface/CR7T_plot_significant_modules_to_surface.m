%% Plot significant modules to surface
purge
%% Base colormap
cmap_base = ones(256,3);
% % Set modules rgb close to Power 2011, Neuron colors
% module_id_rgb = [2,1,0,0;... % DMN, red
%                  6,1,1,0;... % FPN, yellow
%                  10,0,1,1;... % Sensorimotor, cyan
%                  12,0.8,0,1;... % Cingulo/Opercular/Salience, purple
%                  16,0.3,0.3,0.3;... % STS, dark gray
%                  17,1,0.7,1;... % Auditory, light pink
%                  20,1,0.5,0;... % ventral Caudate, orange
%                  21,0,0,1;... % Visual, blue
%                  22,0,0.5,0.5;... % Hippocampus, teal
%                  24,0,1,0;... % Dorsal Attention (ITG/MTG),green
%                  27,0.5,0.5,1;... % Precuneus, pale blue
%                  34,0.3,0,0;... % Basal ganglia, dark brown
%                  35,0,0,0;... % Caudate, black
%                  36,0.6,0,0;...% Thalumus, light brown
%                  256,0.1,0.1,0.1]; % Signal dropout regions
% for ii = 1:length(module_id_rgb)
%     ind = module_id_rgb(ii);
%     cmap_base(ind,:) = module_id_rgb(ii,2:end);
% end

module_id_rgb = {'DMN',2,1,0,0;... % DMN, red
                 'FPN',6,1,1,0;... % FPN, yellow
                 'SMN',10,0,1,1;... % Sensorimotor, cyan
                 'Sal',12,0.8,0,1;... % Cingulo/Opercular/Salience, purple
                 'STS',16,1,0.7,1;... % STS, dark gray
                 'CaudV',20,1,0.5,0;... % ventral Caudate, orange
                 'Vis',21,0,0,1;... % Visual, blue
                 'Hipp',22,0,0.5,0.5;... % Hippocampus, teal
                 'DAN',24,0,1,0;... % Dorsal Attention (ITG/MTG),green
                 'Prec',27,0.5,0.5,1;... % Precuneus, pale blue
                 'BG',34,0.3,0.1,0;... % Basal ganglia, dark brown
                 'CaudD',35,.8,.33,0;... % Caudate, burnt orange
                 'Thal',36,0.5,0.25,0.1;...% Thalumus, light brown
                 'Dropout',256,0.1,0.1,0.1}; % Signal dropout regions
             
for ii = 1:length(module_id_rgb)
    ind = module_id_rgb{ii,2};
    cmap_base(ind,:) = cell2mat(module_id_rgb(ii,3:end));
end

%% Symbolic ---------------------------------------------------------------------
load('workspace_Qc_significance_symbolic.mat','Qc_sig_zval','Qc_sig_pval','modules_true','C');
mods = modules_true;
% mz_all = mean(Qc_sig_zval,2);
% mz_sig = mz_all > 2.33;
mz_all = median(Qc_sig_pval,2);
mz_sig = mz_all < 0.01;
mz_ind = mods(mz_sig);
C_new = C .* ismember(C,mz_ind);

    % Read in data (vector of values corresponding to each ROI in the label/surface files)
    data = CR7T_convert_vec_BN202_to_BN246(C_new);
    % Add in values at 256 for signal dropout regions
    d = CR7T_convert_vec_BN202_to_BN246(C);
    data(d == 0) = 256;
    % Color range (min val, max val)
    colorrange = [1,256];
    % Surface gifti file
    fs_dir = 'fsaverage_BrainnetomeAtlas\164k\' ;
    surface = [fs_dir 'fsaverage.R.very_inflated.164k.surf.gii'];
    % Label gifti file
    label = [fs_dir 'fsaverage.R.BN_Atlas.164k.label.gii'];
    % Output directory
    out_dir = 'Z:\CR7T_Connectivity\matfiles\CR7T_network_analysis\modules_surface';
    % Output file name
    out_names = {'Symbolic_lat_r_90_0';...
                'Symbolic_med_r_270_0';...
                'Symbolic_sup_r_0_90'};           
    % Viewing angle
    view_angs = [90,0;...
                270,0;...
                0,90];
    % ----------- RUN FUNCTION --------------
    for ii = 1:numel(out_names)
        out_name = out_names{ii};
        view_ang = view_angs(ii,:);
        plot_on_indexed_surface(data,cmap_base,colorrange,surface,label,out_dir,out_name,view_ang);
    end
    
    % Surface gifti file
    fs_dir = 'fsaverage_BrainnetomeAtlas\164k\' ;
    surface = [fs_dir 'fsaverage.L.very_inflated.164k.surf.gii'];
    % Label gifti file
    label = [fs_dir 'fsaverage.L.BN_Atlas.164k.label.gii'];
    % Output directory
    out_dir = 'Z:\CR7T_Connectivity\matfiles\CR7T_network_analysis\modules_surface';
    % Output file name
    out_names = {'Symbolic_lat_l_270_0';...
                'Symbolic_med_l_90_0';...
                'Symbolic_sup_l_0_90'};            
    % Viewing angle
    view_angs = [270,0;...
                90,0;...
                0,90];
    % ----------- RUN FUNCTION --------------
    for ii = 1:numel(out_names)
        out_name = out_names{ii};
        view_ang = view_angs(ii,:);
        plot_on_indexed_surface(data,cmap_base,colorrange,surface,label,out_dir,out_name,view_ang);
    end
    
%% Nonsymbolic ---------------------------------------------------------------------
load('workspace_Qc_significance_nonsymbolic.mat','Qc_sig_zval','Qc_sig_pval','modules_true','C');
mods = modules_true;
% mz_all = mean(Qc_sig_zval,2);
% mz_sig = mz_all > 2.33;
mz_all = median(Qc_sig_pval,2);
mz_sig = mz_all < 0.01;
mz_ind = mods(mz_sig);
C_new = C .* ismember(C,mz_ind);
    % Read in data (vector of values corresponding to each ROI in the label/surface files)
    data = CR7T_convert_vec_BN202_to_BN246(C_new);    
    % Add in values at 256 for signal dropout regions
    d = CR7T_convert_vec_BN202_to_BN246(C);
    data(d == 0) = 256;
    % Color range (min val, max val)
    colorrange = [1,256];
    % Surface gifti file
    fs_dir = 'fsaverage_BrainnetomeAtlas\164k\' ;
    surface = [fs_dir 'fsaverage.R.very_inflated.164k.surf.gii'];
    % Label gifti file
    label = [fs_dir 'fsaverage.R.BN_Atlas.164k.label.gii'];
    % Output directory
    out_dir = 'Z:\CR7T_Connectivity\matfiles\CR7T_network_analysis\modules_surface';
    % Output file name
    out_names = {'Nonsymbolic_lat_r_90_0';...
                'Nonsymbolic_med_r_270_0';...
                'Nonsymbolic_sup_r_0_90'};            
    % Viewing angle
    view_angs = [90,0;...
                270,0;...
                0,90];
    % ----------- RUN FUNCTION --------------
    for ii = 1:numel(out_names)
        out_name = out_names{ii};
        view_ang = view_angs(ii,:);
        plot_on_indexed_surface(data,cmap_base,colorrange,surface,label,out_dir,out_name,view_ang);
    end
    
    % Surface gifti file
    fs_dir = 'fsaverage_BrainnetomeAtlas\164k\' ;
    surface = [fs_dir 'fsaverage.L.very_inflated.164k.surf.gii'];
    % Label gifti file
    label = [fs_dir 'fsaverage.L.BN_Atlas.164k.label.gii'];
    % Output directory
    out_dir = 'Z:\CR7T_Connectivity\matfiles\CR7T_network_analysis\modules_surface';
    % Output file name
    out_names = {'Nonsymbolic_lat_l_270_0';...
                'Nonsymbolic_med_l_90_0';...
                'Nonsymbolic_sup_l_0_90'};           
    % Viewing angle
    view_angs = [270,0;...
                90,0;...
                0,90];
    % ----------- RUN FUNCTION --------------
    for ii = 1:numel(out_names)
        out_name = out_names{ii};
        view_ang = view_angs(ii,:);
        plot_on_indexed_surface(data,cmap_base,colorrange,surface,label,out_dir,out_name,view_ang);
    end




