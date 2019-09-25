function [vals,counts] = CR7T_automask_calculate_histogram(data,mask,nbins)
%% Calculate histogram using Afni 3dhistog, send to temporary file, then load in to matlab

unix(['3dhistog -nbins ' num2str(nbins) ' -mask ' mask ' ' data ' > temphist.txt'])
d = importdata('temphist.txt');
vals = d.data(:,1);
counts = d.data(:,2);

unix('rm -f temphist.txt');