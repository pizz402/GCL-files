%% pdf4 

format compact
clear
clc

num_cells = 100; % number of cells = size(M,1) 200
n = 200; % number of genes = size(M,2) 400 

multi_weight = 1; % weights from the interval (0, multi_weight)
multi_start = 1; % initial conditions x0 from the interval (0, multi_start)
T = 40; % relaxation time

n_hold = 30; % hold genes
n_hold_zero = 0; % for ease
n_hold_one = 1; % for ease

% for DOC
threshold = 1e-3;% threshold for gene expression
overlap_type = 1;% Overlap with weights
rlowess_span = 0.1;% smooth

%% Random network - activation + inhibition + activationDim2

k_act = 2;
k_inh = 0;
p_A2 = 0;%0.01
[ A, B, A2 ] = Build_Network_activation_inhibition_activationDim2( k_act, k_inh, p_A2, n, multi_weight );

%%


%% dynamics - nodes:
%%
%% M_nhold - different x0 and holding x0(1:n_hold)
%% 1
M_nhold = M_diff_starts_function_ABA2( num_cells, A, B, A2, n_hold, multi_start, T );
%% 2
M_multiplot_ABA2( M_nhold, threshold, n_hold, rlowess_span, overlap_type)
% saveas(gcf,strcat( 'results genes/080418 summary/diff_starts_',num2str(round(1000*rand)),'.png'))
%%
%% M_diffnhold - like M_nhold but holding *different* nodes for every cell (to do later)
%%
%% M_nockout - different x0 and every cell is for holding one gene = 0 
%% 1
M_nokout = M_nockout_function_ABA2( A, B, A2, multi_start, T );
%% 2
M_multiplot_ABA2( M_nokout,threshold, n_hold, rlowess_span, overlap_type)
% saveas(gcf,strcat( 'results genes/080418 summary/nockout_',num2str(round(1000*rand)),'.png'))
%%

%% dynamics - weights: (we can add nhold to the dynamics by changing n_hold_zero)
%%
%% M_diff_weight - different x0 and noise on weights 
%% 1
sigma_weights = 0.9;% only here - normal noise on weights
M_diff_weight = M_diff_weight_function_ABA2( num_cells, A, B, A2, n_hold_zero, multi_start, T, sigma_weights);
%% 2
M_multiplot_ABA2( M_diff_weight, threshold, n_hold_zero, rlowess_span, overlap_type)
% saveas(gcf,strcat( 'results genes/080418 summary/diff_weights_SigmaWeights',strrep(num2str(sigma_weights),'.',''),'_',num2str(round(1000*rand)),'.png'))
%%
%% M_diffpart_diff_weight - like M_diff_weight but the noise is on *different part* of the weights for every cell
%% 1
fraction_noise = 1; sigma_weights = 0.4;% from here - flat noise
M_diffpart_diff_weight = M_diff_diffpartweight_function_ABA2( num_cells, A, B, A2, n_hold_zero, multi_start, T, sigma_weights, fraction_noise);
%% 2
M_multiplot_ABA2( M_diffpart_diff_weight, threshold, n_hold_zero, rlowess_span, overlap_type)
% saveas(gcf,strcat( 'results genes/080418 summary/diffpart_diff_weights_SigmaWeights',strrep(num2str(sigma_weights),'.',''),'_FractionNoise',strrep(num2str(fraction_noise),'.',''),'_',num2str(round(1000*rand)),'.png'))
%%
%% M_samepart_diff_weight - like M_diff_weight but the noise is on *one part* of the weights
%% 1 - the article choise! % with fraction_noise = 1
fraction_noise = 1; sigma_weights = 0.4;
M_samepart_diff_weight = M_diff_samepartweight_function_ABA2( num_cells, A, B, A2, n_hold_zero, multi_start, T, sigma_weights, fraction_noise);
%% 2
M_multiplot_ABA2( M_samepart_diff_weight, threshold, n_hold_zero, rlowess_span, overlap_type)
% saveas(gcf,strcat( 'results genes/080418 summary/samepart_diff_weights_SigmaWeights',strrep(num2str(sigma_weights),'.',''),'_FractionNoise',strrep(num2str(fraction_noise),'.',''),'_',num2str(round(1000*rand)),'.png'))
%%
%% completely diff weight!
%% 1
M_completely_diff_weight = M_completely_diff_weight_function_ABA2( num_cells, A, B, A2, n_hold_zero, multi_start, T);
%% 2
M_multiplot_ABA2( M_completely_diff_weight, threshold, n_hold_zero, rlowess_span, overlap_type)
%%
%%


%% diff links:
%% 1
k_act = 2;
percent_links_unfixed = 0;
%%
M_diff_links = M_diff_links_function( k_act, percent_links_unfixed, n, multi_weight, num_cells, multi_start, T);
%% 2
M_multiplot_ABA2( M_diff_links, threshold, n_hold_zero, rlowess_span, overlap_type)
%% 3 - with weights noise
sigma_weights = 0.2;
M_difflinks_diffweight = M_difflinks_diffweight_function( k_act, percent_links_unfixed, n, multi_weight, num_cells, multi_start, T, sigma_weights);
%% 4
M_multiplot_ABA2( M_difflinks_diffweight, threshold, n_hold_zero, rlowess_span, overlap_type)
%%


%% measurment noise:
%%
%% M_Null - one cell with noise 
%% 1 - the article choise! 
sigma_results = 0.5;
M_null = M_Null_function_ABA2( num_cells, A, B, A2, n_hold_zero, multi_start, T, sigma_results);
%% 2
M_multiplot_ABA2( M_null, threshold, n_hold_zero, rlowess_span, overlap_type)
% saveas(gcf,strcat( 'results genes/080418 summary/null_SigmaResults',strrep(num2str(sigma_results),'.',''),'_',num2str(round(1000*rand)),'.png'))
%%
%% M_sampl - noise because of sampling
%% 1
num_sampl = 500;
M_sample = M_sampling_ABA2( num_cells, A, B, A2, n_hold_zero, multi_start, T, num_sampl);
%% 2
M_multiplot_ABA2( M_sample, threshold, n_hold_zero, rlowess_span, overlap_type)
% saveas(gcf,strcat( 'results genes/080418 summary/sampling_NumSampl',strrep(num2str(num_sampl),'.',''),'_',num2str(round(1000*rand)),'.png'))
%%
%% M_diff_depth - noise because of different sequencing depth - not good!!!
%% 1
depth_min = 0; depth_max = 1e-2;
M_diff_depth = M_diff_sequencing_depth_ABA2( num_cells, A, B, A2, n_hold_zero, multi_start, T, depth_min, depth_max);
%% 2
M_multiplot_ABA2( M_diff_depth, threshold, n_hold_zero, rlowess_span, overlap_type)
% saveas(gcf,strcat( 'results genes/080418 summary/diff_depth_DepthMin',strrep(num2str(depth_min),'.',''),'_DepthMax',strrep(num2str(depth_max),'.',''),'_',num2str(round(1000*rand)),'.png'))
%%
%%



%% phase diagrams:
%%
%% WeightNoise + MeasureNoise
%% 1
num_itr = 10;
sigma_weights_vec = 0:0.05:1;
sigma_noise_vec = 0:0.05:1;
%% 2
tic
[CorrCellsMatrix, CorrCCCMatrix] = CorrCells_WeightNoise_MeasureNoise( num_cells, A, B, A2, n_hold_zero, multi_start, T, sigma_weights_vec, sigma_noise_vec,num_itr);
toc
%%
% randnum = round(rand()*100);
% sfile_name_Cell = strcat('CorrCellsMatrix',num2str(randnum),'numRuns',num2str(num_itr),'numCells',num2str(num_cells),'numGenes',num2str(n),'.csv');
% csvwrite(file_name_Cell,CorrCellsMatrix)
% file_name_CCC = strcat('CorrCCCMatrix',num2str(randnum),'numRuns',num2str(num_itr),'numCells',num2str(num_cells),'numGenes',num2str(n),'.csv');
% csvwrite(file_name_CCC,CorrCCCMatrix)
%% num_itr = 10; sigma_weights_vec = 0:0.05:1; sigma_noise_vec = 0:0.05:1;
CorrCCCMatrix = csvread("CorrCCCMatrix14numRuns10numCells100numGenes200.csv");
CorrCellsMatrix = csvread("CorrCellsMatrix13numRuns10numCells100numGenes200.csv");
%% 3
contour_corr_cell = 0.9;
plot_CorrCells_CorrCOC(sigma_weights_vec, sigma_noise_vec, 'Sigma Weights', 'Sigma Measure', CorrCellsMatrix, CorrCCCMatrix, contour_corr_cell)
%%


%% DiffLinksNoise + MeasureNoise
%% 1
k_act = 2;
num_itr = 1;%5
percent_links_unfixed_vec = 0.1:0.05:0.5;%0.05 
sigma_noise_vec = 0:0.05:0.8;%0.05
%% 2
tic
[CorrCellsMatrix, CorrCCCMatrix] = CorrCells_DiffLinksNoise_MeasureNoise( k_act, n, multi_weight, num_cells, multi_start, T, percent_links_unfixed_vec, sigma_noise_vec, num_itr);
toc
%%
% randnum = round(rand()*100);
% file_name_Cell = strcat('CorrCellsMatrixLinksMeasure',num2str(randnum),'numRuns',num2str(num_itr),'numCells',num2str(num_cells),'numGenes',num2str(n),'.csv');
% csvwrite(file_name_Cell,CorrCellsMatrix)
% file_name_CCC = strcat('CorrCCCMatrixLinksMeasure',num2str(randnum),'numRuns',num2str(num_itr),'numCells',num2str(num_cells),'numGenes',num2str(n),'.csv');
% csvwrite(file_name_CCC,CorrCCCMatrix)
%% k_act = 2; num_itr = 5; percent_links_unfixed_vec = 0.1:0.05:0.5; sigma_noise_vec = 0:0.05:0.8;
CorrCCCMatrix = csvread("CorrCCCMatrixLinksMeasure80numRuns5numCells100numGenes200.csv");
CorrCellsMatrix = csvread("CorrCellsMatrixLinksMeasure80numRuns5numCells100numGenes200.csv");
%% 3
contour_corr_cell = 0.5;
plot_CorrCells_CorrCOC(percent_links_unfixed_vec, sigma_noise_vec, 'Percent Links Un-Fixed', 'Sigma Measure', CorrCellsMatrix, CorrCCCMatrix, contour_corr_cell)
%%


%% DiffLinksNoise + WeightNoise
%% 1
k_act = 2;
num_itr = 5;
percent_links_unfixed_vec = 0.1:0.05:0.5;
sigma_weights_vec = 0:0.05:0.8;
%% 2
tic
[CorrCellsMatrix, CorrCCCMatrix] = CorrCells_DiffLinksNoise_WeightNoise( k_act, n, multi_weight, num_cells, multi_start, T, percent_links_unfixed_vec, sigma_weights_vec, num_itr);
toc
%%
% randnum = round(rand()*100);
% file_name_Cell = strcat('CorrCellsMatrixLinksWeight',num2str(randnum),'numRuns',num2str(num_itr),'numCells',num2str(num_cells),'numGenes',num2str(n),'.csv');
% csvwrite(file_name_Cell,CorrCellsMatrix)
% file_name_CCC = strcat('CorrCCCMatrixLinksWeight',num2str(randnum),'numRuns',num2str(num_itr),'numCells',num2str(num_cells),'numGenes',num2str(n),'.csv');
% csvwrite(file_name_CCC,CorrCCCMatrix)
%% k_act = 2; num_itr = 5; percent_links_unfixed_vec = 0.1:0.05:0.5; sigma_weights_vec = 0:0.05:0.8;
CorrCCCMatrix = csvread("CorrCCCMatrixLinksWeight26numRuns5numCells100numGenes200.csv");
CorrCellsMatrix = csvread("CorrCellsMatrixLinksWeight26numRuns5numCells100numGenes200.csv");
%% 3
contour_corr_cell = 0.5;
plot_CorrCells_CorrCOC(percent_links_unfixed_vec, sigma_weights_vec, 'Percent Links Un-Fixed', 'Sigma Weights Noise', CorrCellsMatrix, CorrCCCMatrix, contour_corr_cell)
%%


%%
font_size = 18;
set(gca,'fontsize', font_size);

%% save to pdf & png

% str name
str = strcat( '../COC TEX/COC',num2str(round(1000*rand)));

% pdf
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf, 'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print(str,'-dpdf','-fillpage');

% png
saveas(gcf,strcat( str,'.png'));

%%


%%
%% comparing dcorr and correlation matrices - IM HERE

tic

fraction_noise = 1; sigma_weights = 0.9;
sigma_results = 0.5;
num_cells_vec = 100:100:400;
num_iterations = 1;
[ stat_mean_corr_weights , stat_mean_corr_meas, stat_bcR_weights, stat_bcR_meas, stat_COC_weights, stat_COC_meas ] = mean_corr_matrix(num_cells_vec, A, B, A2, n_hold_zero, multi_start, T, sigma_weights, fraction_noise, sigma_results, num_iterations);

toc

%%

font_size = 18;

figure;
errorbar(num_cells_vec, stat_mean_corr_weights.avg, stat_mean_corr_weights.std)
hold on
errorbar(num_cells_vec, stat_mean_corr_meas.avg, stat_mean_corr_meas.std)
hold on
plot(num_cells_vec,1./sqrt(num_cells_vec))
title(strcat('sigma weights = ',num2str(sigma_weights),'sigma results = ',num2str(sigma_results)));
legend('weights','measurment','1/sqrt(numcells)')
xlabel('number of cells')
ylabel('<correlation>')
set(gca,'fontsize', font_size);

figure;
errorbar(num_cells_vec, stat_bcR_weights.avg, stat_bcR_weights.std)
hold on
errorbar(num_cells_vec, stat_bcR_meas.avg, stat_bcR_meas.std)
legend('weights','measurment')
title(strcat('sigma weights = ',num2str(sigma_weights),'sigma results = ',num2str(sigma_results)));
xlabel('number of cells')
ylabel('bcR')
set(gca,'fontsize', font_size);

figure;
errorbar(num_cells_vec, stat_COC_weights.avg, stat_COC_weights.std)
hold on
errorbar(num_cells_vec, stat_COC_meas.avg, stat_COC_meas.std)
legend('weights','measurment')
title(strcat('sigma weights = ',num2str(sigma_weights),'sigma results = ',num2str(sigma_results)));
xlabel('number of cells')
ylabel('COC')
set(gca,'fontsize', font_size);
