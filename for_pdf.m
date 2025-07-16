clear
clc

num_cells = 100; % number of cells = size(M,1) 100
n = 200; % number of genes = size(M,2) 200 

multi_weight = 1; % weights from the interval (0, multi_weight)
multi_start = 1; % initial conditions x0 from the interval (0, multi_start)
T = 10; % relaxation time

n_hold = 30; % hold genes
n_hold_zero = 0;

overlap_type = 1;% Overlap
rlowess_span = 0.1;% for DOC

threshold_vec = 1e-3;% for DOC

%%

k_act = 2;
k_inh = 0;
[ A, B ] = Build_Network_activation_inhibition( k_act, k_inh, n, multi_weight );
% [ A, B ] = Build_SF_Network_activation_inhibition( k_act, k_inh, n);
p_2inter = 0.005;%0.005
C = cat(3, full(multi_weight*sprand(n,n,p_2inter)), full(multi_weight*sprand(n,n,p_2inter)));
C(1:n+1:(end/2)) = 0 ; C(1+(end/2):n+1:end) = 0 ;

%% M_Null - one cell with noise 

sigma_results = 1.5; % the noise

M_Null = M_Null_function( num_cells, A, B, n_hold, multi_start, T, sigma_results);

% figure;
% histogram(M_Null(1,:),50,'Normalization','probability')
% title('Example of a cell')
% saveas(gcf,'Example of cell.png')

%% M_diff_starts - different x0 (including different x0(1:n_hold))

M_diff_starts = M_diff_starts_function( num_cells, A, B, n_hold, multi_start, T);

%% 
new_M_diff_starts = new_M_diff_starts_function( num_cells, A, B, n_hold, multi_start, T, C);

%% M_diff_weight - different x0 and weights with niose 

sigma_weights = 0.25; % the noise

M_diff_weight = M_diff_weight_function( num_cells, A, B, n_hold, multi_start, T, sigma_weights);

%% M_diff_partweight - different x0 and weights with noise for the same part of the weights

fraction_noise = 0.2;
sigma_weights = 0.5; % the noise

M_diff_partweight = M_diff_partweight_function( num_cells, A, B, n_hold, multi_start, T, sigma_weights,fraction_noise);

%% M_diff_partweight - different x0 and weights with noise for different parts of the weights

fraction_noise = 0.2;
sigma_weights = 0.5; % the noise

M_diff_diffpartweight = M_diff_diffpartweight_function( num_cells, A, B, n_hold, multi_start, T, sigma_weights,fraction_noise);

%% M_2noises -  different x0 and weights with niose and measurement noise

% sigma_results = 0;
% sigma_weights = 0.5;
% 
% M_2noises = M_2noises_function( num_cells, A, B, n_hold, multi_start, T, sigma_results, sigma_weights);

%% New Dynamics - canceling percent of the nodes

percent_canceled = 0.1;
sigma_results = 0;

M_New = NewDynamics_function( num_cells, A, B, percent_canceled, multi_start, T, sigma_results);

%% M_sampl - Noise only because of sampling (with the same sequencing depth)

num_sampl = 100; depth_threshold = 0; 

M_sampl = same_sequencing_depth( num_cells, A, B, n_hold, multi_start, T, num_sampl , depth_threshold);

%% M_diff_depth - Noise only because different sequencing depth

depth_min = 0; depth_max = 5e-3;

M_diff_depth = diff_sequencing_depth( num_cells, A, B, n_hold, multi_start, T, depth_min, depth_max);

%% M_bionoise - Molecular noise

c_noise = 1000; sigma = 0.1;

tic
M_bionoise = bionoise_M_diff_starts_function( num_cells, A, B, multi_start, T, c_noise, sigma );
toc

%% tmp
figure;
histogram(M_bionoise(1,:),50,'Normalization','probability')
title('Example of a cell')
saveas(gcf,'Example of cell.png')

%% tmp

threshold = 1e-3;

DOC(M_bionoise, threshold, rlowess_span, overlap_type );

%% %%%%%%%

%% %%%%%%%

%% PCA 

%% PCA for (any) 3 options

plot_scatter_score_pca( M_Null, M_diff_starts, M_diff_weight )
legend('Null', 'different x0','different x0 and weights');
% saveas(gcf,'PCA.png')

%% %%%%%%%

%% Correlations

%% Spearman Correlation between cells

figure;
plot_histogram_corr_Spearman( M_Null )
plot_histogram_corr_Spearman( M_diff_starts )
plot_histogram_corr_Spearman( M_diff_weight )
plot_histogram_corr_Spearman( M_New )
% plot_histogram_corr_Spearman( M_2noises )
legend('Null', 'different x0','different x0 and weights','new')%,'2 noises');
title(strcat('Spearman Correlation - ','cells=',num2str(num_cells),', genes=',num2str(n), ', hold =', num2str(n_hold),', act =',num2str(k_act),', inh =',num2str(k_inh) ));
hold off
% saveas(gcf,'Spearman Correlation.png')

%% between genes
figure;
[rho_null, pval_null] = corr(M_Null(:,n_hold+1:end),'Type','Spearman');
imagesc(rho_null)
figure;
[rho_starts, pval_starts] = corr(M_diff_starts(:,n_hold+1:end),'Type','Spearman');
imagesc(rho_starts)

%% corr on unhold genes
figure;
plot_histogram_corr_Spearman( M_Null(:,n_hold+1:end)' )
plot_histogram_corr_Spearman( M_diff_starts(:,n_hold+1:end)' )
legend('Null', 'different x0');

%% Spearman for noise from different sampling (and no depth)

depth_threshold = 0;
vec_num_sampl = 200:200:400;

plot_Spearman_sampling_diffsampl(vec_num_sampl,num_cells, A, B, n_hold, multi_start, T , depth_threshold)
title(strcat('Spearman Correlation - ','cells=',num2str(num_cells),', genes=',num2str(n), ', hold =', num2str(n_hold),', act =',num2str(k_act),', inh =',num2str(k_inh),', depth = ',num2str(depth_threshold) ));

%% %%%%%%%

%% DOC

%% DOC for all 4 options ('Null','different x0','different x0 and weights','New (canceling)')

rlowess_span = 0.2;
threshold_vec = 1e-3;% for M_New the threshold stays 0

plot_DOC_3M( rlowess_span, threshold_vec, M_Null, M_diff_starts, M_diff_weight , M_New, overlap_type)
title(strcat('DOC with threshold = ',num2str(threshold_vec),', cells=',num2str(num_cells),', genes=',num2str(n), ', hold =', num2str(n_hold),', act =',num2str(k_act),', inh =',num2str(k_inh) ));

%% DOC for different n_hold - only for M_different_x0

rlowess_span = 0.2; threshold_vec = 1e-3;

n_hold_vector = [10,25,50];

plot_DOC_MdiffStarts_diffNhold( n_hold_vector, num_cells, A, B, multi_start, T, rlowess_span, threshold_vec, overlap_type)
title(strcat('DOC with threshold = ',num2str(threshold_vec),', cells=',num2str(num_cells),', genes=',num2str(n), ', hold =', num2str(n_hold),', act =',num2str(k_act),', inh =',num2str(k_inh) ));

%% analyzing DOC for noise in results vs noise in dynamics (long run)

tic

threshold = 1e-3;
sigma_results_vector = linspace(0.4,0.7,20);
num_runs = 20;

noiseR_MhSC( num_cells, A, B, multi_start, T, threshold, n_hold, sigma_results_vector, num_runs , overlap_type)

toc

tic

threshold = 1e-3;
sigma_weights_vector = linspace(0.1,0.8,20);
num_runs = 20;

noiseW_MhSC( num_cells, A, B, multi_start, T, threshold, n_hold, sigma_weights_vector, num_runs , overlap_type)

toc

%%

MhSC = get_MhSC_ff( "noiseR_MhSC24numRuns20numCells400numGenes600nHold50.csv" );

figure;
errorbar(MhSC.mean_M_mean_hist ,MhSC.mean_M_DOC_slope,MhSC.std_M_DOC_slope)
hold on
errorbar(MhSC.mean_M_mean_hist,MhSC.mean_M_DOC_corr,MhSC.std_M_DOC_corr)
% legend('Measurement errors slope','Measurement errors corr')
xlabel('Spearman')

%%

MhSC = get_MhSC_ff( "noiseW_MhSC11numRuns20numCells400numGenes600nHold50.csv" );

% figure;
errorbar(MhSC.mean_M_mean_hist ,MhSC.mean_M_DOC_slope,MhSC.std_M_DOC_slope)
hold on
errorbar(MhSC.mean_M_mean_hist,MhSC.mean_M_DOC_corr,MhSC.std_M_DOC_corr)
% legend('Noise in dynamics slope','Noise in dynamics corr')
legend('Measurement errors slope','Measurement errors corr','Noise in dynamics slope','Noise in dynamics corr')
xlabel('Spearman')

%% DOC for M2noise

tic

threshold = 1e-3;

sigma_results_vec = 0:0.1:0.2;
sigma_weights_vec = 0:0.1:0.5;

[slope_matrix, corr_matrix ] = plot_DOC_diff_M2noises( num_cells, A, B, n_hold, multi_start, T, sigma_results_vec, sigma_weights_vec, threshold, overlap_type);

toc

%% DOC for New Dynamics (canceling insted of zeros)

tic

sigma_results_vec = 0:0.1:0.2;
percent_canceled_vec = 0.05:0.05:0.2;

[slope_matrix, corr_matrix ] = plot_DOC_diff_MNew( num_cells, A, B, multi_start, T, sigma_results_vec, percent_canceled_vec, overlap_type);

toc

%% DOC for diff x0 for diff k_inh & k_act:

n_hold = 0.1*n;
threshold = 1e-3;

k_act_vec = 0:0.5:1;
k_inh_vec = 0.5:0.5:1.5;

plot_DOC_diffstarts_diffkactkinh( num_cells, n, k_act_vec, k_inh_vec, multi_weight, n_hold, multi_start, T, threshold, overlap_type)

%% DOC for the New model (canceling insted of zeros) for diff k_inh & k_act:

percent_canceled = 0.05;
threshold = 0;

k_act_vec = 0:0.5:1;
k_inh_vec = 0.5:0.5:1.5;

plot_DOC_New_diffkactkinh( num_cells, n, k_act_vec, k_inh_vec, multi_weight, percent_canceled, multi_start, T, threshold, overlap_type)


%% DOC2?

DOC2( M_Null)
title('Null')
DOC2(M_diff_starts)
title('Diff Starts')
DOC2(M_New)
title('New')

%% DOC_on_unhold

span = 0.1;
threshold = 1e-3;
threshold_Mnew = 0;

figure;
DOC_on_unhold( M_Null, threshold, span, overlap_type, n_hold)
DOC_on_unhold(M_diff_starts, threshold, span, overlap_type, n_hold)
DOC_on_unhold(M_diff_weight, threshold, span, overlap_type, n_hold)
DOC_on_unhold(M_New, threshold_Mnew, span, overlap_type, n_hold)
legend('Null','Diff Starts','Diff Starts and weights','New (no threshold!)')
title(strcat('DOC on UnHold nodes with threshold = ',num2str(threshold)));

%% DOC_on_hold

span = 0;
threshold = 1e-3;
threshold_Mnew = 0;

figure;
DOC_on_hold( M_Null, threshold, span, overlap_type, n_hold)
DOC_on_hold(M_diff_starts, threshold, span, overlap_type, n_hold)
DOC_on_hold(M_diff_weight, threshold, span, overlap_type, n_hold)
DOC_on_hold(M_New, threshold_Mnew, span, overlap_type, n_hold)
legend('Null','Diff Starts','Diff Starts and weights','New (no threshold!)')
title(strcat('DOC on the nodes that on hold with threshold = ',num2str(threshold)));

%% DOC on netowrks with noise only on part of the weights

rlowess_span = 0.2;
threshold_vec = 1e-3;

figure;
DOC_partweight_curve = DOC(M_diff_partweight, threshold_vec, rlowess_span, overlap_type );
DOC_diffpartweight_curve = DOC(M_diff_diffpartweight, threshold_vec, rlowess_span, overlap_type );
legend('same partweight','diff partweight')

%% DOC for M_sampl (cells with noise only because of sampling with the same sequencing depth)

figure;
span = 0.2; threshold = 1e-3;
DOC_curve = DOC( M_sampl, threshold, span, overlap_type);

%% DOC for M_diff_depth (cells with noise only because different sequencing depth)

figure;
span = 0.2; threshold = 1e-3;
DOC_curve = DOC( M_diff_depth, threshold, span, overlap_type);

%% %%%%%%

%% Networks 

%% corr network (not final)

pval_threshold = 0;
rho_threshold = 0;
num_links_Null = corr_network( M_Null, pval_threshold, rho_threshold );
num_links_diff_starts = corr_network( M_diff_starts, pval_threshold, rho_threshold );
num_links_diff_weight = corr_network( M_diff_weight, pval_threshold, rho_threshold );
num_links_New = corr_network( M_New, pval_threshold, rho_threshold );

%% Orr Levy

genes_minimum_activity = 50;
relevant_genes = sum(M_diff_starts>0)>genes_minimum_activity;
data = M_diff_starts(:,relevant_genes);
[number_of_non_zero_links_net Net_corr strct rand_corr net_params ] = Network_Consistency(data',[100 0.05 0.001 100],{'strings'})

