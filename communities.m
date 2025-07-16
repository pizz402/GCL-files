%% Communities 
% ctrl+enter on the fixed parameters from pdf4

%% Activation communities - k = 2 and flips only between the communities

%% 
k_act = 2; sigma_weights = 0.9;

%% 1
percent_links_unfixed_betweenComm = 0.2;
M_communities_diff_links = M_communities_diff_links_function( k_act, sigma_weights, percent_links_unfixed_betweenComm, n, multi_weight, num_cells, multi_start, T);

%% 2
string_type_corr = 'Spearman';
figure;
corrAcorrB_updown(M_communities_diff_links, rlowess_span, string_type_corr)

%% 3
M_multiplot_ABA2_updown( M_communities_diff_links, threshold, n_hold_zero, rlowess_span, overlap_type)

%% 4 - IM HERE
show_2Communities_COC_and_CorrMatrix(M_communities_diff_links)

%% 5 - IM HERE
tic

num_itr = 5;
k_act = 2; sigma_weights = 0.2;
percent_links_unfixed_betweenComm_vec = 0:0.1:1;

plot_CorrAcorrBUpdown_PercentLinksUnfixedBetweenComm(k_act, sigma_weights, percent_links_unfixed_betweenComm_vec, n, multi_weight, num_cells, multi_start, T, num_itr, rlowess_span)

toc

%% 6 - IM HERE
tic

num_itr = 5;
k_act = 2; percent_links_unfixed_betweenComm = 0.3;
sigma_weights_vec = 0.1:0.1:0.5;

plot_CorrAcorrBUpdown_SigmaWeights(k_act, sigma_weights_vec, percent_links_unfixed_betweenComm, n, multi_weight, num_cells, multi_start, T, num_itr, rlowess_span)

toc






%%
%% Activation communities - k_iner & k_intra & k_intra_C
%% 1
k_inter = 2;
k_intra = 2;
k_intra_C = 2;
% 2 Communities
[ A, B, A2 ] = Build_Network_activation_communities( k_inter, k_intra, n, multi_weight );
% 3 Communities
[ A, B, A2 ] = Build_Network_activation_3communities( k_inter, k_intra, k_intra_C, n, multi_weight );
%% 2 - CCC for 2 communities for different k_intra and k_inter
k_inter_vec = 1:0.5:2; k_intra_vec = 1:0.5:2;
sigma_weights = 0.4;
CorrMatrix = corr_from_ccc(num_cells, k_inter_vec, k_intra_vec, multi_weight, n, n_hold_zero, multi_start, T, sigma_weights);
%% 3 - CCC for 3 communities % for special A
show_CCC_and_CorrMatrix( A, M_diff_weight)

