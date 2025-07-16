%% for_pdf3 

%% %%%%%%
%% comaring with/without C matrix

%% 

figure;
histogram(M_diff_starts(:,:),50,'Normalization','probability')
hold on
histogram(new_M_diff_starts(:,:),50,'Normalization','probability')
legend('different x0', 'new different x0')

%%
% sum(M_diff_starts')% ones(1,num_cells)
% sum(new_M_diff_starts')% ones(1,num_cells)

%%

figure;
histogram(M_diff_starts(:,1:n_hold),50,'Normalization','probability')
figure;
histogram(new_M_diff_starts(:,1:n_hold),50,'Normalization','probability')

%%
% length(find(M_diff_starts(:,1:n_hold))) %~(n_hold/2)*num_cells
% length(find(new_M_diff_starts(:,1:n_hold))) %~(n_hold/2)*num_cells

%% %%%%%
%% DOC and CCC on M with/without C matrix

%% CCC - example

figure;
corrAcorrB(M_Null, rlowess_span)

%% DOC and CCC without C

threshold_vec = 1e-3;
figure;
corrAcorrB(M_diff_starts, rlowess_span)
figure;
DOC(M_diff_starts, threshold_vec, rlowess_span, overlap_type );

%% DOC and CCC with C

threshold_vec = 1e-5;
figure;
corrAcorrB(new_M_diff_starts, rlowess_span)
figure;
DOC(new_M_diff_starts, threshold_vec, rlowess_span, overlap_type );

%% DOC on M with C - with Jaccard ovlerlap

figure;
DOC(new_M_diff_starts, threshold_vec, rlowess_span, 2 );

%% CCC on unhold genes

figure;
corrAcorrB_unhold(M_diff_starts, n_hold, rlowess_span)
figure;
corrAcorrB_unhold(new_M_diff_starts, n_hold, rlowess_span)

%% DOC on unhold genes

figure;
DOC_on_unhold( M_diff_starts, threshold_vec, rlowess_span, overlap_type, n_hold)%, Dis_type)
figure;
DOC_on_unhold( new_M_diff_starts, threshold_vec, rlowess_span, overlap_type, n_hold)%, Dis_type)

%% 

figure;
plot_histogram_corr_Spearman( M_Null )
plot_histogram_corr_Spearman( M_diff_starts )
plot_histogram_corr_Spearman( new_M_diff_starts )
legend('Null', 'different x0', 'new diff starts')


%% For every M - example of cell, all cells, spearman and DOC  
%%
without_C = 0;% without C
figure;
% figure('units','normalized','outerposition',[0 0 1 1])
M_multiplot( M_diff_starts, threshold_vec, num_cells, n, n_hold, k_act, k_inh ,p_2inter, rlowess_span, overlap_type, without_C )
stringwithoutC_for_png = strcat('results genes/withoutC',num2str(round(1000*rand)),'threshold',num2str(threshold_vec),'cells',num2str(num_cells),'genes',num2str(n), 'hold', num2str(n_hold),'act',num2str(k_act),'inh',num2str(k_inh),'p2inter',num2str(p_2inter) );
stringwithoutC_for_png =  strcat(strrep(stringwithoutC_for_png,'.',''),'.png');
% saveas(gcf,stringwithoutC_for_png)

%%
with_C = 1;% with C
figure;
% figure('units','normalized','outerposition',[0 0 1 1])
M_multiplot( new_M_diff_starts, threshold_vec, num_cells, n, n_hold, k_act, k_inh, p_2inter, rlowess_span, overlap_type, with_C )
stringwithC_for_png = strcat('results genes/withC',num2str(round(1000*rand)),'threshold',num2str(threshold_vec),'cells',num2str(num_cells),'genes',num2str(n), 'hold', num2str(n_hold),'act',num2str(k_act),'inh',num2str(k_inh),'p2inter',num2str(p_2inter) );
stringwithC_for_png =  strcat(strrep(stringwithC_for_png,'.',''),'.png');
% saveas(gcf,stringwithC_for_png)


%% multuplot for different strength of 2-activations

with_C = 1;% with C

for p_2inter = 0.005:0.005:0.01
    
    C = cat(3, full(multi_weight*sprand(n,n,p_2inter)), full(multi_weight*sprand(n,n,p_2inter)));
    C(1:n+1:(end/2)) = 0 ; C(1+(end/2):n+1:end) = 0 ;
    
    new_M_diff_starts = new_M_diff_starts_function( num_cells, A, B, n_hold, multi_start, T, C);
    
    figure('units','normalized','outerposition',[0 0 1 1])
    M_multiplot( new_M_diff_starts, threshold_vec, num_cells, n, n_hold, k_act, k_inh, p_2inter, rlowess_span, overlap_type, with_C )
    
end
