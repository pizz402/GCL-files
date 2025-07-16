function [CorrCellsMatrix, CorrCCCMatrix] = CorrCells_DiffLinksNoise_WeightNoise( k_act, n, multi_weight, num_cells, multi_start, T, percent_links_unfixed_vec, sigma_weights_vec, num_itr)

Char_plot_YorN = 'N';% no histogram plots
string_type_corr = 'Spearman';% for CCC
% string_type_corr = 'Pearson';
rlowess_span = 0.2;% for CCC plot

CorrCellsMatrix = zeros(length(sigma_weights_vec),length(percent_links_unfixed_vec));
CorrCCCMatrix = zeros(length(sigma_weights_vec),length(percent_links_unfixed_vec));

index_sigma_weights_noise = 1;

for sigma_weights_noise = sigma_weights_vec
    
    index_percent_links_unfixed = 1;
    
    for percent_links_unfixed = percent_links_unfixed_vec
        
        for i = 1:num_itr
            
            M_diff_links_withnoise = M_difflinks_diffweight_function( k_act, percent_links_unfixed, n, multi_weight, num_cells, multi_start, T, sigma_weights_noise);
            
            hist_prop = plot_histogram_corr_Spearman( M_diff_links_withnoise, Char_plot_YorN );
            CorrCellsMatrix(index_sigma_weights_noise,index_percent_links_unfixed) = CorrCellsMatrix(index_sigma_weights_noise,index_percent_links_unfixed)+ (hist_prop.mean_hist)/num_itr;

            [corrAB, pvalAB] = corrAcorrB(M_diff_links_withnoise, rlowess_span, string_type_corr);
            CorrCCCMatrix(index_sigma_weights_noise,index_percent_links_unfixed) = CorrCCCMatrix(index_sigma_weights_noise,index_percent_links_unfixed) + corrAB/num_itr;
            % here - sometimes NULL - need to be fixed
            
        end
        
        index_percent_links_unfixed = index_percent_links_unfixed+1;
        
    end
    
    index_sigma_weights_noise = index_sigma_weights_noise+1
    
end

%%
figure;
imagesc(percent_links_unfixed_vec, sigma_weights_vec, CorrCellsMatrix)
set(gca,'YDir','normal')
xlabel('Percent Links Un-Fixed')
ylabel('Sigma Weights Noise')
title('CORR CELLS')
colorbar;

%%
figure;
imagesc(percent_links_unfixed_vec, sigma_weights_vec, CorrCCCMatrix)
set(gca,'YDir','normal')
xlabel('Percent Links Un-Fixed')
ylabel('Sigma Weights Noise')
title('CORR CCC')
colorbar;

end