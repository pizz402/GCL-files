function [CorrCellsMatrix, CorrCCCMatrix, bcdMatrix, CM ] = CorrCells_CorrCOC_bcd_CM_WeightNoise_MeasureNoise( num_cells, A, B, A2, n_hold_zero, multi_start, T, sigma_weights_vec, sigma_noise_vec, num_itr)

n = size(A.matrix,2);

fraction_noise_one = 1;% noise on *all* weights
Char_plot_YorN = 'N';% no histogram plots
string_type_corr = 'Spearman';% for CCC
% string_type_corr = 'Pearson';
rlowess_span = 0.2;% for CCC plot

index_sigma_noise = 1;

CorrCellsMatrix = zeros(length(sigma_noise_vec),length(sigma_weights_vec));
CorrCCCMatrix = zeros(length(sigma_noise_vec),length(sigma_weights_vec));
bcdMatrix = zeros(length(sigma_noise_vec),length(sigma_weights_vec));
CM = zeros(length(sigma_noise_vec),length(sigma_weights_vec));

for sigma_noise = sigma_noise_vec
    
    index_sigma_weights = 1;
    
    for sigma_weights = sigma_weights_vec
        
        for i = 1:num_itr
            
            %uniform noise on weights
            M_samepart_diff_weight = M_diff_samepartweight_function_ABA2( num_cells, A, B, A2, n_hold_zero, multi_start, T, sigma_weights, fraction_noise_one);
            
            %normal noise on results
            M_samepart_diff_weight_withnoise = M_samepart_diff_weight.*normrnd(1,sigma_noise,size(M_samepart_diff_weight));
            M_samepart_diff_weight_withnoise(M_samepart_diff_weight_withnoise < 0) = 0;
            M_samepart_diff_weight_withnoise = M_samepart_diff_weight_withnoise./sum(M_samepart_diff_weight_withnoise,2);

            hist_prop = plot_histogram_corr_Spearman( M_samepart_diff_weight_withnoise, Char_plot_YorN );
            CorrCellsMatrix(index_sigma_noise,index_sigma_weights) = CorrCellsMatrix(index_sigma_noise,index_sigma_weights)+ (hist_prop.mean_hist)/num_itr;

            [corrAB, pvalAB] = corrAcorrB(M_samepart_diff_weight_withnoise, string_type_corr);
            CorrCCCMatrix(index_sigma_noise,index_sigma_weights) = CorrCCCMatrix(index_sigma_noise,index_sigma_weights) + corrAB/num_itr;
            
            bcd = new_bcdistcorr(M_samepart_diff_weight_withnoise);
            bcdMatrix(index_sigma_noise,index_sigma_weights) = bcdMatrix(index_sigma_noise,index_sigma_weights) + bcd/num_itr;

            [rho, ~] = corr(M_samepart_diff_weight_withnoise,'Type','Spearman');
            mean_corr = sum(sum(abs(triu(rho,1))))/(n*(n-1)/2);
            CM(index_sigma_noise,index_sigma_weights) = CM(index_sigma_noise,index_sigma_weights)+mean_corr/num_itr;

        
        end
        
        index_sigma_weights = index_sigma_weights+1;
        
    end
    
    index_sigma_noise = index_sigma_noise+1;
    
end

%%
figure;
imagesc(sigma_weights_vec, sigma_noise_vec, CorrCellsMatrix)
set(gca,'YDir','normal')
xlabel('Sigma Weights')
ylabel('Sigma Measure')
title('CORR CELLS')
colorbar;

%%
figure;
imagesc(sigma_weights_vec, sigma_noise_vec, CorrCCCMatrix)
set(gca,'YDir','normal')
xlabel('Sigma Weights')
ylabel('Sigma Measure')
title('CORR CCC')
colorbar;

%%
figure;
imagesc(sigma_weights_vec, sigma_noise_vec, bcdMatrix)
set(gca,'YDir','normal')
xlabel('Sigma Weights')
ylabel('Sigma Measure')
title('bcd')
colorbar;

%%
figure;
imagesc(sigma_weights_vec, sigma_noise_vec, CM)
set(gca,'YDir','normal')
xlabel('Sigma Weights')
ylabel('Sigma Measure')
title('CM')
colorbar;
end