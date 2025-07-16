function [ bcdMatrix, VariabilityCellsMatrix ] = CorrCells_bcd_WeightNoise_MeasureNoise( num_cells, A, B, multi_start, T, sigma_weights_vec, sigma_noise_vec, num_itr, bcd_itr)
% VariabilityCellsMatrix and not CorrCells

% n = size(A,2);
sigma_self = 0;% no noise on B
Char_plot_YorN = 'N';% no histogram plots

% CorrCellsMatrix = zeros(length(sigma_noise_vec),length(sigma_weights_vec));
bcdMatrix = zeros(length(sigma_noise_vec),length(sigma_weights_vec));
VariabilityCellsMatrix = zeros(length(sigma_noise_vec),length(sigma_weights_vec));

index_sigma_noise = 1;

for sigma_noise = sigma_noise_vec
    
    sigma_noise
    
    index_sigma_weights = 1;
    
    for sigma_weights = sigma_weights_vec
        
        for i = 1:num_itr
            
            %uniform noise on weights
            M_diff_noiseweight_B = M_NoiseWeightAndB_function( num_cells, A, B, multi_start, T, sigma_weights, sigma_self);
            
            %normal noise on results
            M_samepart_diff_weight_withnoise = M_diff_noiseweight_B.*normrnd(1,sigma_noise,size(M_diff_noiseweight_B));
            M_samepart_diff_weight_withnoise(M_samepart_diff_weight_withnoise < 0) = 0;
            M_samepart_diff_weight_withnoise = M_samepart_diff_weight_withnoise./sum(M_samepart_diff_weight_withnoise,2);
        
            hist_prop = plot_histogram_corr_Spearman( M_samepart_diff_weight_withnoise, Char_plot_YorN );
%             CorrCellsMatrix(index_sigma_noise,index_sigma_weights) = CorrCellsMatrix(index_sigma_noise,index_sigma_weights)+ (hist_prop.mean_hist)/num_itr;
            VariabilityCellsMatrix(index_sigma_noise,index_sigma_weights) = VariabilityCellsMatrix(index_sigma_noise,index_sigma_weights)+ (1-hist_prop.mean_hist)/num_itr;
            
            bcd = new_bcdistcorr_itr(M_samepart_diff_weight_withnoise, bcd_itr);
            bcdMatrix(index_sigma_noise,index_sigma_weights) = bcdMatrix(index_sigma_noise,index_sigma_weights) + bcd/num_itr;
        
        end
        
        index_sigma_weights = index_sigma_weights+1;
        
    end
    
    index_sigma_noise = index_sigma_noise+1;
    
end


%%
% figure;
% imagesc(sigma_weights_vec, sigma_noise_vec, CorrCellsMatrix)
% set(gca,'YDir','normal')
% xlabel('\sigma_{Dynamics}')
% ylabel('\sigma_{Measure}')
% title('CORR CELLS')
% colorbar;

%%
figure;
imagesc(sigma_weights_vec, sigma_noise_vec, bcdMatrix)
set(gca,'YDir','normal')
xlabel('\sigma_{Dynamics}')
ylabel('\sigma_{Measure}')
title('GCL')
colorbar;

%%
figure;
imagesc(sigma_weights_vec, sigma_noise_vec, VariabilityCellsMatrix)
set(gca,'YDir','normal')
xlabel('\sigma_{Dynamics}')
ylabel('\sigma_{Measure}')
title('Cell-to-Cell Variability')
colorbar;

end