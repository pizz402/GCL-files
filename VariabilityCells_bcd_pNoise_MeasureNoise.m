function [ VariabilityCellsMatrix, bcdMatrix  ] = VariabilityCells_bcd_pNoise_MeasureNoise( num_cells, A, B, T, p_vec, sigma_noise_vec, num_itr, multi_weight, num_holds, bcd_itr)

% n = size(A,2);
% Char_plot_YorN = 'Y';% histogram plots
Char_plot_YorN = 'N';% no histogram plots

bcdMatrix = zeros(length(sigma_noise_vec),length(p_vec));
VariabilityCellsMatrix = zeros(length(sigma_noise_vec),length(p_vec));

index_sigma_noise = 1;

for sigma_noise = sigma_noise_vec
    
    sigma_noise
    
    index_p = 1;
    
    for p = p_vec
        
        p
        
        for i = 1:num_itr
            
            M_NoiseInterSelfWeight = M_NoiseInterSelfWeight_function( num_cells, A, B, T, p, multi_weight, num_holds);
             
            %normal noise on results
            M_NoiseInterSelfWeight_withnoise = M_NoiseInterSelfWeight.*normrnd(1,sigma_noise,size(M_NoiseInterSelfWeight));
            M_NoiseInterSelfWeight_withnoise(M_NoiseInterSelfWeight_withnoise < 0) = 0;
            M_NoiseInterSelfWeight_withnoise = M_NoiseInterSelfWeight_withnoise./sum(M_NoiseInterSelfWeight_withnoise,2);
        
            hist_prop = plot_histogram_corr_Spearman( M_NoiseInterSelfWeight_withnoise, Char_plot_YorN );
%             CorrCellsMatrix(index_sigma_noise,index_sigma_weights) = CorrCellsMatrix(index_sigma_noise,index_sigma_weights)+ (hist_prop.mean_hist)/num_itr;
            VariabilityCellsMatrix(index_sigma_noise,index_p) = VariabilityCellsMatrix(index_sigma_noise,index_p)+ (1-hist_prop.mean_hist)/num_itr;
            
            bcd = new_bcdistcorr_itr(M_NoiseInterSelfWeight_withnoise, bcd_itr);
            bcdMatrix(index_sigma_noise,index_p) = bcdMatrix(index_sigma_noise,index_p) + bcd/num_itr;
        
        end
        
        index_p = index_p+1;
        
    end
    
    index_sigma_noise = index_sigma_noise+1;
    
end

%%
figure;
imagesc(p_vec, sigma_noise_vec, bcdMatrix)
set(gca,'YDir','normal')
xlabel('p')
ylabel('\sigma_{Measure}')
title('GCL')
colorbar;

%%
figure;
imagesc(p_vec, sigma_noise_vec, VariabilityCellsMatrix)
set(gca,'YDir','normal')
xlabel('p')
ylabel('\sigma_{Measure}')
title('Cell-to-Cell Variability')
colorbar;

end