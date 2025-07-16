function [CorrCellsMatrix, CorrCCCMatrix, bcdMatrix, CM ] = CorrCells_CorrCOC_bcd_CM_WeightNoise_BNoise( num_cells, A, B, multi_start, T, sigma_weights_vec, sigma_self_vec, num_itr)

n = size(A,2);

Char_plot_YorN = 'N';% no histogram plots
% string_type_corr = 'Spearman';% for CCC
% string_type_corr = 'Pearson';

CorrCellsMatrix = zeros(length(sigma_self_vec),length(sigma_weights_vec));
CorrCCCMatrix = zeros(length(sigma_self_vec),length(sigma_weights_vec));
bcdMatrix = zeros(length(sigma_self_vec),length(sigma_weights_vec));
CM = zeros(length(sigma_self_vec),length(sigma_weights_vec));

index_sigma_self = 1;

for sigma_self = sigma_self_vec
    
%     sigma_self
    index_sigma_weights = 1;
    
    for sigma_weights = sigma_weights_vec
        
%         sigma_weights
        for i = 1:num_itr
            
            %uniform noise on weights & B
            M = M_NoiseWeightAndB_function( num_cells, A, B, multi_start, T, sigma_weights, sigma_self);
           
            hist_prop = plot_histogram_corr_Spearman( M, Char_plot_YorN );
            CorrCellsMatrix(index_sigma_self,index_sigma_weights) = CorrCellsMatrix(index_sigma_self,index_sigma_weights)+ (hist_prop.mean_hist)/num_itr;

%             [corrAB, ~] = corrAcorrB(M, string_type_corr);
%             CorrCCCMatrix(index_sigma_self,index_sigma_weights) = CorrCCCMatrix(index_sigma_self,index_sigma_weights) + corrAB/num_itr;
            
            bcd = new_bcdistcorr(M);
            bcdMatrix(index_sigma_self,index_sigma_weights) = bcdMatrix(index_sigma_self,index_sigma_weights) + bcd/num_itr;

%             [rho, ~] = corr(M,'Type','Spearman');
%             mean_corr = sum(sum(abs(triu(rho,1))))/(n*(n-1)/2);
%             CM(index_sigma_self,index_sigma_weights) = CM(index_sigma_self,index_sigma_weights)+mean_corr/num_itr;

        
        end
        
        index_sigma_weights = index_sigma_weights+1;
        
    end
    
    index_sigma_self = index_sigma_self+1;
    
end

% figure;
% imagesc(M)
% title('M example')
% colorbar;

%%
figure;
imagesc(sigma_weights_vec, sigma_self_vec, CorrCellsMatrix)
set(gca,'YDir','normal')
xlabel('Sigma Weights')
ylabel('Sigma Self')
title('CORR CELLS')
colorbar;

% %%
% figure;
% imagesc(sigma_weights_vec, sigma_self_vec, CorrCCCMatrix)
% set(gca,'YDir','normal')
% xlabel('Sigma Weights')
% ylabel('Sigma Self')
% title('CORR CCC')
% colorbar;

%%
figure;
imagesc(sigma_weights_vec, sigma_self_vec, bcdMatrix)
set(gca,'YDir','normal')
xlabel('Sigma Weights')
ylabel('Sigma Self')
title('bcd')
colorbar;

% %%
% figure;
% imagesc(sigma_weights_vec, sigma_self_vec, CM)
% set(gca,'YDir','normal')
% xlabel('Sigma Weights')
% ylabel('Sigma Self')
% title('CM')
% colorbar;
end