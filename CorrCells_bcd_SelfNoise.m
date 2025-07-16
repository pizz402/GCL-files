function [CorrCellsVec, bcdVec ] = CorrCells_bcd_SelfNoise( num_cells, A, B, T, p_vec, num_itr, multi_weight, self_inter_cases)

% n = size(A,2);
Char_plot_YorN = 'N';% no histogram plots
 
% CorrCellsVec = zeros(length(p_vec),1);
% bcdVec = zeros(length(p_vec),1);

CorrCellsVec = zeros(length(p_vec),num_itr);
bcdVec = zeros(length(p_vec),num_itr);

index_p_vec = 1;

for p = p_vec
    
    for i = 1:num_itr
        
        i
        
        M = M_NoiseSelfWeight_function( num_cells, A, B, T, p, multi_weight, self_inter_cases);
        
%         hist_prop = plot_histogram_corr_Spearman( M, Char_plot_YorN );
%         CorrCellsVec(index_p_vec) = CorrCellsVec(index_p_vec)+ (hist_prop.mean_hist)/num_itr;
% 
%         bcd = new_bcdistcorr(M);
%         bcdVec(index_p_vec) = bcdVec(index_p_vec) + bcd/num_itr;

        hist_prop = plot_histogram_corr_Spearman( M, Char_plot_YorN );
        CorrCellsVec(index_p_vec,i) = hist_prop.mean_hist;

        bcd = new_bcdistcorr(M);
        bcdVec(index_p_vec,i) = bcd;
        
    end
          
    index_p_vec = index_p_vec+1;
    
%     figure;
%     histogram(mean(M))
%     title('M example')
%     colorbar;

end

% figure;
% imagesc(M)
% title('M example')
% colorbar;

mean_CorrCells = mean(CorrCellsVec');
std_CorrCells = std(CorrCellsVec');

mean_bcd = mean(bcdVec');
std_bcd = std(bcdVec');

%%
figure;
errorbar(p_vec,mean_CorrCells,std_CorrCells)
xlabel('p')
ylabel('CORR CELLS')
title('Self Noise')

%%
figure;
errorbar(p_vec,mean_bcd,std_bcd)
xlabel('p')
ylabel('bcd')
title('Self Noise')

end