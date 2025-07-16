function [VariabilityCellsVec, bcdVec ] = CorrCells_bcd_InterSelfNoise( num_cells, A, Bvec, T, p_vec, num_itr, multi_weight, num_holds, bcd_itr,Font_Size,Line_Width)
% VariabilityCellsVec and not CorrCells

% n = size(A,2);
Char_plot_YorN = 'N';% no histogram plots


% CorrCellsVec = zeros(length(p_vec),num_itr);
VariabilityCellsVec = zeros(length(p_vec),num_itr);
bcdVec = zeros(length(p_vec),num_itr);

index_p_vec = 1;

for p = p_vec
    
    for i = 1:num_itr
        
        i
        
        M = M_NoiseInterSelfWeight_function( num_cells, A, Bvec, T, p, multi_weight, num_holds);
        
        hist_prop = plot_histogram_corr_Spearman( M, Char_plot_YorN );
%         CorrCellsVec(index_p_vec,i) = hist_prop.mean_hist;
        VariabilityCellsVec(index_p_vec,i) = 1 - hist_prop.mean_hist;
              
        bcd = new_bcdistcorr_itr(M, bcd_itr);
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

mean_VariabilityCells = mean(VariabilityCellsVec');
std_VariabilityCells = std(VariabilityCellsVec');

mean_bcd = mean(bcdVec');
std_bcd = std(bcdVec');

%%
% figure;
% errorbar(p_vec,mean_CorrCells,std_CorrCells)
% xlabel('p')
% ylabel('CORR CELLS')
% title('Inter & Self Noise')

%%
figure;
errorbar(p_vec,mean_VariabilityCells,std_VariabilityCells,'LineWidth',Line_Width)
xlabel('p')
ylabel('Cell-to-Cell Variability')
set(gca,'FontSize',Font_Size)
% title('Inter & Self Noise')

%%
figure;
errorbar(p_vec,mean_bcd,std_bcd,'LineWidth',Line_Width)
xlabel('p')
ylabel('GCL')
set(gca,'FontSize',Font_Size)
% title('Inter & Self Noise')


end