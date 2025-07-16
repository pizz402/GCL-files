function [] = plot_CorrCells_CorrCOC(x_vec, y_vec, x_string, y_string, CorrCellsMatrix, CorrCCCMatrix, contour_corr_cell)

% x_string = 'SIGMA Weights' ; % y_string = 'SIGMA Measure';
%%
figure;
imagesc(x_vec, y_vec, CorrCellsMatrix)
set(gca,'YDir','normal')
xlabel(x_string)
ylabel(y_string)
title('CORR CELLS')
colorbar;

%%
% figure
hold on
contour(x_vec, y_vec, CorrCellsMatrix,[contour_corr_cell,contour_corr_cell],'r','LineWidth',3);

%%
figure;
imagesc(x_vec, y_vec, CorrCCCMatrix)
set(gca,'YDir','normal')
xlabel(x_string)
ylabel(y_string)
title('CORR COC')
colorbar;

%%
% figure
hold on
contour(x_vec, y_vec, CorrCellsMatrix,[contour_corr_cell,contour_corr_cell],'r','LineWidth',3);
% A contour plot displays isolines
% figure;
% contour(sigma_weights_vec, sigma_noise_vec, CorrCCCMatrix,'r','LineWidth',3);

%%


end

