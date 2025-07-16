function [] = plot_CorrCells_bcdMatrix_CM(x_vec, y_vec, x_string, y_string, CorrCellsMatrix, bcdMatrix, CM )

% x_string = 'SIGMA Weights' ; % y_string = 'SIGMA Measure';
%% 1
figure;
imagesc(x_vec, y_vec, CorrCellsMatrix)
set(gca,'YDir','normal')
xlabel(x_string)
ylabel(y_string)
title('CORR CELLS')
colorbar;
cl = caxis;

%%
% % figure
% hold on
% contour(x_vec, y_vec, CorrCellsMatrix,[contour_corr_cell,contour_corr_cell],'r','LineWidth',3);
% caxis(cl);

%% 3
figure;
imagesc(x_vec, y_vec, bcdMatrix)
set(gca,'YDir','normal')
xlabel(x_string)
ylabel(y_string)
title('bcd')
colorbar;
cl = caxis;

%%
% % figure
% hold on
% contour(x_vec, y_vec, CorrCellsMatrix,[contour_corr_cell,contour_corr_cell],'r','LineWidth',3);
% caxis(cl);

%% 4
figure;
imagesc(x_vec, y_vec, CM)
set(gca,'YDir','normal')
xlabel(x_string)
ylabel(y_string)
title('CM')
colorbar;
cl = caxis;

%%
% % figure
% hold on
% contour(x_vec, y_vec, CorrCellsMatrix,[contour_corr_cell,contour_corr_cell],'r','LineWidth',3);
% caxis(cl);

end

