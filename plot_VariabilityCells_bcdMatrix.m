function [] = plot_VariabilityCells_bcdMatrix(x_vec, y_vec, x_string, y_string, VariabilityCellsMatrix, bcdMatrix, contour_corr_cell, Font_Size)

% x_string = 'p' ; % y_string = 'SIGMA Measure';

size_red_contourline = 3;
size_contourlines = 0.4;

%% 1 VariabilityCellsMatrix

figure;
imagesc(x_vec, y_vec, VariabilityCellsMatrix)
set(gca,'YDir','normal')
xlabel(x_string)
ylabel(y_string)
title('Cell-to-Cell Variability')
colorbar;
cl = caxis;

%%

% figure
hold on
contour(x_vec, y_vec, VariabilityCellsMatrix,[contour_corr_cell,contour_corr_cell],'r','LineWidth',size_red_contourline);
caxis(cl);
set(gca,'FontSize',Font_Size)

%%

% figure
contour(x_vec, y_vec, VariabilityCellsMatrix,'k--','LineWidth',size_contourlines,'ShowText','on','LabelSpacing',500);
set(gca,'FontSize',Font_Size)


%% 2 bcdMatrix

figure;
imagesc(x_vec, y_vec, bcdMatrix)
set(gca,'YDir','normal')
xlabel(x_string)
ylabel(y_string)
title('GCL')
colorbar;
cl = caxis;

%%

% figure
hold on
contour(x_vec, y_vec, VariabilityCellsMatrix,[contour_corr_cell,contour_corr_cell],'r','LineWidth',size_red_contourline);
caxis(cl);
set(gca,'FontSize',Font_Size)

%%

% smoothing bcdMatrix for contourlines
h = fspecial('average');
bcdMatrix_filter = imfilter(bcdMatrix,h,'replicate');

% figure;
% imagesc(x_vec, y_vec, bcdMatrix_filter);
% set(gca,'FontSize',Font_Size)

% figure
contour(x_vec, y_vec, bcdMatrix_filter,'k--','LineWidth',size_contourlines,'ShowText','on','LabelSpacing',500);
set(gca,'FontSize',Font_Size)


end

