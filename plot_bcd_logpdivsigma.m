function [] = plot_bcd_logpdivsigma(x_vec, y_vec, x_string, y_string, VariabilityCellsMatrix, bcdMatrix, contour_corr_cell_vec, Font_Size)

% x_string = 'p' ; % y_string = 'SIGMA Measure';

scater_marker = ['o';'x';'s'];

size_red_contourline = 3;
length_contour_corr_cell_vec = length(contour_corr_cell_vec);

randnum_fig1 = round(rand()*100);
randnum_fig2 = round(rand()*100);

for i = 1: length_contour_corr_cell_vec
    
    contour_corr_cell = contour_corr_cell_vec(i);
    
    figure(randnum_fig1)
    [a,~] = contour(x_vec, y_vec, VariabilityCellsMatrix,[contour_corr_cell,contour_corr_cell],'r','LineWidth',size_red_contourline);
    xlabel(x_string)
    ylabel(y_string)
    title('Cell-to-Cell Variability')
    set(gca,'FontSize',Font_Size)
    hold on 
    
    figure(randnum_fig2)
    xvec_specific = a(1,2:end-1);
    yvec_specific = a(2,2:end-1);
    % VariabilityCells_specific = interp2(x_vec,y_vec,VariabilityCellsMatrix,xvec_specific,yvec_specific)
    bcd_specific = interp2(x_vec,y_vec,bcdMatrix,xvec_specific,yvec_specific);
    plot(log(xvec_specific./yvec_specific),bcd_specific,scater_marker(i),'DisplayName',strcat('ccVar =  ',num2str(contour_corr_cell_vec(i))));
    xlabel('log( p/ \sigma)')
    ylabel('GCL')
    set(gca,'FontSize',Font_Size)
    hold on 
    
end

legend('Location','best')

hold off

end

