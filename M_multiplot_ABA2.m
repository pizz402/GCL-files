function [  ] = M_multiplot_ABA2( M, threshold, n_hold, rlowess_span, overlap_type)
% ( M, threshold, num_cells, n, n_hold, k_act, k_inh, p_A2, rlowess_span, overlap_type)

num_cells = size(M,1);
n = size(M,2);

string_type_corr = 'Spearman';
% string_type_corr = 'Pearson';

if n_hold==0
    
    figure;

    subplot(2,3,1)
    histogram(M(1,:),50,'Normalization','probability')        
    title('Example of cell')

    subplot(2,3,2)
    histogram(M(:,:),50,'Normalization','probability')        
    title('All cells')

    subplot(2,3,3)
    new_plot_histogram_corr_Spearman( M, n_hold )
    title('Spearman between cells')

    a = subplot(2,3,4);
    imagesc(M')
    pos1 = get(a,'Position');
    colorbar
    set(a,'Position',pos1)
    title('All cells')

    subplot(2,3,5)
    DOC(M, threshold, rlowess_span, overlap_type );
    title('DOC')

    subplot(2,3,6)
    [corrAB, pvalAB] = corrAcorrB(M, string_type_corr,'Y');
%     str = {strcat('corr = ',num2str(corrAB)),strcat('pval = ',num2str(pvalAB))};
%     text(0.1,0.8,str)
    str = strcat(' COC = ',num2str(corrAB),' pval = ',num2str(pvalAB));
    title(str)

else
    
    M_unhold = M(:,(n_hold+1):end);

    figure;

    subplot(2,4,1)
    histogram(M(1,:),50,'Normalization','probability')        
    title('Example of cell')

    subplot(2,4,2)
    histogram(M(:,:),50,'Normalization','probability')        
    title('All cells')

    subplot(2,4,3)
    new_plot_histogram_corr_Spearman( M, n_hold )
    title('Spearman between cells')

    a = subplot(2,4,4);
    imagesc(M')
    pos1 = get(a,'Position');
    colorbar
    set(a,'Position',pos1)
    title('All cells')

    subplot(2,4,5)
    DOC(M, threshold, rlowess_span, overlap_type );
    title('DOC')

    subplot(2,4,6)
    DOC(M_unhold, threshold, rlowess_span, overlap_type );
    title('DOC on unhold nodes')

    subplot(2,4,7)
    [corrAB, pvalAB] = corrAcorrB(M, string_type_corr);
    str = strcat('corr = ',num2str(corrAB),', pval = ', num2str(pvalAB));
    dim = [.5 .5 0 0];
    annotation('textbox',dim,'String',str,'FitBoxToText','on');
    title('CorrCorr')

    subplot(2,4,8)
    [corrAB, pvalAB] = corrAcorrB(M_unhold, string_type_corr);
    str = {strcat('corr = ',num2str(corrAB)),strcat('pval = ',num2str(pvalAB))};
    text(0.1,0.8,str)
    title('CorrCorr on unhold nodes')

end

suptitle(strcat('cells=',num2str(num_cells),', genes=',num2str(n), ', hold =', num2str(n_hold)));
% suptitle(strcat('cells=',num2str(num_cells),', genes=',num2str(n), ', hold =', num2str(n_hold),', act =',num2str(k_act),', inh =',num2str(k_inh) ,', 2inter =',num2str(p_A2) , ', threshold = ',num2str(threshold)));

end

