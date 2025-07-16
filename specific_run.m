function [] = specific_run()
% ORR LEVY - 
% M_switching_links = the data 
% note: number of cells = size(M,1) & number of genes = size(M,2)
% percent_damage_vec = vector size 4 for the 4 ages (can be [1,2,3,4])
% num_itr_hist = 400 (for the histogram of COC)
% rlowess_span = 0 (not used for now)
% string_type_corr = 'Spearman'
% font_size = 12

figure;

for i = 1:length(percent_damage_vec)%1:4
    
    p = percent_damage_vec(i);
    % ORR LEVY - next line in comment because in your case it is the data
%     M_switching_links = M_damagelinks_completelydiffweight_function( k_act, p, n, multi_weight, num_cells, multi_start, T);
    
    subplot(2,4,i)
    [corrAB, pvalAB] = corrAcorrB(M_switching_links, rlowess_span, string_type_corr,1);%plot_YN = 1 for ploting 
    axis square
    str_title = num2str(p);
    title(str_title)
    box on
    set(gca,'fontsize', font_size);
        
    subplot(2,4,4+i)
    COC = corrAcorrB_hist(M_switching_links, rlowess_span, string_type_corr,num_itr_hist);
    mean_hist = mean(COC);
    histogram(COC,'Normalization','probability');
    axis square
    hold on
    line([mean_hist, mean_hist], ylim, 'LineWidth', 2, 'Color', 'r','LineStyle','--');    
    set(gca,'fontsize', font_size);
    
end

end