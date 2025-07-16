function [ ] = plotall_compositionaldata(num_cells, n, D_vec, num_iterations, sp_D_vec)

num_sp_D = length(sp_D_vec);
num_rows = 5;%or 6 with COC
num_column = num_sp_D;%3?

[ stat_mean_corr , stat_bcR, stat_COC ] = MeanCorr_bcR_COC_FromNormExp(num_cells, n, D_vec, num_iterations,'N');

num_fig = 1;

figure;

ha = tight_subplot(num_rows,num_column,[0 0],[0 0],[0 0]);

for D = sp_D_vec
    
    D
    
    % cell without noise
    one_cell = rand(1,n).^(-D);
    [~,index] = sort(one_cell);

    % M matrix:
    % same cell - 
    M = repmat(one_cell,num_cells,1);
    M = M.*abs(normrnd(1,0.2.*M,[num_cells,n]));
    % normalization - 
    M = M./sum(M,2);
    
    %mean matrix corr:
%     [rho, ~] = corr(M,'Type','Spearman');
    [rho, ~] = corr(M(:,index),'Type','Spearman');
    
    subplot(num_rows,num_column,num_fig)
    imagesc(rho);
    title(strcat('D=', num2str(D)))
    colorbar
    caxis([-1 1])

    subplot(num_rows,num_column,num_fig+num_column)
    plot(sort(one_cell));
    xlim([0,n])

    subplot(num_rows,num_column,num_fig+2*num_column)
    histogram(one_cell,10,'Normalization','probability');
    ylim([0,1])
    
    num_fig = num_fig+1;

end

num_fig = num_fig + 2*num_column;

subplot(num_rows,num_column,num_fig:num_fig+num_column-1);
plot_errorbar(D_vec, stat_mean_corr, '', 'CM')
ylim([0,1])
subplot(num_rows,num_column,num_fig+num_column:num_fig+2*num_column-1);
plot_errorbar(D_vec, stat_bcR, 'D', 'bcD')
% subplot(num_rows,num_column,num_fig+2*num_column:num_fig+3*num_column-1);
% plot_errorbar(D_vec, stat_COC, 'D', 'COC')



end