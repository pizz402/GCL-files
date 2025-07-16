function [COC, bcR_data, mean_corr] = MeanCorr_bcR_COC_FromM(M,num_iterations)

n = size(M,2);

COC = zeros(num_iterations,1);
bcR_data = zeros(num_iterations,1);

for n_it = 1 : num_iterations

    COC(n_it) = corrAcorrB(M, 'Spearman');
    bcR_data(n_it) = new_bcdistcorr(M);
    
end

% [rho, ~] = corr(M,'Type','Spearman');

[~,index] = sort(mean(M));%sort(M(:,1));
[rho, ~] = corr(M(:,index),'Type','Spearman');

mean_corr = sum(sum(abs(triu(rho,1))))/(n*(n-1)/2);

% figure;
% imagesc(rho);
% colorbar

end

