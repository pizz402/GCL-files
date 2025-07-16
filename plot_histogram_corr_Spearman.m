function [ hist_prop ] = plot_histogram_corr_Spearman( M, Char_plot_YorN )

string_type_corr = 'Spearman';

[rho, ~] = corr(M','Type',string_type_corr);
% [imax,~] = size(rho);

[ii,jj] = meshgrid(1:size(M,1),1:size(M,1));
rho_vector = rho(jj>ii);

% rho_vector = zeros(1,imax*(imax-1)/2);
% k = 0;
% for i = 1:imax-1
%     for j = i+1:imax
%         k = k + 1;
%         rho_vector(k) = rho(i,j);
%     end
% end


hist_prop.median_hist = median(rho_vector);
hist_prop.mean_hist = mean(rho_vector);

if Char_plot_YorN == 'Y'
histogram(rho_vector,50,'Normalization','probability');
% histogram(triu_rho_vector,50)
hold on
end

end

