function [ ] = show_CCC_and_CorrMatrix(A, M)
% for 3 comm

n = size(M,2);
string_type_corr = 'Spearman';

figure;

subplot(2,3,1)
imagesc(A.matrix)
title('Activation Matrix');

subplot(2,3,2)
histogram(M(:,:),50,'Normalization','probability')        
title('All cells')

a = subplot(2,3,3);
imagesc(M')
pos1 = get(a,'Position');
colorbar
set(a,'Position',pos1)
title('All cells')

subplot(2,3,4)
imagesc(corr(M,'Type',string_type_corr));
title('Corr Matrix');
colorbar; caxis([-1 1]);

subplot(2,3,5)
corrAcorrB_updown(M(:,1:2*n/3), string_type_corr)%1-2
 
% corrAcorrB_updown(M(:,(n/3+1):end), rlowess_span, string_type_corr)%2-3

subplot(2,3,6)
corrAcorrB_updown([M(:,1:n/3),M(:,(2*n/3+1):end)], string_type_corr)%1-3


end

