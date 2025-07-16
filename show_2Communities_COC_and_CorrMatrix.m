function [ ] = show_2Communities_COC_and_CorrMatrix(M)
% for 2 comm

span = 0.2;% unrelevant for now
string_type_corr = 'Spearman';

figure;

subplot(2,2,1)
histogram(M(:,:),50,'Normalization','probability')        
title('All cells')

a = subplot(2,2,2);
imagesc(M')
pos1 = get(a,'Position');
colorbar
set(a,'Position',pos1)
title('All cells')

subplot(2,2,3)
imagesc(corr(M,'Type',string_type_corr));
title('Corr Matrix');
colorbar; caxis([-1 1]);

subplot(2,2,4)
corrAcorrB_updown(M, span, string_type_corr)%1-2


end

