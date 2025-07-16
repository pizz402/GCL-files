%% For every M - plot DOC  
threshold_vec = 1e-3;

%%
figure;
DOC_diff_starts_curve = DOC(M_diff_starts, threshold_vec, rlowess_span, overlap_type );
title('Different starting conditions')
set(gca,'FontSize',32)

%%
figure;
DOC_diff_starts_curve = DOC(M_diff_weight, threshold_vec, rlowess_span, overlap_type );
title('Different starting conditions with noise on weights')
set(gca,'FontSize',32)

%%
figure;
DOC_diff_starts_curve = DOC(M_diff_diffpartweight, threshold_vec, rlowess_span, overlap_type );
title('Different starting conditions with noise on different weights')
set(gca,'FontSize',32)

%%
figure;
DOC_diff_starts_curve = DOC(M_Null, threshold_vec, rlowess_span, overlap_type );
title('The same cell with measuring noises')
set(gca,'FontSize',32)

%%
figure;
DOC_diff_starts_curve = DOC(M_sampl, threshold_vec, rlowess_span, overlap_type );
title('The same cell with noises from sampling')
set(gca,'FontSize',32)

%% 
figure;
DOC_diff_starts_curve = DOC(M_New, 0, rlowess_span, overlap_type );
title('Different starting conditions with cancelling without holding')
set(gca,'FontSize',32)
