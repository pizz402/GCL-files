function [] = plot_CorrAcorrBUpdown_SigmaWeights(k_act, sigma_weights_vec, percent_links_unfixed_betweenComm, n, multi_weight, num_cells, multi_start, T, num_itr, rlowess_span)

string_type_corr = 'Spearman';

s = 1;
corrAB_vec = zeros(1,length(sigma_weights_vec));

for sigma_weights = sigma_weights_vec
    
    for i = 1:num_itr
        
        M_communities_diff_links = M_communities_diff_links_function( k_act, sigma_weights, percent_links_unfixed_betweenComm, n, multi_weight, num_cells, multi_start, T);
        [corrAB, pvalAB] = corrAcorrB_updown(M_communities_diff_links, rlowess_span, string_type_corr);
        corrAB_vec(s) = corrAB_vec(s)+ corrAB/num_itr;
        
    end
    
    s = s+1;
    
end

figure;
plot(sigma_weights_vec, corrAB_vec);

end

