function [] = plot_CorrAcorrBUpdown_PercentLinksUnfixedBetweenComm(k_act, sigma_weights, percent_links_unfixed_betweenComm_vec, n, multi_weight, num_cells, multi_start, T, num_itr, rlowess_span)

string_type_corr = 'Spearman';

p = 1;
corrAB_vec = zeros(1,length(percent_links_unfixed_betweenComm_vec));

for percent_links_unfixed_betweenComm = percent_links_unfixed_betweenComm_vec
    
    for i = 1:num_itr
        
        M_communities_diff_links = M_communities_diff_links_function( k_act, sigma_weights, percent_links_unfixed_betweenComm, n, multi_weight, num_cells, multi_start, T);
        [corrAB, pvalAB] = corrAcorrB_updown(M_communities_diff_links, rlowess_span, string_type_corr);
        corrAB_vec(p) = corrAB_vec(p)+ corrAB/num_itr;
        
    end
    
    p = p+1;
    
end

figure;
plot(percent_links_unfixed_betweenComm_vec,corrAB_vec);

end

