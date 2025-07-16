function [COC] = corrAcorrB_updown_hist(M_matrix, string_type_corr, num_itr_hist)

for i = 1:num_itr_hist
    
    COC(i) = corrAcorrB_updown(M_matrix, string_type_corr);
    
end

end

