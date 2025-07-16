function [ stat_mean_corr_weights , stat_mean_corr_meas, stat_bcR_weights, stat_bcR_meas, stat_COC_weights, stat_COC_meas ] = MeanCorr_bcR_COC_FromSim(num_cells_vec, A, B, A2, n_hold_zero, multi_start, T, sigma_weights, fraction_noise, sigma_results, num_iterations)
     
    n = size(A.matrix,2);

    for n_it = 1 : num_iterations
    
        for n_k = 1:length(num_cells_vec)

            num_cells = num_cells_vec(n_k);

            M_weights = M_diff_samepartweight_function_ABA2( num_cells, A, B, A2, n_hold_zero, multi_start, T, sigma_weights, fraction_noise);
            M_meas = M_Null_function_ABA2( num_cells, A, B, A2, n_hold_zero, multi_start, T, sigma_results);  

%             COC_weights(n_k,n_it) = corrAcorrB_updown(M_weights, 'Spearman');
%             COC_meas(n_k,n_it) = corrAcorrB_updown(M_meas, 'Spearman');
 
%             bcR_weights(n_k,n_it) = bcdistcorr(M_weights(:,1:n/2),M_weights(:,(n/2+1):n));
%             bcR_meas(n_k,n_it) = bcdistcorr(M_meas(:,1:n/2),M_meas(:,(n/2+1):n));
            
            COC_weights(n_k,n_it) = corrAcorrB(M_weights, 'Spearman');
            COC_meas(n_k,n_it) = corrAcorrB(M_meas, 'Spearman');
            
            bcR_weights(n_k,n_it) = new_bcdistcorr(M_weights);
            bcR_meas(n_k,n_it) = new_bcdistcorr(M_meas);

            [rho_weights, ~] = corr(M_weights,'Type','Spearman');
            [rho_meas, ~] = corr(M_meas,'Type','Spearman');
            
            %mean matrix corr:
            mean_corr_weights(n_k,n_it) = sum(sum(abs(triu(rho_weights,1))))/(n*(n-1)/2);
            mean_corr_meas(n_k,n_it) = sum(sum(abs(triu(rho_meas,1))))/(n*(n-1)/2);

        end
        
    end
    
    stat_mean_corr_weights.avg = mean(mean_corr_weights,2);
    stat_mean_corr_weights.std = std(mean_corr_weights,0,2);
    
    stat_mean_corr_meas.avg = mean(mean_corr_meas,2);
    stat_mean_corr_meas.std = std(mean_corr_meas,0,2);
    
    stat_bcR_weights.avg = mean(bcR_weights,2);
    stat_bcR_weights.std = std(bcR_weights,0,2);
    
    stat_bcR_meas.avg = mean(bcR_meas,2);
    stat_bcR_meas.std = std(bcR_meas,0,2);
    
    stat_COC_weights.avg = mean(COC_weights,2);
    stat_COC_weights.std = std(COC_weights,0,2);
    
    stat_COC_meas.avg = mean(COC_meas,2);
    stat_COC_meas.std = std(COC_meas,0,2);
    
    
end

