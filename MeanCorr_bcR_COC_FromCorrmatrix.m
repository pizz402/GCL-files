function [ stat_mean_corr , stat_bcR, stat_COC ] = MeanCorr_bcR_COC_FromCorrmatrix(n, CorrMatrixStd, num_cells_vec, num_iterations)
     
    mean_corr = zeros(length(num_cells_vec),num_iterations);
    bcR = zeros(length(num_cells_vec),num_iterations);
    COC = zeros(length(num_cells_vec),num_iterations);

    for n_it = 1 : num_iterations
    
        for n_k = 1:length(num_cells_vec)

            num_cells = num_cells_vec(n_k);
            
            M = data_from_corrmatrix(n, CorrMatrixStd, num_cells);
         
            COC(n_k,n_it) = corrAcorrB_updown(M, 'Spearman');
            bcR(n_k,n_it) = bcdistcorr(M(:,1:n/2),M(:,(n/2+1):n));

            [rho, ~] = corr(M,'Type','Spearman');
            %mean matrix corr:
            mean_corr(n_k,n_it) = sum(sum(abs(triu(rho,1))))/(n*(n-1)/2);

        end
        
    end
    
    stat_mean_corr.avg = mean(mean_corr,2);
    stat_mean_corr.std = std(mean_corr,0,2);
    
    stat_bcR.avg = mean(bcR,2);
    stat_bcR.std = std(bcR,0,2);

    stat_COC.avg = mean(COC,2);
    stat_COC.std = std(COC,0,2);
        
end

