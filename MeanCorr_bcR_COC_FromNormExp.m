function [ stat_mean_corr , stat_bcR, stat_COC ] = MeanCorr_bcR_COC_FromNormExp(num_cells, n, D_vec, num_iterations ,YNstring)
    
    mean_corr = zeros(length(D_vec),num_iterations);
    bcR = zeros(length(D_vec),num_iterations);
    COC = zeros(length(D_vec),num_iterations);

    for i = 1:length(D_vec)
        D = D_vec(i);
        for j = 1 : num_iterations 
            [ mean_corr(i,j), bcR(i,j), COC(i,j)] = MeanCorr_bcR_COC_Exp( num_cells, n, D, YNstring);
        end
    end

    stat_mean_corr.avg = mean(mean_corr,2);
    stat_mean_corr.std = std(mean_corr,0,2);

    stat_bcR.avg = mean(bcR,2);
    stat_bcR.std = std(bcR,0,2);

    stat_COC.avg = mean(COC,2);
    stat_COC.std = std(COC,0,2);
        
end

