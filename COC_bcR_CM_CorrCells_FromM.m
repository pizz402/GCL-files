function [COC_prop, bcR_prop, MC, CorrCells] = COC_bcR_CM_CorrCells_FromM(M,num_iterations)

n = size(M,2);
num_cells = size(M,1);

COC = zeros(num_iterations,1);
bcR = zeros(num_iterations,1);

for n_it = 1 : num_iterations

    COC(n_it) = corrAcorrB(M, 'Spearman');
    bcR(n_it) = new_bcdistcorr(M);
    
end

COC_prop.mean = mean(COC);
COC_prop.std = std(COC);

bcR_prop.mean = mean(bcR);
bcR_prop.std = std(bcR);


string_type_corr = 'Spearman';

% CorrCells:(with abs!!)

[rho_CorrCells, ~] = corr(M','Type',string_type_corr);
CorrCells = sum(sum(abs(triu(rho_CorrCells,1))))/(num_cells*(num_cells-1)/2);

% MC:

[rho_MC, ~] = corr(M,'Type',string_type_corr);
MC = sum(sum(abs(triu(rho_MC,1))))/(n*(n-1)/2);

end

