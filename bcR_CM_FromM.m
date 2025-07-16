function [bcR, CM] = bcR_CM_FromM(M)

string_type_corr = 'Spearman';
n = size(M,2);

bcR = new_bcdistcorr(M);

% CM:
[rho_CM, ~] = corr(M,'Type',string_type_corr);
CM = sum(sum(abs(triu(rho_CM,1))))/(n*(n-1)/2);

end

