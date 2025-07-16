function [corrAB,pvalAB,dmax] = corrAcorrB(M_matrix, string_type_corr,plot_YN)

n = size(M_matrix,2);% number of genes = size(M,2)
all_indexes = 1:n;
half_n = round(n/2);

indexesA = randsample(all_indexes,half_n);
all_indexes(indexesA)=0;
indexesB = find(all_indexes);

M_A = M_matrix(:,indexesA);
M_B = M_matrix(:,indexesB);

[ii,jj] = meshgrid(1:size(M_A,1),1:size(M_A,1));

[rho_A, pval_A] = corr(M_A','Type',string_type_corr);
triu_rho_vector_A = rho_A(jj>ii);

[rho_B, pval_B] = corr(M_B','Type',string_type_corr);
triu_rho_vector_B = rho_B(jj>ii);

[corrAB,pvalAB] = corr(triu_rho_vector_A,triu_rho_vector_B,'Type',string_type_corr);

if nargin==3
        
    [values, centers] = hist3([triu_rho_vector_A(:) triu_rho_vector_B(:)],[51 51]);
    dmax = max(values(:));
    imagesc(centers{:}, values.')
    set(gca,'YDir','normal')
    
end

end

