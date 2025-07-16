function [ bcR ] = new_bcdistcorr(M_matrix)

n = size(M_matrix,2);% number of genes = size(M,2)
all_indexes = 1:n;
half_n = round(n/2);

indexesA = randsample(all_indexes,half_n);
all_indexes(indexesA)=0;
indexesB = all_indexes>0;

bcR = bcdistcorr( M_matrix(:,indexesA),M_matrix(:,indexesB));


end

