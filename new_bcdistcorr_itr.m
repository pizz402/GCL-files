function [ bcR ] = new_bcdistcorr_itr(M_matrix,bcd_itr)

bcRvec = zeros(1,bcd_itr);

for i = 1: bcd_itr
    i
    bcRvec(i) =  new_bcdistcorr(M_matrix);
    
end

bcR = mean(bcRvec);

end

