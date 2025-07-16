function [ A, B, A2 ] = Build_Network_activation_communities( k_inter, k_intra, n, multi_weight )

% two communities with the same size (so 2 must devide n)
if ~~mod(n,2)
    'choose other n'
    return
end

half_n = floor(n/2);

A.n = n;
A.k_act = k_intra+k_inter;
A.multi_weight = multi_weight;

p_intra = k_intra/half_n;
p_inter = k_inter/half_n;

A.matrix = multi_weight*blkdiag(sprand(half_n, half_n, p_intra), sprand(half_n, half_n, p_intra))...
    + multi_weight*fliplr(blkdiag(sprand(half_n, half_n, p_inter), sprand(half_n, half_n, p_inter)));

A.matrix(1:n+1:end)=0;% no loops


% from here - for the sake of propriety

k_inh = 0;
p_inh = k_inh/n;%0
B.matrix = sprand(n,n,p_inh);%zeros

A2.genes_with_activationDim2 = false(n,1);
A2.firstvec = ones(n,1);
A2.secondvec = ones(n,1);

end
