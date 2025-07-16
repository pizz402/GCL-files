function [ A, B, A2 ] = Build_Network_activation_3communities( k_inter, k_intra,k_intra_C, n, multi_weight )

% three communities with the same size (so 3 must devide n)
if ~~mod(n,3)
    'choose other n'
    return
end

n_third = floor(n/3);

A.n = n;
A.k_act = k_intra+k_inter;
A.multi_weight = multi_weight;

p_intra = k_intra/n_third;
p_intra_C = k_intra_C/n_third;
p_inter = k_inter/n_third;

A.matrix = multi_weight*[sprand(n_third, n_third, p_intra), sprand(n_third, n_third, p_inter), sprand(n_third, n_third, 0);...
    sprand(n_third, n_third, p_inter), sprand(n_third, n_third, p_intra), sprand(n_third, n_third, 0);...
    sprand(n_third, n_third, 0), sprand(n_third, n_third, 0), sprand(n_third, n_third, p_intra_C)];

A.matrix(1:n+1:end)=0;% no loops


% from here - for the sake of propriety

k_inh = 0;
p_inh = k_inh/n;%0
B.matrix = sprand(n,n,p_inh);%zeros

A2.genes_with_activationDim2 = false(n,1);
A2.firstvec = ones(n,1);
A2.secondvec = ones(n,1);

end
