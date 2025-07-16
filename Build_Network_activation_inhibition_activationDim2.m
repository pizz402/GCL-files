function [ A, B, A2 ] = Build_Network_activation_inhibition_activationDim2( k_act, k_inh, p_A2, n, multi_weight )

A.n = n;
A.k_act = k_act;
A.multi_weight = multi_weight;

B.n = n;
B.k_inh = k_inh;
B.multi_weight = multi_weight;

A2.p_A2 = p_A2;

p_act = k_act/n;
p_inh = k_inh/n;

A.matrix = multi_weight*sprand(n,n,p_act);% 0 < aij < 1
A.matrix(1:n+1:end)=0;% no loops

B.matrix = multi_weight*sprand(n,n,p_inh);% 0 < bij < 1
B.matrix(1:n+1:end)=0;% no loops


A2.genes_with_activationDim2 = rand(n,1)<p_A2;

A2.firstvec = randperm(n)';
A2.secondvec = randperm(n)';
indexes_of_genes_with_activationDim2 = find(A2.genes_with_activationDim2);

for i = 1 : length(indexes_of_genes_with_activationDim2)
    
    gene = indexes_of_genes_with_activationDim2(i);
    
    if (A2.firstvec(gene) == gene) || (A2.secondvec(gene) == gene)
        
        A2.genes_with_activationDim2(gene) = 0;

    end
    
end

end
