function [ M_diff_links ] = M_diff_links_function( k_act, percent_links_unfixed, n, multi_weight, num_cells, multi_start, T)
%% B & A2 are zeros:
k_inh = 0;
p_A2 = 0;
n_hold_zero = 0;
B.n = n;
B.k_inh = k_inh;
B.multi_weight = multi_weight;
p_inh = k_inh/n;
B.matrix = multi_weight*sprand(n,n,p_inh);% 0 < bij < 1
B.matrix(1:n+1:end)=0;% no loops
A2.p_A2 = p_A2;
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

%% from here - the parameters are relevant
A.n = n;
A.k_act = k_act;
A.multi_weight = multi_weight;

p_act_unfixed = percent_links_unfixed*k_act/n;
p_act_fixed = k_act/n - p_act_unfixed;% USE THIS!

A_fixed = multi_weight*sprand(n,n,p_act_fixed);% 0 < aij < 1
A_fixed(1:n+1:end)=0;% no loops

for i = 1:num_cells
    
    A_unfixed = multi_weight*sprand(n,n,p_act_unfixed);
    A_unfixed(1:n+1:end)=0; 
    
    A.matrix = A_fixed+A_unfixed;
    
    x0 = multi_start*rand(1,n);
     
    [t,x] = ode45(@(t,x) hold_nodes_ABA2(t,x,A.matrix,B.matrix,A2,n_hold_zero), [0,T],x0);
    
    M_diff_links(i,:) = x(length(t),:);
    
end

M_diff_links = M_diff_links./sum(M_diff_links,2);
M_diff_links(M_diff_links<0)=0;

end

