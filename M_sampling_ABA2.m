function [ M_sample ] = M_sampling_ABA2( num_cells, A, B, A2, n_hold, multi_start, T, num_sampl)

n = A.n;% number of genes

x0 = multi_start*rand(1,n);
% x0(1:n_hold) = x0(1:n_hold)>rand();

[t,x] = ode45(@(t,x) hold_nodes_ABA2(t,x,A.matrix,B.matrix,A2,n_hold), [0,T],x0);

one_cell = x(length(t),:);
one_cell = one_cell/sum(one_cell);

M_sample = zeros(num_cells,n);

for i = 1:num_cells
    
    genes_sampl(i,:) = discretesample(one_cell,num_sampl);
    
    for j = 1:n
        M_sample(i,j) = sum(genes_sampl(i,:)==j);
    end
 
end

M_sample = M_sample./sum(M_sample,2);

end

