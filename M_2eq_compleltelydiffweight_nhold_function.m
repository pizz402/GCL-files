function [ M_diff_starts ] = M_2eq_compleltelydiffweight_nhold_function( num_cells, n_hold, A_template, multi_start, T)

n = size(A_template,2);

M_diff_starts = zeros(num_cells,n);

for i = 1:num_cells
    
    %%% A:
    A_new = A_template;
    A_new(A_new>0) = rand([1,nnz(A_new)])';
    A_new(sign(A_new)<0) = 0;
    
    %%% ODE:
    
    z0 = multi_start*rand(1,2*n);
    
    permx = randperm(n);
    holdnodes = permx(1:n_hold);
    
    z0(holdnodes) = 0;

    [t,z] = ode45(@(t,z) twoeq_mrna_protein(t,z,A_new,holdnodes), [0,T],z0);

    x = z(:,1:n);

    M_diff_starts(i,:) = x(length(t),:);
    
end

M_diff_starts(M_diff_starts < 0) = 0;

M_diff_starts = M_diff_starts./sum(M_diff_starts,2);

end

