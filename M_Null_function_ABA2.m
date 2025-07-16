function [ M_Null ] = M_Null_function_ABA2( num_cells, A, B, A2, n_hold, multi_start, T, sigma_results)

n = A.n;

x0 = multi_start*rand(1,n);
% x0(1:n_hold) = x0(1:n_hold)>rand();

[t,x] = ode45(@(t,x) hold_nodes_ABA2(t,x,A.matrix,B.matrix,A2,n_hold), [0,T],x0);

original_cell = x(length(t),:);

M_Null = repmat(original_cell,num_cells,1).*normrnd(1,sigma_results,[num_cells,length(original_cell)]);

M_Null(M_Null < 0) = 0;

M_Null = M_Null./sum(M_Null,2);


end

