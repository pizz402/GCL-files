function [ M_diff_depth ] = M_diff_sequencing_depth_ABA2( num_cells, A, B, A2, n_hold, multi_start, T, depth_min, depth_max)

n = A.n;% number of genes

x0 = multi_start*rand(1,n);
% x0(1:n_hold) = x0(1:n_hold)>rand();

[t,x] = ode45(@(t,x) hold_nodes_ABA2(t,x,A.matrix,B.matrix,A2,n_hold), [0,T], x0);

one_cell = x(length(t),:);
one_cell = one_cell/sum(one_cell);%normalize

M_diff_depth = zeros(num_cells,n);
depth_vec = linspace(depth_min,depth_max,num_cells);

for i = 1:num_cells
    
    depth_threshold = depth_vec(i);
    one_cell(one_cell < depth_threshold) = 0;
    M_diff_depth(i,:) = one_cell;
    
end

M_diff_depth = M_diff_depth./sum(M_diff_depth,2);

end

