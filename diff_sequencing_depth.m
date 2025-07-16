function [ M ] = diff_sequencing_depth( num_cells, A, B, n_hold, multi_start, T, depth_min, depth_max)

n = A.n;% number of genes

x0 = multi_start*rand(1,n);
x0(1:n_hold) = x0(1:n_hold)>rand();
[t,x] = ode45(@(t,x) hold_nodes(t,x,A.matrix,B.matrix,n_hold), [0,T],x0);

one_cell = x(length(t),:);
one_cell = one_cell/sum(one_cell);%normalize

M = zeros(num_cells,n);
depth_vec = linspace(depth_min,depth_max,num_cells);

for i = 1:num_cells
    
    depth_threshold = depth_vec(i);
    one_cell(one_cell < depth_threshold) = 0;
    M(i,:) = one_cell;
    
end

M = M./sum(M,2);

end

