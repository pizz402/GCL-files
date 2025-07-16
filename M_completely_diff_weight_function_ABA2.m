function [ M_diff_weight ] = M_completely_diff_weight_function_ABA2( num_cells, A, B, A2, n_hold, multi_start, T)
% weight = weight*(1 + epsilon)

n = A.n;

for i = 1:num_cells
    
    x0 = multi_start*rand(1,n);
%     x0(1:n_hold) = x0(1:n_hold)>rand();
       
    A_new = A.matrix;
    A_new(A_new>0) = rand([1,nnz(A_new)])';
    A_new(sign(A_new)<0) = 0;
    
    B_new = B.matrix;
    B_new(B_new>0) = rand([1,nnz(B_new)])';
    B_new(sign(B_new)<0) = 0;
    
    [t,x] = ode45(@(t,x) hold_nodes_ABA2(t,x,A_new,B_new,A2,n_hold), [0,T],x0);
    
    M_diff_weight(i,:) = x(length(t),:);
    
end

M_diff_weight = M_diff_weight./sum(M_diff_weight,2);

M_diff_weight(M_diff_weight<0)=0;

end

