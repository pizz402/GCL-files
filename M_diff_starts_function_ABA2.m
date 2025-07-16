function [ M_nhold ] = M_diff_starts_function_ABA2( num_cells, A, B, A2, n_hold, multi_start, T )

n = A.n;

for i = 1:num_cells
    
    x0 = multi_start*rand(1,n);
%     x0(1:n_hold) = x0(1:n_hold)>rand();
    
    [t,x] = ode45(@(t,x) hold_nodes_ABA2(t,x,A.matrix,B.matrix,A2,n_hold), [0,T],x0);
    
    M_nhold(i,:) = x(length(t),:);
          
%     figure;
%     hold on
%     for j = (n_hold+1):n
% %     for j = 1:n
%         plot(t,x(:,j),'-o')
%     end
%     hold off

end

M_nhold(M_nhold < 0) = 0;

M_nhold = M_nhold./sum(M_nhold,2);

end

