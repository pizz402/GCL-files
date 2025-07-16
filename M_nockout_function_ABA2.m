function [ M_nokout ] = M_nockout_function_ABA2( A, B, A2, multi_start, T )

n = A.n;
num_cells = n;

for i = 1:num_cells
    
    x0 = multi_start*rand(1,n);
    
    x0(i) = 0;
    
    [t,x] = ode45(@(t,x) nockout_nodes_ABA2(t,x,A.matrix,B.matrix,A2,i), [0,T],x0);
    
    M_nokout(i,:) = x(length(t),:);
          
%     figure;
%     hold on
%     for j = 1:n
%         plot(t,x(:,j),'-o')
%     end
%     hold off

end

M_nokout(M_nokout < 0) = 0;

M_nokout = M_nokout./sum(M_nokout,2);

end

