function [ M_diff_partweight ] = M_NoiseWeightAndSelf_function( num_cells, A, multi_start, T, sigma_weights, sigma_self)
    
n = size(A,2);

M_diff_partweight = zeros(num_cells,n);
A_diagonal = zeros(size(A));

for i = 1:num_cells
        
    x0 = multi_start*rand(1,n);
         
    A_new = A;
    A_new(1:n+1:end)=0;
    
    A_vec_noise_weights = 1+ sigma_weights*(2*rand(nnz(A_new),1)-1);
    A_new(A_new>0) = A_new(A_new>0).*A_vec_noise_weights;
%     A_new(sign(A_new)<0) = 0;

    A_diagonal(1:n+1:end) = A(1:n+1:end);
    A_vec_noise_self = 1+ sigma_self*(2*rand(nnz(A_diagonal),1)-1);
    A_diagonal(A_diagonal>0) = A_diagonal(A_diagonal>0).*A_vec_noise_self;
    
    A_new = A_new + A_diagonal;
    
    [t,x] = ode45(@(t,x) hold_nodes_ABA2(t,x,A_new), [0,T],x0);
    M_diff_partweight(i,:) = x(length(t),:);
    
%     figure;
%     hold on
%     for j = 1:n
%         plot(t,x(:,j),'-o')
%     end
%     title('x')
%     hold off
    
end

M_diff_partweight = M_diff_partweight./sum(M_diff_partweight,2);

M_diff_partweight(M_diff_partweight<0)=0;

end

