function [ M_diff_partweight ] = M_NoiseWeightAndB_function( num_cells, A, B, multi_start, T, sigma_weights, sigma_self)
    
n = size(A,2);

M_diff_partweight = zeros(num_cells,n);

for i = 1:num_cells
        
%     x0 = multi_start*rand(1,n);
    x0 = 0.5*ones(1,n);
         
    A_new = A;
    A_vec_noise_weights = 1+ sigma_weights*(2*rand(nnz(A_new),1)-1);
%     A_vec_noise_weights = normrnd(1,sigma_weights,[1,nnz(A_new)])';
    A_new(A_new>0) = A_new(A_new>0).*A_vec_noise_weights;
    A_new(sign(A_new)<0) = 0;
    
    B_new = B;
    B_vec_noise_weights = 1+ sigma_self*(2*rand(1,nnz(B_new))-1);
    B_new = B_new.*B_vec_noise_weights;
    B_new(sign(B_new)<0) = 0;
    
    [t,x] = ode45(@(t,x) hold_nodes_ABA2(t,x,A_new,B_new), [0,T],x0);
    M_diff_partweight(i,:) = x(length(t),:);
    
%     figure;
%     hold on
%     for j = 1:n
%         plot(t,x(:,j),'-o')
%     end
%     title('x')
%     hold off
    
end

% figure; plot(x);

M_diff_partweight = M_diff_partweight./sum(M_diff_partweight,2);

M_diff_partweight(M_diff_partweight<0)=0;

end

