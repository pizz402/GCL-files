function [ M_diff_starts ] = M_2eq_function_weightshold( num_cells, n_hold, sigma_weights, A_template, multi_start, T)

n = size(A_template,2);

M_diff_starts = zeros(num_cells,n);
 
for i = 1:num_cells
    
    z0 = multi_start*rand(1,2*n);
    
    % nhold noise
    permx = randperm(n);
    holdnodes = permx(1:n_hold);
    z0(holdnodes) = 0;
    
    % sigma weights noise on all weights
    A_new = A_template;
    A_vec_noise_weights = 1 + sigma_weights*(2*rand(nnz(A_new),1)-1);
    A_new(A_template>0) = A_new(A_new>0).*A_vec_noise_weights;

    [t,z] = ode45(@(t,z) twoeq_mrna_protein(t,z,A_new,holdnodes), [0,T],z0);

    x = z(:,1:n);

    M_diff_starts(i,:) = x(length(t),:);
          
%     figure;
%     hold on
%     for j = 1:n
%         plot(t,x(:,j),'-o')
%     end
%     title('x')
%     hold off

%     y = z(:,n+1:2*n);
%     figure;
%     hold on
%     for j = 1:n
%         plot(t,y(:,j),'-o')
%     end
%     title('y')
%     hold off

end

M_diff_starts(M_diff_starts < 0) = 0;

M_diff_starts = M_diff_starts./sum(M_diff_starts,2);

end

