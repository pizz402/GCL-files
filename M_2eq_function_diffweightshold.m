function [ M_diff_starts ] = M_2eq_function_diffweightshold( num_cells, n_hold, p_weights, A_template, multi_start, multi_weight, T)

% multi_weight = 1;

n = size(A_template,2);

M_diff_starts = zeros(num_cells,n);

num_damage = round(p_weights*nnz(A_template));
index_template = find(A_template>0);

for i = 1:num_cells
    
    z0 = multi_start*rand(1,2*n);
    
    % nhold noise
    permx = randperm(n);
    holdnodes = permx(1:n_hold);
    z0(holdnodes) = 0;
    
    % changing the values of p_weights of the links
    A_new = A_template;
    damege_links = randsample(index_template,num_damage);
    A_new(damege_links) = multi_weight.*rand([1,num_damage])';

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

