function [ M_diff_starts ] = M_1eq_function_diffweightshold( num_cells, n_hold, p_weights, A_template, B_template, A2_template, A2_twovec, multi_start, multi_weight, T, c_hill)

n = size(A_template,2);

M_diff_starts = zeros(num_cells,n);

num_damage_A = round(p_weights*nnz(A_template));
index_template_A = find(A_template>0);

num_damage_B = round(p_weights*nnz(B_template));
index_template_B = find(B_template>0);

num_damage_A2 = round(p_weights*nnz(A2_template));
index_template_A2 = find(A2_template>0);

for i = 1:num_cells
    
    x0 = multi_start*rand(1,n);
%     x0 = 0.5*ones(1,n);
    
    % nhold noise
    holdnodes = randsample(1:n,n_hold);
    x0(holdnodes) = 0;
    
    % changing the values of p_weights of the links
    A_new = A_template;
    damege_links_A = randsample(index_template_A,num_damage_A);
    A_new(damege_links_A) = multi_weight.*rand([1,num_damage_A])';
    
    % changing the values of p_weights of the links
    B_new = B_template;
    damege_links_B = randsample(index_template_B,num_damage_B);
    B_new(damege_links_B) = multi_weight.*rand([1,num_damage_B])';
    
    % changing the values of p_weights of the links
    A2_new = A2_template;
    damege_links_A2 = randsample(index_template_A2,num_damage_A2);
    A2_new(damege_links_A2) = multi_weight.*rand([1,num_damage_A2])';
     
    [t,x] = ode45(@(t,x) hold_nodes_A_B_A2(t,x,A_new,B_new, A2_new, A2_twovec, holdnodes, c_hill), [0,T],x0);
%     [t,x] = ode45(@(t,x) hold_nodes_A(t,x,A_new,holdnodes), [0,T],x0);%     old version - no B_template

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

