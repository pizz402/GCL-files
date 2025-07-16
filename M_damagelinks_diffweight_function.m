function [ M_switching_links ] = M_damagelinks_diffweight_function( k_act, percent_damage, n, multi_weight, num_cells, multi_start, T, sigma_weights)
% switching links randomly

A.n = n;
A.k_act = k_act;
A.multi_weight = multi_weight;

p_act = k_act/n;
A_template = multi_weight*sprand(n,n,p_act);% 0 < aij < 1
A_template(1:n+1:end)=0;% no loops

index_template = find(A_template>0);
num_damage = round(percent_damage*length(index_template));

for i = 1:num_cells
    
    A_vec_noise_weights = (1-sigma_weights) + 2*sigma_weights*rand(nnz(A_template),1);
    
    A.matrix = A_template;
    A.matrix(A.matrix>0) = A.matrix(A.matrix>0).*A_vec_noise_weights;
    
    % nodes_added:
    whole_vec = 1:n^2;
    whole_vec(index_template)=[];
    nodes_added = randsample(whole_vec,num_damage);
    % nodes_removed:
    nodes_removed = randsample(index_template,num_damage);
    
    A.matrix(nodes_removed)= 0;
    A.matrix(nodes_added) = rand(size(nodes_added));
 
%     A.matrix = A_template;
    
    x0 = multi_start*rand(1,n);
     
    [t,x] = ode45(@(t,x) hold_nodes_ABA2(t,x,A.matrix), [0,T],x0);
    
    M_switching_links(i,:) = x(length(t),:);
    
end

M_switching_links(M_switching_links<0)=0;
M_switching_links = M_switching_links./sum(M_switching_links,2);

end

