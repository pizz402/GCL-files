function [ M_diff_links ] = M_difflinks_diffweight_function( k_act, percent_links_unfixed, n, multi_weight, num_cells, multi_start, T, sigma_weights)
% same core and adding diff links for every cell

A.n = n;
A.k_act = k_act;
A.multi_weight = multi_weight;

p_act_unfixed = percent_links_unfixed*k_act/n;
p_act_fixed = k_act/n - p_act_unfixed;

A_fixed = multi_weight*sprand(n,n,p_act_fixed);% 0 < aij < 1
A_fixed(1:n+1:end)=0;% no loops

for i = 1:num_cells
    
    A_vec_noise_weights = (1-sigma_weights) + 2*sigma_weights*rand(nnz(A_fixed),1);%%%%%%%%%%%%% IM HERE!
    A_fixed(A_fixed>0) = A_fixed(A_fixed>0).*A_vec_noise_weights;
    A_fixed(sign(A_fixed)<0) = 0;
    
    A_unfixed = multi_weight*sprand(n,n,p_act_unfixed);
    A_unfixed(1:n+1:end)=0; 
    
    A.matrix = A_fixed+A_unfixed;
    
    x0 = multi_start*rand(1,n);
     
    [t,x] = ode45(@(t,x) hold_nodes_ABA2(t,x,A.matrix), [0,T],x0);
    
    M_diff_links(i,:) = x(length(t),:);
    
end

M_diff_links(M_diff_links<0)=0;
M_diff_links = M_diff_links./sum(M_diff_links,2);

end

