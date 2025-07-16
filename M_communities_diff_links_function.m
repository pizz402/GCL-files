function [ M_communities_diff_links ] = M_communities_diff_links_function( k_act, sigma_weights, percent_links_unfixed_betweenComm, n, multi_weight, num_cells, multi_start, T)
%% B & A2 are zeros:
% k_inh = 0;
% p_A2 = 0;
% n_hold_zero = 0;
% B.n = n;
% B.k_inh = k_inh;
% B.multi_weight = multi_weight;
% p_inh = k_inh/n;
% B.matrix = multi_weight*sprand(n,n,p_inh);% 0 < bij < 1
% B.matrix(1:n+1:end)=0;% no loops
% A2.p_A2 = p_A2;
% A2.genes_with_activationDim2 = rand(n,1)<p_A2;
% A2.firstvec = randperm(n)';
% A2.secondvec = randperm(n)';
% indexes_of_genes_with_activationDim2 = find(A2.genes_with_activationDim2);
% for i = 1 : length(indexes_of_genes_with_activationDim2)
%     gene = indexes_of_genes_with_activationDim2(i);
%     if (A2.firstvec(gene) == gene) || (A2.secondvec(gene) == gene)
%         A2.genes_with_activationDim2(gene) = 0;
%     end
% end

%% two communities with the same size (so 2 must devide n)

if ~~mod(n,2)
    'choose other n'
    return
end

half_n = floor(n/2);

%% IM HERE

A.n = n;
A.k_act = k_act;
A.multi_weight = multi_weight;

% 1/2 degree InComm and 1/2 BetweenComm:
p_InComm = k_act/n;
p_BetweenComm = k_act/n;

percent_links_fixed_betweenComm = 1 - percent_links_unfixed_betweenComm;

A_fixed_InComm = multi_weight*blkdiag(sprand(half_n, half_n, p_InComm), sprand(half_n, half_n, p_InComm));
A_fixed_BeteenComm = multi_weight*fliplr(blkdiag(sprand(half_n, half_n, percent_links_fixed_betweenComm*p_BetweenComm),...
    sprand(half_n, half_n, percent_links_fixed_betweenComm*p_BetweenComm)));

A_fixed = A_fixed_InComm + A_fixed_BeteenComm;

for i = 1:num_cells
    
    A_unfixed = multi_weight*fliplr(blkdiag(sprand(half_n, half_n, percent_links_unfixed_betweenComm*p_BetweenComm),...
    sprand(half_n, half_n, percent_links_unfixed_betweenComm*p_BetweenComm)));
    
    A.matrix = A_fixed + A_unfixed;
    A.matrix(1:n+1:end)=0;% no loops
    
    A_new = A.matrix;
%     A_new(A_new>0) = A_new(A_new>0).*normrnd(1,sigma_weights,[1,nnz(A_new)])';
    A_vec_noise_weights = (1-sigma_weights) + 2*sigma_weights*rand(nnz(A_new),1);%%%%%%%%%%%%% IM HERE!
    A_new(A_new>0) = A_new(A_new>0).*A_vec_noise_weights;
    A_new(sign(A_new)<0) = 0;
    
    x0 = multi_start*rand(1,n);
    
%     [t,x] = ode45(@(t,x) hold_nodes_ABA2(t,x,A_new,B.matrix,A2,n_hold_zero), [0,T],x0);
    [t,x] = ode45(@(t,x) hold_nodes_ABA2(t,x,A_new), [0,T],x0);
    
    M_communities_diff_links(i,:) = x(length(t),:);
    
end

M_communities_diff_links = M_communities_diff_links./sum(M_communities_diff_links,2);
M_communities_diff_links(M_communities_diff_links<0)=0;

end

