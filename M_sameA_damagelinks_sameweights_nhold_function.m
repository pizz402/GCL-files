function [ M_switching_links ] = M_sameA_damagelinks_sameweights_nhold_function( A_template, percent_damage, n_hold, num_cells, multi_start, T)
% switching links randomly

n = size(A_template,1);

index_template = find(A_template>0);
num_damage = round(percent_damage*length(index_template));

M_switching_links =  zeros(num_cells,n);

for i = 1:num_cells

    A.matrix = A_template;
%     A.matrix(A.matrix>0) = rand([1,nnz(A.matrix)])';
%     A.matrix(sign(A.matrix)<0) = 0;
    
    % links_added:
    whole_vec = 1:n^2;
    whole_vec(index_template)=0;%unconnected links in the original network
    whole_vec(1:n+1:end)=0;
    whole_vec = whole_vec(whole_vec>0);
    links_added = randsample(whole_vec,num_damage);%adding links from the unconnected links
    
    % links_removed:
    links_removed = randsample(index_template,num_damage);
    
    A.matrix(links_removed)= 0;
    A.matrix(links_added) = rand(size(links_added));
    
    x0 = multi_start*rand(1,n);
    
    permx = randperm(n);
    holdnodes = permx(1:n_hold);
    
    x0(holdnodes) = 0;
     
    [t,x] = ode45(@(t,x) hold_nodes_A(t,x,A.matrix,holdnodes), [0,T],x0);
    
    M_switching_links(i,:) = x(length(t),:);
    
end

M_switching_links(M_switching_links<0)=0;
M_switching_links = M_switching_links./sum(M_switching_links,2);

end

