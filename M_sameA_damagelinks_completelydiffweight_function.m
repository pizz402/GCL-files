function [ M_switching_links ] = M_sameA_damagelinks_completelydiffweight_function( A_template, percent_damage, num_cells, multi_start, T)
% switching links randomly

n = size(A_template,1);

index_template = find(A_template>0);
num_damage = round(percent_damage*length(index_template));

for i = 1:num_cells

    A.matrix = A_template;
    A.matrix(A.matrix>0) = rand([1,nnz(A.matrix)])';
    A.matrix(sign(A.matrix)<0) = 0;
    
    % nodes_added:
    whole_vec = 1:n^2;
    whole_vec(index_template)=[];%unconnected links in the original network
    links_added = randsample(whole_vec,num_damage);%adding links from the unconnected links
    % nodes_removed:
    links_removed = randsample(index_template,num_damage);
    
    A.matrix(links_removed)= 0;
    A.matrix(links_added) = rand(size(links_added));
 
%     A.matrix = A_template;
    
    x0 = multi_start*rand(1,n);
     
    [t,x] = ode45(@(t,x) hold_nodes_ABA2(t,x,A.matrix), [0,T],x0);
    
    M_switching_links(i,:) = x(length(t),:);
    
end

M_switching_links(M_switching_links<0)=0;
M_switching_links = M_switching_links./sum(M_switching_links,2);

end

