function [ M_diff_starts ] = M_2eq_damagelinks_function( num_cells, percent_damage, n_hold, A_template, multi_start, T)

n = size(A_template,2);

index_template = find(A_template>0);
num_damage = round(percent_damage*length(index_template));

M_diff_starts = zeros(num_cells,n);

for i = 1:num_cells
    
    %%% A:
    A_new = A_template;
    
    % links_added:
    whole_vec = 1:n^2;
    whole_vec(index_template)=0;%unconnected links in the original network
    whole_vec(1:n+1:end)=0;
    whole_vec = whole_vec(whole_vec>0);
    links_added = randsample(whole_vec,num_damage);%adding links from the unconnected links
    
    % links_removed:
    links_removed = randsample(index_template,num_damage);
    
    A_new(links_removed)= 0;
    A_new(links_added) = rand(size(links_added));
    
    %%% ODE:
    
    z0 = multi_start*rand(1,2*n);
    
    permx = randperm(n);
    holdnodes = permx(1:n_hold);
    
    z0(holdnodes) = 0;

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
% 
end

M_diff_starts(M_diff_starts < 0) = 0;

M_diff_starts = M_diff_starts./sum(M_diff_starts,2);

end

