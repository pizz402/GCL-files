function [ M_diff_partweight ] = M_NoiseInterSelfWeight_function( num_cells, A, Bvec, T, p_noise, multi_weight, num_holds)%, self_inter_cases)

n = size(A,2);

M_diff_partweight = zeros(num_cells,n);
diagonal = 1:n+1:n^2;

A_no_diagonal = A;
A_no_diagonal(diagonal) = 0;
index_template = find(A_no_diagonal>0);

num_damage_1 = floor(p_noise*length(index_template));%no diagonal?
num_damage_2 = floor(p_noise*n);%diagonal?

for i = 1:num_cells
    
%     x0 = multi_start*rand(1,n);
    x0 = 0.5*ones(1,n);
    
    n_hold_vec = randsample(1:n,num_holds);
    x0(n_hold_vec) = 0;
    
    damege_links_1 = randsample(index_template,num_damage_1);
    damege_links_2 = randsample(diagonal,num_damage_2);
    
    A_new = A;
            
    A_new(damege_links_1) = multi_weight*rand(num_damage_1,1);
    A_new(damege_links_2) = multi_weight*rand(1,num_damage_2);
            
    [t,x] = ode45(@(t,x) hold_nodes_ABA2(t,x,A_new,Bvec,n_hold_vec), [0,T],x0);
    M_diff_partweight(i,:) = x(length(t),:);

end

% figure;plot(x)

M_diff_partweight = M_diff_partweight./sum(M_diff_partweight,2);

M_diff_partweight(M_diff_partweight<0)=0;

end

