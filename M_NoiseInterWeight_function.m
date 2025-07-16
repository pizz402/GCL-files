function [ M_diff_partweight ] = M_NoiseInterWeight_function( num_cells, A, B, T, p_noise, multi_weight, self_inter_cases)
n = size(A,2);

M_diff_partweight = zeros(num_cells,n);
index_template = find(A>0);
num_damage = floor(p_noise*length(index_template));

for i = 1:num_cells
    
%     clear x t
    
%     x0 = multi_start*rand(1,n);
    x0 = 0.5*ones(1,n);

    damege_links = randsample(index_template,num_damage);
        
    A_new = A;
%     A_new(damege_links) = rand(num_damage,1);
%     A_new(damege_links) = 1+rand(num_damage,1);

    
    switch self_inter_cases
        
        case 2 % only interactions 0 < aij < multi_weight
             A_new(damege_links) = multi_weight*rand(num_damage,1);
    
        case 3 % only interactions 1 < aij < 1 + multi_weight
            A_new(damege_links) = 1 + multi_weight*rand(num_damage,1);
        
    end
    
    [t,x] = ode45(@(t,x) hold_nodes_ABA2(t,x,A_new,B), [0,T],x0);
    M_diff_partweight(i,:) = x(length(t),:);

end

figure;plot(x)

M_diff_partweight = M_diff_partweight./sum(M_diff_partweight,2);

M_diff_partweight(M_diff_partweight<0)=0;

end

