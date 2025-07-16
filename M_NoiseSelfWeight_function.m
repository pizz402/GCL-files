function [ M_diff_partweight ] = M_NoiseSelfWeight_function( num_cells, A, B, T, p_noise, multi_weight, self_inter_cases)

n = size(A,2);

M_diff_partweight = zeros(num_cells,n);
num_damage = floor(p_noise*n);
diagonal = 1:n+1:n^2;

for i = 1:num_cells
    
%     x0 = multi_start*rand(1,n);
    x0 = 0.5*ones(1,n);
    
    damege_links = randsample(diagonal,num_damage);
    
    A_new = A;
%     A_new(damege_links) = rand(1,num_damage);
%     A_new(damege_links) =  1+rand(1,num_damage);

    switch self_inter_cases

        case 0 % only loops 0 < aii < multi_weight
            A_new(damege_links) = multi_weight*rand(1,num_damage);
      
        case 1 % only loops 1 < aii < 1 + multi_weight
            A_new(damege_links) = 1 + multi_weight*rand(1,num_damage);
            
    end

    [t,x] = ode45(@(t,x) hold_nodes_ABA2(t,x,A_new,B), [0,T],x0);
    M_diff_partweight(i,:) = x(length(t),:);

end

figure;plot(x)

M_diff_partweight = M_diff_partweight./sum(M_diff_partweight,2);

M_diff_partweight(M_diff_partweight<0)=0;

end

