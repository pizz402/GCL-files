function [ M_diff_partweight ] = M_diff_diffpartweight_function_ABA2( num_cells, A, B, A2, n_hold, multi_start, T, sigma_weights, fraction_noise)
 % 0 < sigma_weights < 1 - uniform noise
fraction_samesamples = 1-fraction_noise;%fraction of wieghts that doesn't change

A_num_samesamples = round(length(A.matrix(A.matrix>0))*fraction_samesamples);
B_num_samesamples = round(length(B.matrix(B.matrix>0))*fraction_samesamples);

n = A.n;

for i = 1:num_cells
        
    A_rand_indexes_samesamples_unsorted = randsample(1:length(A.matrix(A.matrix>0)),A_num_samesamples);
    B_rand_indexes_samesamples_unsorted = randsample(1:length(B.matrix(B.matrix>0)),B_num_samesamples);

    x0 = multi_start*rand(1,n);
     
%      weight = weight*(1 + epsilon):
    
    A_new = A.matrix;
    A_vec_noise_weights = (1-sigma_weights) + 2*sigma_weights*rand(nnz(A_new),1);%%%%%%%%%%%%% IM HERE!
%     A_vec_noise_weights = normrnd(1,sigma_weights,[1,nnz(A_new)])';
    A_vec_noise_weights(A_rand_indexes_samesamples_unsorted) = ones(length(A_rand_indexes_samesamples_unsorted),1);
    A_new(A_new>0) = A_new(A_new>0).*A_vec_noise_weights;
    A_new(sign(A_new)<0) = 0;
    
    B_new = B.matrix;
    B_vec_noise_weights = (1-sigma_weights) + 2*sigma_weights*rand(nnz(B_new),1);%%%%%%%%%%%%% IM HERE!
%     B_vec_noise_weights = normrnd(1,sigma_weights,[1,nnz(B_new)])';
    B_vec_noise_weights(B_rand_indexes_samesamples_unsorted) = ones(length(B_rand_indexes_samesamples_unsorted),1);
    B_new(B_new>0) = B_new(B_new>0).*B_vec_noise_weights;
    B_new(sign(B_new)<0) = 0;

    [t,x] = ode45(@(t,x) hold_nodes_ABA2(t,x,A_new,B_new,A2,n_hold), [0,T],x0);
    M_diff_partweight(i,:) = x(length(t),:);
    
    
end

M_diff_partweight = M_diff_partweight./sum(M_diff_partweight,2);

M_diff_partweight(M_diff_partweight<0)=0;

end

