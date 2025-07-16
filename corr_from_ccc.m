function [CorrMatrix] = corr_from_ccc(num_cells, k_inter_vec, k_intra_vec, multi_weight, n, n_hold_zero, multi_start, T, sigma_weights)

string_type_corr = 'Spearman';

index_inter = 1;

for k_inter = k_inter_vec
    
    index_intra = 1;
    
    for k_intra = k_intra_vec
        
        [ A, B, A2 ] = Build_Network_activation_communities( k_inter, k_intra, n, multi_weight );
        M_diff_weight = M_diff_weight_function_ABA2( num_cells, A, B, A2, n_hold_zero, multi_start, T, sigma_weights);
        [corrAB,pvalAB] = corrAcorrB_updown(M_diff_weight, string_type_corr);
        CorrMatrix(index_inter,index_intra) = corrAB;
        
        index_intra = index_intra+1;
        
    end
    
    index_inter = index_inter+1;
    
end

%%
figure;
imagesc(k_intra_vec, k_inter_vec, CorrMatrix)
xlabel('INTRA')
ylabel('INTER')
title('CORR CCC')
colorbar;

%%
figure;
for index_intra = 1:length(k_intra_vec)
    plot(k_inter_vec,CorrMatrix(:,index_intra),'DisplayName',strcat('INTRA = ',num2str(k_intra_vec(index_intra))));
    hold on
end
xlabel('INTER')
title('CORR CCC')
legend show

%%
figure;
for index_inter = 1:length(k_inter_vec)
    plot(k_intra_vec,CorrMatrix(index_inter,:),'DisplayName',strcat('INTER = ',num2str(k_inter_vec(index_inter))));
    hold on
end
xlabel('INTRA')
title('CORR CCC')
legend show


end

