function [bcdVec, CMVec ] = bcd_CM_InterSelfNoise( num_cells_vec, A, B, T, p, num_itr, multi_weight, num_holds, bcd_itr)

n = size(A,2);

bcdVec = zeros(length(num_cells_vec),num_itr);
CMVec = zeros(length(num_cells_vec),num_itr);

index_num_cells_vec = 1;

for num_cells = num_cells_vec
    
    num_cells
    
    for i = 1:num_itr
        
%         i
        
        M = M_NoiseInterSelfWeight_function( num_cells, A, B, T, p, multi_weight, num_holds);

        [rho, ~] = corr(M,'Type','Spearman');
        mean_corr = sum(sum(abs(triu(rho,1))))/(n*(n-1)/2);
        CMVec(index_num_cells_vec,i) = mean_corr;
        
        bcd = new_bcdistcorr_itr(M, bcd_itr);
%         bcd = new_bcdistcorr(M);
        bcdVec(index_num_cells_vec,i) = bcd;
        
    end
          
    index_num_cells_vec = index_num_cells_vec+1;
    
%     figure;
%     histogram(mean(M))
%     title('M example')
%     colorbar;

end

% figure;
% imagesc(M)
% title('M example')
% colorbar;

% mean_bcd = mean(bcdVec');
% std_bcd = std(bcdVec');
% 
% mean_CM = mean(CMVec');
% std_CM = std(CMVec');
% 

% %%
% figure;
% errorbar(num_cells_vec,mean_bcd,std_bcd)
% xlabel('num cells')
% ylabel('GCL')
% % title('Inter & Self Noise')
% 
% %%
% figure;
% errorbar(num_cells_vec,mean_CM,std_CM)
% xlabel('num cells')
% ylabel('CM')
% % title('Inter & Self Noise')

end