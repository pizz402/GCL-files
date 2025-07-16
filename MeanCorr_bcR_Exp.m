function [ bcdVec, CMVec ] = MeanCorr_bcR_Exp( num_cells, n, D_vec, bcd_itr, std_measurement, num_itr, norm_yesno)
    
bcdVec = zeros(length(D_vec),num_itr);
CMVec = zeros(length(D_vec),num_itr);

index_D_vec = 1;

for D = D_vec
    
    D
    
    for i = 1:num_itr
        
%         i
        
        % cell without noise
        one_cell = rand(1,n).^(-D);
%         [~,index] = sort(one_cell);

        % M matrix:
        % same cell - 
        M = repmat(one_cell,num_cells,1);
        M = M.*abs(normrnd(1,std_measurement.*M,[num_cells,n]));  

        if (norm_yesno==1)
            
            % normalization - 
            M = M./sum(M,2);
            
        end

        
        [rho, ~] = corr(M, 'Type','Spearman');
        mean_corr = sum(sum(abs(triu(rho,1))))/(n*(n-1)/2);
        CMVec(index_D_vec,i) = mean_corr;
        
        bcd = new_bcdistcorr_itr(M, bcd_itr);
        bcdVec(index_D_vec,i) = bcd;
        
    end
          
    index_D_vec = index_D_vec+1;
    
%     figure;
%     histogram(mean(M))
%     title('M example')
%     colorbar;

end

%     % cell without noise
%     one_cell = rand(1,n).^(-D);
%     [~,index] = sort(one_cell);
% 
%     % M matrix:
%     % same cell - 
%     M = repmat(one_cell,num_cells,1);
%     M = M.*abs(normrnd(1,std_measurement.*M,[num_cells,n]));
% 
%     % normalization - 
%     M = M./sum(M,2);
%     
%     % dcorr
%     bcR = new_bcdistcorr_itr(M, bcd_itr);
%     
%     %mean matrix corr:
%     [rho, ~] = corr(M(:,index),'Type','Spearman');
%     mean_corr = sum(sum(abs(triu(rho,1))))/(n*(n-1)/2);
    
   
end

