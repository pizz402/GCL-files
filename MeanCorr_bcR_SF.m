function [ bcdVec, CMVec ] = MeanCorr_bcR_SF( num_cells, n, gamma_vec, bcd_itr, std_measurement, num_itr, norm_yesno, k0)
    

% k0 = 1;% arbitrary minimum 

bcdVec = zeros(length(gamma_vec),num_itr);
CMVec = zeros(length(gamma_vec),num_itr);

index_gamma_vec = 1;

for gamma = gamma_vec
    
    gamma
    
    for i = 1:num_itr
        
%         i
        
        % cell without noise
        one_cell = power_law_dist(k0,n,gamma);
%         one_cell = rand(1,n).^(-gamma);
%         [~,index] = sort(one_cell);

        % M matrix:
        % same cell - 
        M = repmat(one_cell,num_cells,1);
        M = M.*abs(normrnd(1,std_measurement,[num_cells,n]));  

        if (norm_yesno==1)
            
            % normalization - 
            M = M./sum(M,2);
            
        end

        
        [rho, ~] = corr(M, 'Type','Spearman');
        mean_corr = sum(sum(abs(triu(rho,1))))/(n*(n-1)/2);
        CMVec(index_gamma_vec,i) = mean_corr;
        
        bcd = new_bcdistcorr_itr(M, bcd_itr);
        bcdVec(index_gamma_vec,i) = bcd;
        
    end
          
    index_gamma_vec = index_gamma_vec+1;
    
%     figure;
%     histogram(mean(M))
%     title('M example')
%     colorbar;
% 
%     figure;
%     plot(sort(M));
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

