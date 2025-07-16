function [ mean_corr, bcR, COC ] = MeanCorr_bcR_COC_Exp( num_cells, n, D, YNstring  )
    
    % cell without noise
    one_cell = rand(1,n).^(-D);
    [~,index] = sort(one_cell);

    % M matrix:
    % same cell - 
    M = repmat(one_cell,num_cells,1);
    M = M.*abs(normrnd(1,0.2.*M,[num_cells,n]));
    % normalization - 
    M = M./sum(M,2);
    
    % coc and dcorr
    COC = corrAcorrB_updown(M, 'Spearman');%COC
    bcR = bcdistcorr(M(:,1:n/2),M(:,(n/2+1):n));%dcorr
    
    %mean matrix corr:
%     [rho, ~] = corr(M,'Type','Spearman');
    [rho, ~] = corr(M(:,index),'Type','Spearman');
    mean_corr = sum(sum(abs(triu(rho,1))))/(n*(n-1)/2);
    
    if YNstring == 'Y'
        
        f = figure;
        suptitle(strcat('D=', num2str(D)))
        
        subplot(1,3,1)
        imagesc(rho);
%         title(strcat('D=', num2str(D)))
        colorbar
        
%         figure;
        subplot(1,3,2)
        plot(sort(one_cell));
        
%         figure;
        subplot(1,3,3)
        histogram(one_cell,10,'Normalization','probability');

        
    end

   
end

