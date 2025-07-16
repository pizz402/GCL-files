function [ DOC_curve ] = DOC( M_matrix, threshold, span, overlap_type)%, Dis_type)

DOC_curve.span = span;
DOC_curve.threshold = threshold;

M_relevant = (M_matrix>threshold).*M_matrix;% values under threshold = 0
c = combnk(1:size(M_relevant,1),2);% all combinations of cells (i,j)

switch overlap_type
    
    case 1 % Overlap
        
        overlap_matrix = (M_relevant(c(:,1),:) & M_relevant(c(:,2),:)).*(M_relevant(c(:,1),:) + M_relevant(c(:,2),:))./2;
        overlap_vector = sum(overlap_matrix')';
    
    case 2 % Jaccard (no weights)
        
        overlap_matrix = (M_relevant(c(:,1),:) & M_relevant(c(:,2),:));
        overlap_vector = sum(overlap_matrix')'./sum((M_relevant(c(:,1),:) | M_relevant(c(:,2),:))')';
    
    otherwise % new
        
        overlap_matrix = not(xor(M_relevant(c(:,1),:) , M_relevant(c(:,2),:)));
        overlap_vector = sum(overlap_matrix')'./size(M_relevant,2);
              
end
        
if isempty(find(overlap_vector>0))

    'no overlap vector!'

end

disimilarity_vector = 2*ones(length(overlap_vector),1);

for relevant_i = find(overlap_vector>0)'

    v1 = M_relevant(c(relevant_i,1),find(overlap_matrix(relevant_i,:)>0));
    v1 = v1/sum(v1);
    v2 = M_relevant(c(relevant_i,2),find(overlap_matrix(relevant_i,:)>0));
    v2 = v2/sum(v2); 
    
    unzero_idexes = v1&v2;
    v1 = v1(unzero_idexes);
    v2 = v2(unzero_idexes);

    [rho(relevant_i),pval_disimilarity(relevant_i)] = corr(v1',v2','Type','Spearman');

    disimilarity_vector(relevant_i) = 1-rho(relevant_i);

end

disimilarity_vector(isnan(disimilarity_vector)) = 0;

new_disimilarity_vector = disimilarity_vector(disimilarity_vector<2);
new_overlap_vector = overlap_vector(disimilarity_vector<2);
plot(new_overlap_vector,new_disimilarity_vector,'.')
hold on 

[xx,ind] = sort(new_overlap_vector);
yy = smooth(new_disimilarity_vector(ind),span,'rlowess');
plot(xx,yy,'r','LineWidth',2)
hold on

corr_DOC = corr(xx, new_disimilarity_vector(ind));

[coefficients,~] = polyfit(xx, new_disimilarity_vector(ind), 1);
slope = coefficients(1);

DOC_curve.slope = slope;
DOC_curve.corr = corr_DOC;

xlabel('Overlap')
ylabel('Disimilarity')

end

