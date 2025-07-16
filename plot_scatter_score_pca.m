function [  ] = plot_scatter_score_pca( M_1, M_2, M_3 )

figure;

NumOfCommunities = 3;
Colors = colormap(parula(NumOfCommunities));

M = [M_1; M_2; M_3];
[coeff,score] = pca(M);

indexes = [1, 1+size(M_1,1), 1+size(M_1,1)+size(M_2,1),size(M_1,1)+size(M_2,1)+size(M_3,1)];

for i = 1:NumOfCommunities

    scatter(score(indexes(i):indexes(i+1),1),score(indexes(i):indexes(i+1),2),'filled','MarkerEdgeColor',Colors(i,:));
    hold on     
    
end

title('PCA')
hold off

end