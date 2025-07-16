function [  ] = ploting_GCL_CM(x_vec,x_title,bcR,CM)

%% ploting GCL
ylimdown = 0;
figure;
errorbar(x_vec,mean(bcR'),std(bcR'),'LineWidth',2)
limsy=get(gca,'YLim');
ylimdown = min(ylimdown,limsy(1));
set(gca,'Ylim',[ylimdown limsy(2)]);
xlabel(x_title)
ylabel('Gene-to-Gene Coordination, GCL')
set(gca,'FontSize', 18)

%% ploting CM
ylimdown = 0;
figure;
errorbar(x_vec,mean(CM'),std(CM'),'LineWidth',2)
limsy=get(gca,'YLim');
ylimdown = min(ylimdown,limsy(1));
set(gca,'Ylim',[ylimdown limsy(2)]);
xlabel(x_title)
ylabel('Average of Pairwise Correlations')
set(gca,'FontSize', 18)

end

