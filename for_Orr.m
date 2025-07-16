%% for_Orr
% ctrl+enter on the fixed parameters from pdf4

%% 
% M_difflinks_diffweight_function - same core and adding links randomly. all in one FIG.
for same_core_and_adding_links_randomly = 1:1
    
k_act = 2; 
sigma_weights = 0.5;%uniform dist
percent_links_unfixed_vec = [0.1,0.4,0.9];
string_type_corr = 'Spearman';

figure;
labels = {};
for percent_links_unfixed = percent_links_unfixed_vec

    M_diff_links = M_difflinks_diffweight_function( k_act, percent_links_unfixed, n, multi_weight, num_cells, multi_start, T, sigma_weights);

    [corrAB, pvalAB] = corrAcorrB(M_diff_links, rlowess_span, string_type_corr);
    hold on
    
    str_label = strcat('damage = ', num2str(percent_links_unfixed),', corr = ',num2str(corrAB),', pval = ', num2str(pvalAB));
    labels{end+1} = str_label;
    
end
legend(labels, 'Location', 'NorthWest')
title('COC')

end

%%
%% 
% M_damagelinks_diffweight_function - switching links randomly. all in one FIG.
for switching_links_randomly = 1:1
    
k_act = 2; 
sigma_weights = 0.2;%uniform dist
percent_damage_vec = [0.1,0.4,0.7,1];
string_type_corr = 'Spearman';

figure;
labels = {};
for percent_damage = percent_damage_vec

    M_switching_links = M_damagelinks_diffweight_function( k_act, percent_damage, n, multi_weight, num_cells, multi_start, T, sigma_weights);

    [corrAB, pvalAB] = corrAcorrB(M_switching_links, rlowess_span, string_type_corr);
    hold on
    
    str_label = strcat('damage = ', num2str(percent_damage),', corr = ',num2str(corrAB),', pval = ', num2str(pvalAB));
    labels{end+1} = str_label;
    
end
legend(labels, 'Location', 'NorthWest')
title('COC')

end

%%
%% 3 FIGs for different damage of switching - DRAFT FIG 1

k_act = 2; 
percent_damage_vec = [0.005,0.1,0.5];
num_itr_hist = 400;
num_fig = 3;

string_type_corr = 'Spearman';
rlowess_span = 0.05;% smooth

figure;
min_hist = 1; max_hist = 0;
for i = 1:length(percent_damage_vec)%1:4
    
    subplot(2,num_fig,i)
    p = percent_damage_vec(i);
    M_switching_links = M_damagelinks_completelydiffweight_function( k_act, p, n, multi_weight, num_cells, multi_start, T);
    [corrAB, pvalAB] = corrAcorrB(M_switching_links, rlowess_span, string_type_corr,1);%plot_YN = 1 for ploting 
%     limit1 = min([xlim,ylim]); limit2 = max([xlim,ylim]);
%     axis([limit1 limit2 limit1 limit2])
    axis square
    str_title = num2str(p);
    title(str_title)
    box on
    set(gca,'fontsize', font_size);
        
    subplot(2,num_fig,num_fig+i)
    COC = corrAcorrB_hist(M_switching_links, rlowess_span, string_type_corr,num_itr_hist);
%     mean_hist = mean(COC);
    min_hist = min(min_hist,min(COC));
    max_hist = max(max_hist,max(COC));
%     histogram(COC,'Normalization','probability');
    histogram(COC,20,'Normalization','probability');
    axis square
    hold on
%     quiver([mean_hist,0],[mean_hist,0.1])
%     annotation('textarrow',[mean_hist,0],[mean_hist,0.1],'Color', 'r');
%     line([mean_hist, mean_hist], ylim, 'LineWidth', 2, 'Color', 'r','LineStyle','--');
    
end

for i = 1:num_fig
    subplot(2,num_fig,num_fig+i)
    xlim([min_hist-0.05,max_hist+0.05])
end

% put_letters_on_corner_of_subfigures(gcf)%for 4*2 Figs

%%
%% COC as function of damage p - DRAFT FIG 2
tic

string_type_corr = 'Spearman';
rlowess_span = 0.05;% smooth

k_act = 2;
num_itr_hist = 400;

percent_damage_vec = 0.005:0.005:0.5;
% percent_damage_vec = [0.005,0.05,0.5];
mean_hist = zeros(size(percent_damage_vec));
err = zeros(size(percent_damage_vec));

for i = 1:length(percent_damage_vec)%1:4
    
    p = percent_damage_vec(i);
    M_switching_links = M_damagelinks_completelydiffweight_function( k_act, p, n, multi_weight, num_cells, multi_start, T);
    COC = corrAcorrB_hist(M_switching_links, rlowess_span, string_type_corr,num_itr_hist);
    mean_hist(i) = mean(COC);
    err(i) = std(COC);
    
end

figure;
errorbar(percent_damage_vec,mean_hist,err)
xlabel('damage \textit{p}','Interpreter', 'LaTeX')
ylabel('COC')

toc

%%
%% Phase diagram
for i = 1:1
k_act = 2; 
sigma_weights_vec = 0.1:0.1:1;
percent_damage_vec = 0.01:0.01:0.2;
num_itr_hist = 100;
string_type_corr = 'Spearman';
COC_matrix = zeros(length(sigma_weights_vec),length(percent_damage_vec));

for s_i = 1:length(sigma_weights_vec)
    
    sigma_weights = sigma_weights_vec(s_i);
    
    for p_i = 1:length(percent_damage_vec)

        p = percent_damage_vec(p_i);      
        M_switching_links = M_damagelinks_diffweight_function( k_act, p, n, multi_weight, num_cells, multi_start, T, sigma_weights);
        COC = corrAcorrB_hist(M_switching_links, rlowess_span, string_type_corr,num_itr_hist);        
        mean_hist = mean(COC);        
        COC_matrix(p_i,s_i) = mean_hist;
        
    end
    
end

%%

figure;
imagesc(sigma_weights_vec, percent_damage_vec, COC_matrix)
set(gca,'YDir','normal')
xlabel('Sigma Weights Noise')
ylabel('Percent Damage Noise')
title('CORR CCC')
colorbar;

%%
% contour_corr_cell = 0.5;
hold on
% contour(sigma_weights_vec, percent_damage_vec, COC_matrix',[contour_corr_cell,contour_corr_cell],'r','LineWidth',3);
contour(sigma_weights_vec, percent_damage_vec, COC_matrix','r','LineWidth',3);

end


%% before saveing
font_size = 24;
set(gca,'fontsize', font_size);
%% before saving -2
subplot(2,3,1)
xticks(0.6:0.2:0.8)
%% before saving -3
font_size = 16;
for i = 1:num_fig*2
    subplot(2,num_fig,i)
    set(gca,'fontsize', font_size);
end

%% save to pdf & png

str = strcat( 'results genes/Fig for Orr/COC',num2str(round(1000*rand)),...
    '_NumCells',num2str(num_cells),'_NumGenes',num2str(n),'_NumItrHist',num2str(num_itr_hist));

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf, 'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);

print(str,'-dpdf','-fillpage');
saveas(gcf,strcat( str,'.png'));
