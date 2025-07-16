%% HEY ORR! GOOD LUCK! CTRL+ENTER EVERYTHING (:

format compact
clear
clc

num_cells = 200; % number of cells = size(M,1)%200
n = 400; % number of genes = size(M,2)%400

multi_weight = 1; % weights from the interval (0, multi_weight)
multi_start = 1; % initial conditions x0 from the interval (0, multi_start)
T = 40; % relaxation time

% for DOC
string_type_corr = 'Spearman';
overlap_type = 1;% Overlap with weights
threshold = 1e-3;% threshold for gene expression
rlowess_span = 0.1;% smooth (not used for now)

%% THE NETWORK

k_act = 2;
p_act = k_act/n;
A_template = multi_weight*sprand(n,n,p_act);% 0 < aij < 1
A_template(1:n+1:end)=0;% no loops

%% FIGS PARAMETERS 
% ORR YOU CAN CHANGE WHATEVER YOU WANT

percent_damage_vec = [0.005,0.05,0.5];% 3 values for 3 damages
num_itr_hist = 500; % the number of divitions for the COC histogram % set num_itr_hist = 100 for short runs

font_size = 14;% size font for figs

%% 6 FIGS FOR 3 DAMAGES (for num_itr_hist = 100 ~18 sec run, for num_itr_hist = 500 ~40 sec run)
% ORR YOU CAN CHANGE XLABEL & YLABLE & XTICKS

tic

num_figs = length(percent_damage_vec);%3 columns
num_rows = 2;%2 rows (can be 3)

figure;
axis square

min_hist = 1; max_hist = 0;
maxy_hist = 0;
for i = 1:num_figs
    
    subplot(num_rows,num_figs,i)
    p = percent_damage_vec(i);
    M_switching_links  = M_sameA_damagelinks_completelydiffweight_function( A_template, p, num_cells, multi_start, T);
    [corrAB, pvalAB] = corrAcorrB(M_switching_links, string_type_corr,1);%plot_YN = 1 for ploting 
    str_title = strcat('p = '," ",num2str(p));
    title(str_title)
    axis square
    xlabel('A correlation','Interpreter', 'LaTeX');%xlabel
    ylabel('B correlation','Interpreter', 'LaTeX');%ylabel
    xticks(-0.2:0.1:0.9)%xticks
    yticks(-0.2:0.1:0.9)%yticks
    set(gca,'fontsize', font_size);%fontsize
        
    subplot(num_rows,num_figs,num_figs+i)
    COC = corrAcorrB_hist(M_switching_links, string_type_corr,num_itr_hist);
    min_hist = min(min_hist,min(COC));
    max_hist = max(max_hist,max(COC));
    histogram(COC,'Normalization','probability');
    yl = get(gca,'ylim');
    maxy_hist = max(maxy_hist,yl(2));
    axis square
    xlabel('COC','Interpreter', 'LaTeX');%xlabel
    set(gca,'fontsize', font_size);%fontsize

end

% the same xlim and ylim for the histograms
for i = 1:num_figs
    subplot(num_rows,num_figs,num_figs+i)
    xlim([min_hist-0.05,max_hist+0.05])
    ylim([0,maxy_hist])
end

toc

%% ONE FIG - COC as function of damage p (for num_itr_hist = 100 ~10 minuts run)

tic

percent_damage_vec = 0.005:0.005:0.5;
% percent_damage_vec = [0.005,0.05,0.5];

mean_hist = zeros(size(percent_damage_vec));
err = zeros(size(percent_damage_vec));

for i = 1:length(percent_damage_vec)
    
    p = percent_damage_vec(i);
    M_switching_links  = M_sameA_damagelinks_completelydiffweight_function( A_template, p, num_cells, multi_start, T);
    COC = corrAcorrB_updown_hist(M_switching_links, string_type_corr,num_itr_hist);
    mean_hist(i) = mean(COC);
    err(i) = std(COC);
    
end

figure;
% subplot(num_rows,num_figs,7:9)%instead of "figure;" for num_rows=3
errorbar(percent_damage_vec,mean_hist,err)
xlabel('damage \textit{p}','Interpreter', 'LaTeX')
ylabel('COC','Interpreter', 'LaTeX');
set(gca,'fontsize', font_size);

toc

%% save to pdf & png

% str name
str = strcat( 'results genes/Fig for Orr/COC',num2str(round(1000*rand)),...
    '_NumCells',num2str(num_cells),'_NumGenes',num2str(n),'_NumItrHist',num2str(num_itr_hist));

% pdf
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf, 'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print(str,'-dpdf','-fillpage');

% png
saveas(gcf,strcat( str,'.png'));




