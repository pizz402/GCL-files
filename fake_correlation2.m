%% fake correlation2 (without gene-regulation) 
% (1) compositional data & (2) data from corr matrix
format compact
clear
clc

%%
num_cells = 100; % number of cells = size(M,1)
n = 200; % number of genes = size(M,2)

%%
%% (1) - compositional data

%% all in one plot
tic
num_iterations = 20;%20
D_vec = 0.1:0.05:1;
sp_D_vec = [0.1,0.5,1];
plotall_compositionaldata(num_cells, n, D_vec, num_iterations, sp_D_vec)
toc
fontsize = 16; linewidth = 2;
set(findall(gcf,'-property','FontSize'),'FontSize',fontsize)
set(findall(gcf,'-property','LineWidth'),'LineWidth',linewidth)


%%
path_string = '/Users/user/Google Drive/Amir Bashan/COC TEX/';
str = strcat(path_string,'compositional_data',num2str(round(1000*rand)),'_iterations',num2str(num_iterations));
saving_png_pdf(str)
% saveas(gcf,strcat(path_string,'compositional_data',num2str(round(1000*rand)),'_iterations',num2str(num_iterations),'.png'));

%%



%% (2) - data from corr matrix
tic
num_iterations = 2;
num_cells_vec = 10:10:100;
CorrMatrixStd = 20;
[ stat_mean_corr , stat_bcR, stat_COC ] = MeanCorr_bcR_COC_FromCorrmatrix(n, CorrMatrixStd, num_cells_vec, num_iterations);
toc

%%
figure
subplot(2,1,1);
plot_errorbar(num_cells_vec, stat_mean_corr, 'num cells', 'CM')
xlabel([]); xticks([]);
% subplot(2,1,2);
% plot_errorbar(num_cells_vec, stat_COC, 'num cells', 'COC')
subplot(2,1,2);
plot_errorbar(num_cells_vec, stat_bcR, 'num cells', 'bcD')
fontsize = 16; linewidth = 2;
set(findall(gcf,'-property','FontSize'),'FontSize',fontsize)
set(findall(gcf,'-property','LineWidth'),'LineWidth',linewidth)

%%
path_string = '/Users/user/Google Drive/Amir Bashan/COC TEX/';
str = strcat(path_string,'CMnumcells',num2str(round(1000*rand)),'_iterations',num2str(num_iterations));
saving_png_pdf(str)

%% save to pdf & png

saving_png_pdf(str)
