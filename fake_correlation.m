%% fake correlation (without gene-regulation) 
% (1) compositional data & (2) data from corr matrix
format compact
clear
clc

%%
num_cells = 100; % number of cells = size(M,1)
n = 200; % number of genes = size(M,2)

%%
%% (1) - compositional data
% here I use corrAcorrB_updown but independent blocks of data.
tic
num_iterations = 1;%20
D_vec = 0.5:1.5:3.5;%0.1:0.1:4
YNstring = 'N';
% YNstring = 'Y';%for ploting CorrMatrix
[ stat_mean_corr , stat_bcR, stat_COC ] = MeanCorr_bcR_COC_FromNormExp(num_cells, n, D_vec, num_iterations,YNstring);
toc
%%
fontsize = 14;
figure
subplot(2,1,1);
plot_errorbar(D_vec, stat_mean_corr, 'D', 'CM')
xlabel([]); xticks([]);
subplot(2,1,2);
plot_errorbar(D_vec, stat_COC, 'D', 'COC')
% subplot(2,1,2);
% plot_errorbar(D_vec, stat_bcR, 'D', 'bcR')
set(findall(gcf,'-property','FontSize'),'FontSize',fontsize)
%%
path_string = '/Users/user/Google Drive/Amir Bashan/COC TEX/';
str = strcat(path_string,'CMD',num2str(round(1000*rand)),'_iterations',num2str(num_iterations));
% saveas(gcf,strcat(path_string,'CMD',num2str(round(1000*rand)),'_iterations',num2str(num_iterations),'.png'));

%% all in one plot
tic
num_iterations = 2;%20
D_vec = 0.1:0.1:2;
sp_D_vec = [0.1,1,2.5];
plotall_compositionaldata(num_cells, n, D_vec, num_iterations, sp_D_vec)
toc

%%



%% (2) - data from corr matrix
tic
num_iterations = 20;
num_cells_vec = 10:10:100;
CorrMatrixStd = 20;
[ stat_mean_corr , stat_bcR, stat_COC ] = MeanCorr_bcR_COC_FromCorrmatrix(n, CorrMatrixStd, num_cells_vec, num_iterations);
toc

%%
fontsize = 14;
figure
subplot(2,1,1);
plot_errorbar(num_cells_vec, stat_mean_corr, 'num cells', 'CM')
xlabel([]); xticks([]);
subplot(2,1,2);
plot_errorbar(num_cells_vec, stat_COC, 'num cells', 'COC')
% subplot(2,1,2);
% plot_errorbar(num_cells_vec, stat_bcR, 'num cells', 'bcR')
set(findall(gcf,'-property','FontSize'),'FontSize',fontsize)

%%
path_string = '/Users/user/Google Drive/Amir Bashan/COC TEX/';
str = strcat(path_string,'CMnumcells',num2str(round(1000*rand)),'_iterations',num2str(num_iterations));

%% save to pdf & png

saving_png_pdf(str)
