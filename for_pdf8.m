%% pdf8 - comparing to no interactions 
% + bcdVec and CMVec for different num_cells
% + bcdVec and CMVec for normalization with D

format compact
clear
clc

multi_weight = 2; % weights from the interval (0, multi_weight)
multi_start = 1; % initial conditions x0 from the interval (0, multi_start)
T = 100; % relaxation time

%% creating A_template and Bvec_template

% num_cells = 50; % number of cells = size(M,1) 200
n = 200; % number of genes = size(M,2) 400 

k_act = 2;%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
p_act = k_act/n;

self_inter_cases = 4;%%%%%%%%%%%%%%%%%%%%%%%

switch self_inter_cases
      
    case 2 % only interactions 0 < aij < multi_weight
        A_template = multi_weight*sprand(n,n,p_act);
        A_template(1:n+1:end) = 0;
    
    case 4 % loops 0 < aii < multi_weight & interactions 0 < aij < multi_weight
        A_template = multi_weight*sprand(n,n,p_act);
        A_template(1:n+1:end) = multi_weight*rand(n,1);
    
end

Bvec_template = ones(1,n);%1

%% CorrCellsVec and bcdVec for different p noise
% for here self_inter_cases = 4 and k_act = 0/2

tic

num_cells = 100;
num_holds = 10; 
num_itr = 20;%20
bcd_itr = 5;%5

p_inter_vec = linspace(0,1,10);

Font_Size = 30;
Line_Width = 2;

[VariabilityCellsVec, bcdVec ] = CorrCells_bcd_InterSelfNoise( num_cells, A_template, Bvec_template, T, p_inter_vec, num_itr, multi_weight, num_holds, bcd_itr, Font_Size, Line_Width);

toc

% %% saving
% 
% randnum = round(rand()*100);
% 
% file_name_Cell = strcat('VariabilityCellsVec',num2str(randnum),'numRuns',num2str(num_itr),'numitbcd',num2str(bcd_itr),'numCells',num2str(num_cells),'numGenes',num2str(n),...
%     'kact',num2str(k_act),'multiweight',num2str(multi_weight),'numholds',num2str(num_holds),'.csv');
% csvwrite(file_name_Cell,VariabilityCellsVec)
% 
% file_name_Cell = strcat('bcdVec',num2str(randnum),'numRuns',num2str(num_itr),'numitbcd',num2str(bcd_itr),'numCells',num2str(num_cells),'numGenes',num2str(n),...
%     'kact',num2str(k_act),'multiweight',num2str(multi_weight),'numholds',num2str(num_holds),'.csv');
% csvwrite(file_name_Cell,bcdVec)

%% reading from saved data

VariabilityCellsVec66 = csvread("VariabilityCellsVec66numRuns20numitbcd5numCells100numGenes200kact2multiweight2numholds10.csv");%%%run here
bcdVec66 = csvread("bcdVec66numRuns20numitbcd5numCells100numGenes200kact2multiweight2numholds10.csv");%%%run here
VariabilityCellsVec43 = csvread("VariabilityCellsVec43numRuns20numitbcd5numCells100numGenes200kact0multiweight2numholds10.csv");
bcdVec43 = csvread("bcdVec43numRuns20numitbcd5numCells100numGenes200kact0multiweight2numholds10.csv");

%% specific run 

p_vec = linspace(0,1,10);

VariabilityCells_interself = VariabilityCellsVec66;
bcdVec_interself = bcdVec66;

VariabilityCells_self = VariabilityCellsVec43;
bcdVec_self = bcdVec43;

mean_VariabilityCells_interself = mean(VariabilityCells_interself,2);
std_VariabilityCells_interself = std(VariabilityCells_interself,0,2);
mean_bcd_interself = mean(bcdVec_interself,2);
std_bcd_interself = std(bcdVec_interself,0,2);

mean_VariabilityCells_self = mean(VariabilityCells_self,2);
std_VariabilityCells_self = std(VariabilityCells_self,0,2);
mean_bcd_self = mean(bcdVec_self,2);
std_bcd_self = std(bcdVec_self,0,2);

%% ploting specific run

Font_Size = 30;
Line_Width = 2;

figure;
errorbar(p_vec,mean_VariabilityCells_interself,std_VariabilityCells_interself,'LineWidth',Line_Width);
xlabel('p')
ylabel('Cell-to-Cell Variability')
set(gca,'FontSize', Font_Size)

figure;
errorbar(p_vec,mean_VariabilityCells_self,std_VariabilityCells_self,'LineWidth',Line_Width);
xlabel('p')
ylabel('Cell-to-Cell Variability')
set(gca,'FontSize', Font_Size)

figure;
errorbar(p_vec,mean_bcd_interself,std_bcd_interself,'LineWidth',Line_Width);
xlabel('p')
ylabel('GCL')
set(gca,'FontSize', Font_Size)

figure;
errorbar(p_vec,mean_bcd_self,std_bcd_self,'LineWidth',Line_Width);
xlabel('p')
ylabel('GCL')
%ylim([-0.05,0.05])
set(gca,'FontSize', Font_Size)


%% bcdVec and CMVec for different num_cells
% for here self_inter_cases = 4 and k_act = 2

tic

num_holds = 20; 
num_itr = 20;
bcd_itr = 5;

% num_cells_vec = 10;
num_cells_vec = linspace(20,110,10);
p = 0;

[CorrCellsVec, bcdVec, CMVec ] = CorrCells_bcd_CM_InterSelfNoise( num_cells_vec, A_template, Bvec_template, T, p, num_itr, multi_weight, num_holds, bcd_itr);

toc

% % saving 
% 
% randnum = round(rand()*100);
% 
% file_name_Cell = strcat('CMVec',num2str(randnum),'numRuns',num2str(num_itr),'numitbcd',num2str(bcd_itr),'numGenes',num2str(n),...
%     'kact',num2str(k_act),'multiweight',num2str(multi_weight),'numholds',num2str(num_holds),'.csv');
% csvwrite(file_name_Cell,CMVec)
% 
% file_name_Cell = strcat('bcdVec',num2str(randnum),'numRuns',num2str(num_itr),'numitbcd',num2str(bcd_itr),'numGenes',num2str(n),...
%     'kact',num2str(k_act),'multiweight',num2str(multi_weight),'numholds',num2str(num_holds),'.csv');
% csvwrite(file_name_Cell,bcdVec)

%% reading from saved data - specific run

CMVec = csvread("VariabilityCellsVec66numRuns20numitbcd5numCells100numGenes200kact2multiweight2numholds10.csv");%%%run here
bcdVec = csvread("bcdVec66numRuns20numitbcd5numCells100numGenes200kact2multiweight2numholds10.csv");%%%run here

%%

Font_Size = 30;
Line_Width = 2;

num_cells_vec = linspace(20,110,10);

mean_bcd = mean(bcdVec');
std_bcd = std(bcdVec');

mean_CM = mean(CMVec');
std_CM = std(CMVec');

figure;
errorbar(num_cells_vec,mean_bcd,std_bcd,'LineWidth',Line_Width);
xlabel('num cells')
ylabel('GCL')
set(gca,'FontSize', Font_Size)

figure;
errorbar(num_cells_vec,mean_bcd,std_bcd,'LineWidth',Line_Width);
xlabel('num cells')
ylabel('GCL')
set(gca,'FontSize', Font_Size)

figure;
errorbar(num_cells_vec,mean_CM,std_CM,'LineWidth',Line_Width);
xlabel('num cells')
ylabel('\langle  Corr(i,j) \rangle')
set(gca,'FontSize', Font_Size)


%% bcdVec and CMVec for normalization with D

num_cells = 100; % number of cells = size(M,1)
n = 200; % number of genes = size(M,2)

num_iterations = 20;%20
D_vec = 0.1:0.05:1;
bcd_itr = 5;
std_measurement = 0.2;

% [ stat_mean_corr , stat_bcR, stat_COC ] = MeanCorr_bcR_COC_FromNormExp(num_cells, n, D_vec, num_iterations,'N');
[ stat_mean_corr , stat_bcR ] = MeanCorr_bcR_FromNormExp(num_cells, n, D_vec, num_iterations ,bcd_itr, std_measurement);

figure
plot_errorbar(D_vec, stat_mean_corr, 'heterogeneity, D', '\langle  Corr(i,j) \rangle')
ylim([0,1])
figure
plot_errorbar(D_vec, stat_bcR, 'heterogeneity, D', 'bcD')
