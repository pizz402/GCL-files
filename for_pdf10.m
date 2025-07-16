%% pdf10 - improving pdf9 by increasing m
% figs 1 2 3 
% fig 2 - k = 2, sigma noise, sigma weigths=>p,
% => add hold nodes, add loops

format compact
clear
clc

multi_weight = 2; % weights from the interval (0, multi_weight)
multi_start = 1; % initial conditions x0 from the interval (0, multi_start)
T = 100; % relaxation time

%% creating A_template and Bvec_template

% num_cells = 50; % number of cells = size(M,1)
n = 200; % number of genes = size(M,2)

k_act = 2;%%%%%%%%%%%%%%%%%%%%%%%%%%%%% for figure 1 also k_act = 0;
p_act = k_act/n;

A_template = multi_weight*sprand(n,n,p_act);
A_template(1:n+1:end) = multi_weight*rand(n,1);

Bvec_template = ones(1,n);%degregation

%% figure 1 - VariabilityCellsVec and bcdVec for different p noise - for k=0,2

tic

num_cells = 100;
num_holds = 10; 
num_itr = 20;
bcd_itr = 50;

p_inter_vec = linspace(0,1,10);

Font_Size = 30;
Line_Width = 2;

[VariabilityCellsVec, bcdVec ] = CorrCells_bcd_InterSelfNoise( num_cells, A_template, Bvec_template, T, p_inter_vec, num_itr, multi_weight, num_holds, bcd_itr, Font_Size, Line_Width);

toc

% % saving
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

VariabilityCellsVec70 = csvread("VariabilityCellsVec70numRuns20numitbcd50numCells100numGenes200kact2multiweight2numholds10.csv");
bcdVec70 = csvread("bcdVec70numRuns20numitbcd50numCells100numGenes200kact2multiweight2numholds10.csv");
VariabilityCellsVec25 = csvread("VariabilityCellsVec25numRuns20numitbcd50numCells100numGenes200kact0multiweight2numholds10.csv");%run
bcdVec25 = csvread("bcdVec25numRuns20numitbcd50numCells100numGenes200kact0multiweight2numholds10.csv");%run

%% import specific run 

p_vec = linspace(0,1,10);

VariabilityCells_interself = VariabilityCellsVec70;
bcdVec_interself = bcdVec70;

VariabilityCells_self = VariabilityCellsVec25;
bcdVec_self = bcdVec25;

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

% figure;
% errorbar(p_vec,mean_VariabilityCells_interself,std_VariabilityCells_interself,'LineWidth',Line_Width);
% xlabel('p')
% ylabel('Cell-to-Cell Variability')
% set(gca,'FontSize', Font_Size)

figure;
errorbar(p_vec,mean_VariabilityCells_self,std_VariabilityCells_self,'LineWidth',Line_Width);
xlabel('p')
ylabel('Cell-to-Cell Variability')
set(gca,'FontSize', Font_Size)

% figure;
% errorbar(p_vec,mean_bcd_interself,std_bcd_interself,'LineWidth',Line_Width);
% xlabel('p')
% ylabel('GCL')
% set(gca,'FontSize', Font_Size)

figure;
errorbar(p_vec,mean_bcd_self,std_bcd_self,'LineWidth',Line_Width);
xlabel('p')
ylabel('GCL')
ylim([-0.05,0.05])
set(gca,'FontSize', Font_Size)


%% figure 2 - comparing sigma noise and p noise

tic 

num_cells = 100;
num_holds = 10;
num_itr = 20;
bcd_itr = 50;

p_inter_vec = linspace(0.25,1,20);%20
sigma_noise_vec = linspace(0,0.5,21);%21

% sigma_noise_vec_1 = sigma_noise_vec(1:10);
% sigma_noise_vec_2 = sigma_noise_vec(11:end);
% sigma_noise_vec = sigma_noise_vec_1;%%% RUN- OPTION 1
% sigma_noise_vec = sigma_noise_vec_2;%%% RUN - OPTION 2

[ VariabilityCellsMatrix, bcdMatrix  ] = VariabilityCells_bcd_pNoise_MeasureNoise( num_cells, A_template, Bvec_template, T, p_inter_vec, sigma_noise_vec, num_itr, multi_weight, num_holds, bcd_itr);

toc

% % saving
% 
% randnum = round(rand()*100);
% 
% file_name_Cell = strcat('VariabilityCellsMatrix',num2str(randnum),'numRuns',num2str(num_itr),'numitbcd',num2str(bcd_itr),'numCells',num2str(num_cells),'numGenes',num2str(n),...
%     'kact',num2str(k_act),'multiweight',num2str(multi_weight),'numholds',num2str(num_holds),'.csv');
% csvwrite(file_name_Cell,VariabilityCellsMatrix)
% 
% file_name_Cell = strcat('bcdMatrix',num2str(randnum),'numRuns',num2str(num_itr),'numitbcd',num2str(bcd_itr),'numCells',num2str(num_cells),'numGenes',num2str(n),...
%     'kact',num2str(k_act),'multiweight',num2str(multi_weight),'numholds',num2str(num_holds),'.csv');
% csvwrite(file_name_Cell,bcdMatrix)

%% import specific run numholds=10 % old pdf

p_inter_vec = linspace(0,1,20);
sigma_noise_vec = linspace(0,0.5,20);
% contour_corr_cell = 0.5;

VariabilityCellsMatrix = csvread("VariabilityCellsMatrix91numRuns15numitbcd5numCells100numGenes200kact2multiweight2numholds10.csv");
bcdMatrix = csvread("bcdMatrix91numRuns15numitbcd5numCells100numGenes200kact2multiweight2numholds10.csv");

VariabilityCellsMatrix = VariabilityCellsMatrix(:,6:end);
bcdMatrix = bcdMatrix(:,6:end);
p_inter_vec = p_inter_vec(6:end);

%% import specific run numholds=10 % new pdf 

p_inter_vec = linspace(0.25,1,20);%20
sigma_noise_vec = linspace(0,0.5,21);%21
sigma_noise_vec_1 = sigma_noise_vec(1:10);
sigma_noise_vec_2 = sigma_noise_vec(11:end);
sigma_noise_vec = [sigma_noise_vec_1,sigma_noise_vec_2];

VariabilityCellsMatrix_1 = csvread("VariabilityCellsMatrix77numRuns20numitbcd50numCells100numGenes200kact2multiweight2numholds10.csv");
bcdMatrix_1 = csvread("bcdMatrix77numRuns20numitbcd50numCells100numGenes200kact2multiweight2numholds10.csv");

VariabilityCellsMatrix_2 = csvread("VariabilityCellsMatrix47numRuns20numitbcd50numCells100numGenes200kact2multiweight2numholds10.csv");
bcdMatrix_2 = csvread("bcdMatrix47numRuns20numitbcd50numCells100numGenes200kact2multiweight2numholds10.csv");

VariabilityCellsMatrix = [VariabilityCellsMatrix_1;VariabilityCellsMatrix_2];
bcdMatrix = [bcdMatrix_1;bcdMatrix_2];

%% ploting specific run 

Font_Size = 22;

contour_corr_cell = 0.5;
plot_VariabilityCells_bcdMatrix(p_inter_vec, sigma_noise_vec,'p', '\sigma', VariabilityCellsMatrix, bcdMatrix, contour_corr_cell, Font_Size)

contour_corr_cell_vec = [0.45,0.5,0.55];
plot_bcd_logpdivsigma(p_inter_vec, sigma_noise_vec,'p', '\sigma', VariabilityCellsMatrix, bcdMatrix, contour_corr_cell_vec, Font_Size)


%% figure 3 - bcdVec and CMVec for different num_cells

tic

num_holds = 10; 
num_itr = 20;
bcd_itr = 50;

num_cells_vec = linspace(20,150,14);
p = 0.5;

[bcdVec, CMVec ] = bcd_CM_InterSelfNoise( num_cells_vec, A_template, Bvec_template, T, p, num_itr, multi_weight, num_holds, bcd_itr);

toc

% % saving 
% 
% randnum = round(rand()*100);
% 
% file_name_Cell = strcat('CMVec',num2str(randnum),'numRuns',num2str(num_itr),'numitbcd',num2str(bcd_itr),'numGenes',num2str(n),...
%     'tenp',num2str(10*p),'kact',num2str(k_act),'multiweight',num2str(multi_weight),'numholds',num2str(num_holds),'.csv');
% csvwrite(file_name_Cell,CMVec)
% 
% file_name_Cell = strcat('bcdVec',num2str(randnum),'numRuns',num2str(num_itr),'numitbcd',num2str(bcd_itr),'numGenes',num2str(n),...
%     'tenp',num2str(10*p),'kact',num2str(k_act),'multiweight',num2str(multi_weight),'numholds',num2str(num_holds),'.csv');
% csvwrite(file_name_Cell,bcdVec)

%% import specific run - old pdf

num_cells_vec = linspace(20,150,14);

CMVec = csvread("CMVec81numRuns20numitbcd5numGenes200tenp5kact2multiweight2numholds10.csv");
bcdVec = csvread("bcdVec81numRuns20numitbcd5numGenes200tenp5kact2multiweight2numholds10.csv");

%% import specific run - new pdf

num_cells_vec = linspace(20,150,14);

CMVec = csvread("CMVec15numRuns20numitbcd50numGenes200tenp5kact2multiweight2numholds10.csv");
bcdVec = csvread("bcdVec15numRuns20numitbcd50numGenes200tenp5kact2multiweight2numholds10.csv");

%% ploting import specific run 

Font_Size = 36;
Line_Width = 2;

mean_bcd = mean(bcdVec');
std_bcd = std(bcdVec');
mean_CM = mean(CMVec');
std_CM = std(CMVec');

figure;
errorbar(num_cells_vec,mean_bcd,std_bcd,'LineWidth',Line_Width);
xlabel('number of cells')
ylabel('GCL')
set(gca,'FontSize', Font_Size)

figure;
errorbar(num_cells_vec,mean_CM,std_CM,'LineWidth',Line_Width);
xlabel('number of cells')
ylabel('\langle C \rangle')
set(gca,'FontSize', Font_Size)


%% figure 3 - bcdVec and CMVec for normalization with GAMMA

tic 

num_cells = 100; % number of cells = size(M,1)
n = 200; % number of genes = size(M,2)
k0 = 1;% arbitrary minimum 

num_itr = 20;%%%%%%%20
bcd_itr = 50;%%%%%%50

% D_vec = linspace(0.1,1,14);%%%%%%% I WAS HERE!!!
gamma_vec = linspace(1.1,2,20);
std_measurement = 0.2;

norm_yesno = 1;%%%%% 0/1

[ bcdVec, CMVec ] = MeanCorr_bcR_SF( num_cells, n, gamma_vec, bcd_itr, std_measurement, num_itr, norm_yesno, k0);
% [ bcdVec, CMVec ] = MeanCorr_bcR_Exp( num_cells, n, D_vec, bcd_itr,
% std_measurement, num_itr, norm_yesno);%%%%%%%%%%%%%%%%% I WAS HERE!!!

toc

mean_bcd = mean(bcdVec');
std_bcd = std(bcdVec');
mean_CM = mean(CMVec');
std_CM = std(CMVec');

Font_Size = 30;
Line_Width = 2;

figure;
errorbar(gamma_vec,mean_bcd,std_bcd,'LineWidth',Line_Width);
% ylim([-0.05,0.05]);
xlabel( '\gamma')
ylabel('GCL')
set(gca,'FontSize', Font_Size)

figure;
errorbar(gamma_vec,mean_CM,std_CM,'LineWidth',Line_Width);
% ylim([-0.05,0.05]);
xlabel( ' \gamma')
ylabel('\langle C \rangle')
set(gca,'FontSize', Font_Size)

% 
% % saving 
% 
% randnum = round(rand()*100);
% 
% file_name_Cell = strcat('CMVec',num2str(randnum),'numRuns',num2str(num_itr),'numitbcd',num2str(bcd_itr),'numCells',num2str(num_cells),'numGenes',num2str(n),...
%     'tenstdmeasurement',num2str(10*std_measurement),'k0',num2str(k0),'normyesno',num2str(norm_yesno),'.csv');
% csvwrite(file_name_Cell,CMVec)
% 
% file_name_Cell = strcat('bcdVec',num2str(randnum),'numRuns',num2str(num_itr),'numitbcd',num2str(bcd_itr),'numCells',num2str(num_cells),'numGenes',num2str(n),...
%     'tenstdmeasurement',num2str(10*std_measurement),'k0',num2str(k0),'normyesno',num2str(norm_yesno),'.csv');
% csvwrite(file_name_Cell,bcdVec)


%% import specific run - old pdf

D_vec = linspace(0.1,1,14);

% with normalization
CMVec_withnorm = csvread("CMVec51numRuns20numitbcd50numCells100numGenes200tenstdmeasurement2normyesno1.csv");
bcdVec_withnorm = csvread("bcdVec51numRuns20numitbcd50numCells100numGenes200tenstdmeasurement2normyesno1.csv");

% no normalization
CMVec_nonorm = csvread("CMVec97numRuns20numitbcd50numCells100numGenes200tenstdmeasurement2normyesno0.csv");
bcdVec_nonorm = csvread("bcdVec97numRuns20numitbcd50numCells100numGenes200tenstdmeasurement2normyesno0.csv");

%% import specific run - new pdf

gamma_vec = linspace(1.1,2,20);

% with normalization
CMVec_withnorm = csvread("CMVec85numRuns20numitbcd50numCells100numGenes200tenstdmeasurement5k01normyesno1.csv");
bcdVec_withnorm = csvread("bcdVec85numRuns20numitbcd50numCells100numGenes200tenstdmeasurement5k01normyesno1.csv");

% no normalization
CMVec_nonorm = csvread("CMVec53numRuns20numitbcd50numCells100numGenes200tenstdmeasurement5k01normyesno0.csv");
bcdVec_nonorm = csvread("bcdVec53numRuns20numitbcd50numCells100numGenes200tenstdmeasurement5k01normyesno0.csv");

%%

Font_Size = 30;
Line_Width = 2;

mean_bcd_withnorm = mean(bcdVec_withnorm');
std_bcd_withnorm = std(bcdVec_withnorm');
mean_CM_withnorm = mean(CMVec_withnorm');
std_CM_withnorm = std(CMVec_withnorm');

mean_bcd_nonorm = mean(bcdVec_nonorm');
std_bcd_nonorm = std(bcdVec_nonorm');
mean_CM_nonorm = mean(CMVec_nonorm');
std_CM_nonorm = std(CMVec_nonorm');

figure;
errorbar(gamma_vec,mean_bcd_withnorm,std_bcd_withnorm,'LineWidth',Line_Width);
hold on
errorbar(gamma_vec,mean_bcd_nonorm,std_bcd_nonorm,'LineWidth',Line_Width);
ylim([-0.05,0.05]);
xlabel( '\gamma')
ylabel('GCL')
legend('with normalization','no normalization','Location','best')
set(gca,'FontSize', Font_Size)

figure;
errorbar(gamma_vec,mean_CM_withnorm,std_CM_withnorm,'LineWidth',Line_Width);
hold on
errorbar(gamma_vec,mean_CM_nonorm,std_CM_nonorm,'LineWidth',Line_Width);
xlabel( '\gamma')
ylabel('\langle C \rangle')
legend('with normalization','no normalization','Location','best')
set(gca,'FontSize', Font_Size)

