%% pdf7 - all noises 
format compact
clear
clc

multi_weight = 2; % weights from the interval (0, multi_weight)%%%%%%%%% 1?
multi_start = 1; % initial conditions x0 from the interval (0, multi_start)
T = 100; % relaxation time

%% creating A_template

% num_cells = 50; % number of cells = size(M,1) 200
n = 100; % number of genes = size(M,2) 400 

k_act = 2;
p_act = k_act/n;
A_template = multi_weight*sprand(n,n,p_act);% 0 < aij < 1
A_template(1:n+1:end) = 0;% No self activation

%% WeightNoise + MeasureNoise

%% CorrCells and bcd phase diagrams

Bvec_template = ones(1,n);
num_cells = 50;

num_itr = 10;
bcd_itr = 5;
sigma_weights_vec = 0:0.05:1; %0.05
sigma_noise_vec = 0:0.05:1; %0.05

tic

[bcdMatrix, VariabilityCellsMatrix ] = CorrCells_bcd_WeightNoise_MeasureNoise( num_cells, A_template, Bvec_template, multi_start, T, sigma_weights_vec, sigma_noise_vec, num_itr, bcd_itr);

toc

% %% saving data 
% 
% randnum = round(rand()*100);
% file_name_Cell = strcat('VariabilityCellsMatrix',num2str(randnum),'numRuns',num2str(num_itr),'bcdRuns',num2str(bcd_itr),'numCells',num2str(num_cells),'numGenes',num2str(n),'.csv');
% csvwrite(file_name_Cell,VariabilityCellsMatrix)
% file_name_bcd = strcat('bcdMatrix',num2str(randnum),'numRuns',num2str(num_itr),'bcdRuns',num2str(bcd_itr),'numCells',num2str(num_cells),'numGenes',num2str(n),'.csv');
% csvwrite(file_name_bcd,bcdMatrix)

%% import specific data #40min

sigma_weights_vec = 0:0.05:1; %0.05
sigma_noise_vec = 0:0.05:1; %0.05
VariabilityCellsMatrix = csvread("VariabilityCellsMatrix78numRuns10bcdRuns5numCells50numGenes100.csv");
bcdMatrix = csvread("bcdMatrix78numRuns10bcdRuns5numCells50numGenes100.csv");

Font_Size = 22;
contour_corr_cell = 0.1;
plot_VariabilityCells_bcdMatrix(sigma_weights_vec, sigma_noise_vec,'\sigma_{dynamics}', '\sigma_{measure}', VariabilityCellsMatrix, bcdMatrix, contour_corr_cell, Font_Size)



%% %%%%%%%%%%%%%%%%%%%%%% old!

%% import specific data % 3min?

% n = 100; Bvec_template = ones(1,n); num_cells = 50; num_itr = 2;
% sigma_weights_vec = 0:0.1:1; sigma_noise_vec = 0:0.1:1;
CorrCellsMatrix = csvread("CorrCellsMatrix_____.csv");
bcdMatrix = csvread("bcdMatrix______.csv");

%% ploting phase diagrams 

contour_corr_cell = 0.92;
plot_CorrCells_bcdMatrix(sigma_weights_vec, sigma_self_vec, '\sigma_{weights}', '\sigma_{measure}', CorrCellsMatrix, bcdMatrix, contour_corr_cell)

%% WeightNoise + SelfActivationNoise_Bmatrix

%% one M

Bvec_template = ones(1,n);%1
% Bvec_template = 0.5+rand(1,n);%2
% Bvec_template = 2*rand(1,n);%3

num_cells = 10;
sigma_weights = 0.5;
sigma_self = 0.5;

M_diff_noiseweight_B = M_NoiseWeightAndB_function( num_cells, A_template, Bvec_template, multi_start, T, sigma_weights, sigma_self);

% figure;
% plot(M_diff_noiseweight_B')
% title('M example')
% colorbar;

%% CorrCells and bcd phase diagrams

Bvec_template = ones(1,n);%1
% Bvec_template = 0.5+rand(1,n);%2
% Bvec_template = 2*rand(1,n);%3

num_cells = 50;
sigma_weights_vec = 0:0.1:1; 
sigma_self_vec = 0:0.1:1;
num_itr = 20;

tic
[CorrCellsMatrix, bcdMatrix ] = CorrCells_bcd_WeightNoise_BNoise( num_cells, A_template, Bvec_template, multi_start, T, sigma_weights_vec, sigma_self_vec, num_itr);
toc

%% different A different results? %50 min 

n = 50;
k_act = 2;
p_act = k_act/n;
    
num_cells = 30;

sigma_weights_vec = 0:0.1:1; 
sigma_self_vec = 0:0.1:1;
num_itr = 20;

num_runs = 1:3;

Bvec_template = ones(1,n);%1
% Bvec_template = 0.5+rand(1,n);%2
% Bvec_template = 2*rand(1,n);%3

tic

for i = num_runs
    
    i
    
    A_template = multi_weight*sprand(n,n,p_act);% 0 < aij < 1
    % A_template(1:n+1:end) = rand(n,1);% for all genes - old self activation
    A_template(1:n+1:end) = 0;% No self activation
    
    [CorrCellsMatrix, bcdMatrix ] = CorrCells_bcd_WeightNoise_BNoise( num_cells, A_template, Bvec_template, multi_start, T, sigma_weights_vec, sigma_self_vec, num_itr);

end

toc

