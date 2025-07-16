%% pdf6 - WeightNoise + SelfActivationNoise_Bmatrix
format compact
clear
clc

% cd \Users\user\'Google Drive'\'Amir Bashan'\'Gene Matlab'\
% path_string = '/Users/user/Google Drive/Amir Bashan/COC TEX/';

% num_cells = 50; % number of cells = size(M,1)%200
% n = 100; % number of genes = size(M,2)%400

multi_weight = 1; % weights from the interval (0, multi_weight)
multi_start = 1; % initial conditions x0 from the interval (0, multi_start)
T = 40; % relaxation time

font_size = 14;% size font for figs

%%

%% Random network - activation

n = 100;

k_act = 2;
p_act = k_act/n;
A_template = multi_weight*sprand(n,n,p_act);% 0 < aij < 1
% A_template(1:n+1:end) = rand(n,1);% for all genes - old self activation
A_template(1:n+1:end) = 0;% No self activation

% Bvec_template = 1;%1
Bvec_template = 0.5+rand(1,n);%2
% Bvec_template = 2*rand(1,n);%3


%% one M

sigma_weights = 0.8;
sigma_self = 0;
num_cells = 15;%%%
M_diff_noiseweight_B = M_NoiseWeightAndB_function( num_cells, A_template, Bvec_template, multi_start, T, sigma_weights, sigma_self);

figure;
plot(M_diff_noiseweight_B')
title('M example')
colorbar;

%% one phase diagram
tic

% sigma_weights_vec = 0:0.05:1;
% sigma_self_vec = 0:0.05:1;
% num_itr = 10;

num_cells = 50;
sigma_weights_vec = 0:0.1:1; 
sigma_self_vec = 0:0.1:1;
num_itr = 2;

[CorrCellsMatrix, ~, bcdMatrix, ~ ] = CorrCells_CorrCOC_bcd_CM_WeightNoise_BNoise( num_cells, A_template, Bvec_template, multi_start, T, sigma_weights_vec, sigma_self_vec, num_itr);

toc

%% diff num_cells and n
%~0.5h
tic
 
sigma_weights_vec = 0:0.1:1; 
sigma_self_vec = 0:0.1:1;
num_itr = 2;
k_act = 2;

for n = 100:50:200   
    n
    
    p_act = k_act/n;
    A_template = multi_weight*sprand(n,n,p_act);% 0 < aij < 1
    A_template(1:n+1:end)=0;% No self activation

    Bvec_template = 0.5+rand(1,n);
    
    for num_cells = 50:25:100
        num_cells
        
        [CorrCellsMatrix,~, bcdMatrix, ~ ] = CorrCells_CorrCOC_bcd_CM_WeightNoise_BNoise( num_cells, A_template, Bvec_template, multi_start, T, sigma_weights_vec, sigma_self_vec, num_itr);        
        title(strcat('bcd for n = ',num2str(n),' and numcells = ',num2str(num_cells)))

    end
end

toc

%% runs with the same parameters
%~1.5h
tic
 
sigma_weights_vec = 0:0.5:1; 
sigma_self_vec = 0:0.5:1;
num_itr = 2;

n = 50;
num_cells = 6;%%%

k_act = 2;
p_act = k_act/n;
    
for loop_itr = 1:2
    
    loop_itr
    
    A_template = multi_weight*sprand(n,n,p_act);% 0 < aij < 1
    A_template(1:n+1:end)=0;% No self activation

    Bvec_template = 0.5+rand(1,n);
%     Bvec_template = ones(1,n);
    % Bvec_template = 2*rand(1,n);

    [CorrCellsMatrix, ~, bcdMatrix, ~ ] = CorrCells_CorrCOC_bcd_CM_WeightNoise_BNoise( num_cells, A_template, Bvec_template, multi_start, T, sigma_weights_vec, sigma_self_vec, num_itr);        
    title(strcat('itr = ',num2str(loop_itr),' bcd for n = ',num2str(n),' and numcells = ',num2str(num_cells)))

end

toc


%%
path_string = '/Users/user/Google Drive/Amir Bashan/COC TEX/diff runs n and numcells/';

figs = get(0,'Children');
str_rand_num = num2str(round(100*rand));
path_string = strcat(path_string,str_rand_num);

for i = 1:length(figs)
    newpath = strcat(path_string,num2str(i));
    saveas(figs(i),newpath,'png')
end

%%
contour_corr_cell = 0.92;
plot_CorrCells_bcdMatrix(sigma_weights_vec, sigma_self_vec, 'Sigma Weights', 'Sigma Self', CorrCellsMatrix, bcdMatrix, contour_corr_cell)
% plot_CorrCells_CorrCOC_bcdMatrix_CM(sigma_weights_vec, sigma_self_vec, 'Sigma Weights', 'Sigma Self', CorrCellsMatrix, CorrCCCMatrix, bcdMatrix, CM, contour_corr_cell)

%% saving
set(gca,'FontSize', 18)


