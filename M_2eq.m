%% create cells from 2-eq

% dx/dt = mi*fi(y)-lambda_rna*xi
% dy/dt = ri*xi-lambda_prot*yi

% x - mRNA concentration
% y - protein concentration

% M = x matrix
clear
clc

%%

num_cells = 100; % number of cells = size(M,1)%100
n = 200; % number of genes = size(M,2)%200

multi_weight = 1; % weights from the interval (0, multi_weight)
multi_start = 1; % initial conditions x0 from the interval (0, multi_start)
T = 40; % relaxation time

%%

k_act = 2;
p_act = k_act/n;
A_template = multi_weight*sprand(n,n,p_act);% 0 < aij < 1
% A_template(1:n+1:end)=rand(n,1);% for all genes - self activation
A_template(1:n+1:end)=0;% No self activation


%% bcd(sigma_weights,n_hold) 

tic

num_runs = 5;%5

sigma_weights_vec = 0:0.1:1;
num_weights = length(sigma_weights_vec);

n_hold_vec = 0:2:6;%0:1:20
num_nhold = length(n_hold_vec);

CorrCellsMatrix = zeros(num_nhold,num_weights);
CM = zeros(num_nhold,num_weights);
bcd_matrix = zeros(num_nhold,num_weights);

for nhold_index = 1:num_nhold
    
    nhold_index
    
    n_hold = n_hold_vec(nhold_index);

    for sigma_weights_index = 1:num_weights

        sigma_weights_index
        
        sigma_weights = sigma_weights_vec(sigma_weights_index);

        A_new = A_template;
        A_vec_noise_weights = 1+ sigma_weights*(2*rand(nnz(A_new),1)-1);
        A_new(A_new>0) = A_new(A_new>0).*A_vec_noise_weights;

        [ M_diff_starts ] = M_2eq_function( num_cells, n_hold, A_new, multi_start, T );
        
        [rho, ~] = corr(M_diff_starts','Type','Spearman');
        CorrCells = sum(sum(abs(triu(rho,1))))/(num_cells*(num_cells-1)/2);
        CorrCellsMatrix(nhold_index,sigma_weights_index) = CorrCells;

        [rho, ~] = corr(M_diff_starts,'Type','Spearman');
        mean_corr = sum(sum(abs(triu(rho,1))))/(n*(n-1)/2);
        CM(nhold_index,sigma_weights_index) = mean_corr;

        for run = 1:num_runs
        
            bcd = new_bcdistcorr(M_diff_starts);
            bcd_matrix(nhold_index,sigma_weights_index) = bcd_matrix(nhold_index,sigma_weights_index)+bcd./num_runs;
            
        end

    end
end

toc


%% bcd(p_damage,n_hold) %322 sec

tic

num_runs = 5;%5

p_vec = 0:0.1:1;
num_p = length(p_vec);

n_hold_vec = 0:1:20;%0:1:20
num_nhold = length(n_hold_vec);

CorrCellsMatrix = zeros(num_nhold,num_p);
CM = zeros(num_nhold,num_p);
bcd_matrix = zeros(num_nhold,num_p);

for nhold_index = 1:num_nhold
    
    n_hold = n_hold_vec(nhold_index);

    for p_index = 1:num_p
        
        p_index
        p = p_vec(p_index);

        [ M_diff_starts ] = M_2eq_damagelinks_function( num_cells, p, n_hold, A_template, multi_start, T);
        
        [rho, ~] = corr(M_diff_starts','Type','Spearman');
        CorrCells = sum(sum(abs(triu(rho,1))))/(num_cells*(num_cells-1)/2);
        CorrCellsMatrix(nhold_index,p_index) = CorrCells;

        [rho, ~] = corr(M_diff_starts,'Type','Spearman');
        mean_corr = sum(sum(abs(triu(rho,1))))/(n*(n-1)/2);
        CM(nhold_index,p_index) = mean_corr;
            
        for run = 1:num_runs
            
            bcd = new_bcdistcorr(M_diff_starts);
            bcd_matrix(nhold_index,p_index) = bcd_matrix(nhold_index,p_index)+bcd./num_runs;
            
        end

    end
end

toc

%% saving data

% randnum = round(rand()*100);
% file_name_Cell = strcat('CorrCellsMatrix',num2str(randnum),'numRuns',num2str(num_runs),'numCells',num2str(num_cells),'numGenes',num2str(n),'.csv');
% csvwrite(file_name_Cell,CorrCellsMatrix)
% file_name_bcd = strcat('bcdMatrix',num2str(randnum),'numRuns',num2str(num_runs),'numCells',num2str(num_cells),'numGenes',num2str(n),'.csv');
% csvwrite(file_name_bcd,bcd_matrix)
% file_name_CM = strcat('CM',num2str(randnum),'numRuns',num2str(num_runs),'numCells',num2str(num_cells),'numGenes',num2str(n),'.csv');
% csvwrite(file_name_CM,CM)

%% import specific data % 4850 sec

k_act = 3; num_itr = 5; sigma_weights_vec = 0:0.05:1; n_hold_vec = 0:1:20;

CorrCellsMatrix = csvread("CorrCellsMatrix78numRuns5numCells100numGenes200.csv");
bcd_matrix = csvread("bcdMatrix78numRuns5numCells100numGenes200.csv");
CM = csvread("CM78numRuns5numCells100numGenes200.csv");

%%

plot_CorrCells_bcdMatrix_CM(sigma_weights_vec, n_hold_vec, 'Sigma Weights', 'n hold', CorrCellsMatrix, bcd_matrix, CM )

%%

plot_CorrCells_bcdMatrix_CM(p_vec, n_hold_vec, 'Precent Damage', 'n hold', CorrCellsMatrix, bcd_matrix, CM )

%%
ylimdown = 0;

figure;
plot(p_vec,CM(10,:))
xlabel('p damage');ylabel('CM');
limsy=get(gca,'YLim');
set(gca,'Ylim',[ylimdown limsy(2)]);

figure;
plot(p_vec,bcd_matrix(10,:))
xlabel('p damage');ylabel('bcd');
limsy=get(gca,'YLim');
set(gca,'Ylim',[ylimdown limsy(2)]);

figure;
plot(p_vec,CorrCellsMatrix(10,:))
xlabel('p damage');ylabel('Corr Cells');
limsy=get(gca,'YLim');
set(gca,'Ylim',[ylimdown limsy(2)]);


%% 
%%%%%%%%%%%%%%%%% 28.05.19

%% first noise -  holding genes on the same value
tic
n_hold_vec = 0:2:100;% values for n_hold
num_iterations = 5; % creating new cells every run
bcR = zeros(length(n_hold_vec),num_iterations);
CM = zeros(length(n_hold_vec),num_iterations);
for j = 1 : num_iterations
    j
    for i = 1:length(n_hold_vec)
        i
        n_hold = n_hold_vec(i);
        M_2eq_compleltelydiffweight_nhold = M_2eq_compleltelydiffweight_nhold_function( num_cells, n_hold, A_template, multi_start, T);
        [bcR(i,j), CM(i,j)] = bcR_CM_FromM(M_2eq_compleltelydiffweight_nhold);
    end
end
toc

%% saving
randnum = round(rand()*100);
file_name_bcd = strcat('bcdMatrix',num2str(randnum),'numRuns',num2str(num_iterations),'numCells',num2str(num_cells),'numGenes',num2str(n),'.csv');
csvwrite(file_name_bcd,bcR)
file_name_CM = strcat('CM',num2str(randnum),'numRuns',num2str(num_iterations),'numCells',num2str(num_cells),'numGenes',num2str(n),'.csv');
csvwrite(file_name_CM,CM)

%% import specific data % 820 sec num 31
addpath('\Users\user\Google Drive\Amir Bashan\results genes 280519\hold nodes\')
k_act = 2;
n_hold_vec = 0:2:100;% values for n_hold
num_iterations = 5; % creating new cells every run
bcR = csvread("bcdMatrix31numRuns5numCells100numGenes200.csv");
CM = csvread("CM31numRuns5numCells100numGenes200.csv");

%%
ploting_GCL_CM(n_hold_vec,'n hold',bcR,CM)

%% second noise -  link rewiring
tic
percent_damage_vec = 0:0.02:1; % values for damages
num_iterations = 5; % creating new cells every run
bcR = zeros(length(percent_damage_vec),num_iterations);
CM = zeros(length(percent_damage_vec),num_iterations);
for j = 1 : num_iterations
    j
    for i = 1:length(percent_damage_vec)
        i
        p = percent_damage_vec(i);
        M_2eq_compleltelydiffweight_damagelinks = M_2eq_compleltelydiffweight_damagelinks_function( num_cells, p, A_template, multi_start, T);
        [bcR(i,j), CM(i,j)] = bcR_CM_FromM(M_2eq_compleltelydiffweight_damagelinks);
    end
end
toc

%% saving (again)
randnum = round(rand()*100);
file_name_bcd = strcat('bcdMatrix',num2str(randnum),'numRuns',num2str(num_iterations),'numCells',num2str(num_cells),'numGenes',num2str(n),'.csv');
csvwrite(file_name_bcd,bcR)
file_name_CM = strcat('CM',num2str(randnum),'numRuns',num2str(num_iterations),'numCells',num2str(num_cells),'numGenes',num2str(n),'.csv');
csvwrite(file_name_CM,CM)

%% import specific data % 250 sec num 96
addpath('\Users\user\Google Drive\Amir Bashan\results genes 280519\rewiring\')
k_act = 2;
percent_damage_vec = 0:0.02:1; % values for damages
num_iterations = 5; % creating new cells every run
bcR = csvread("bcdMatrix96numRuns5numCells100numGenes200.csv");
CM = csvread("CM96numRuns5numCells100numGenes200.csv");

%%
ploting_GCL_CM(percent_damage_vec,'p',bcR,CM)
