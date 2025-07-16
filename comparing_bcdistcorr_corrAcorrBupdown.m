%% comparing bcdistcorr & corrAcorrB_updown

format compact
clear
clc

num_cells = 300; % number of cells = size(M,1)%200
n = 400; % number of genes = size(M,2)%400

multi_weight = 1; % weights from the interval (0, multi_weight)
multi_start = 1; % initial conditions x0 from the interval (0, multi_start)
T = 40; % relaxation time

% for DOC
string_type_corr = 'Spearman';
% string_type_corr = 'Pearson';
overlap_type = 1;% Overlap with weights
threshold = 1e-3;% threshold for gene expression
rlowess_span = 0.1;% smooth (not used for now)

%% THE NETWORK

k_act = 2;
p_act = k_act/n;
A_template = multi_weight*sprand(n,n,p_act);% 0 < aij < 1
A_template(1:n+1:end)=0;% no loops

%%
font_size = 14;% size font for figs

%%
tic

% percent_damage_vec = 0.01:0.01:0.5;
percent_damage_vec = [0.005,0.05,0.5];

COC = zeros(size(percent_damage_vec));
bcR = zeros(size(percent_damage_vec));

for i = 1:length(percent_damage_vec)
    
    p = percent_damage_vec(i);
    M_switching_links  = M_sameA_damagelinks_completelydiffweight_function( A_template, p, num_cells, multi_start, T);
    COC(i) = corrAcorrB_updown(M_switching_links, string_type_corr);
    [bcR(i),~,~,~] = bcdistcorr(M_switching_links(:,1:n/2),M_switching_links(:,(n/2 + 50):n));
    
end

figure;
plot(percent_damage_vec,COC);
hold on
plot(percent_damage_vec,bcR);
xlabel('damage \textit{p}','Interpreter', 'LaTeX')
legend('COC','bcR')
set(gca,'fontsize', font_size);

toc
