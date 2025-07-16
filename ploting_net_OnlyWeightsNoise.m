function [] = ploting_net_OnlyWeightsNoise(A_template, percent_damage)

n = size(A_template,1);
index_template = find(A_template>0);
num_damage = round(percent_damage*length(index_template));

dt = 2*pi/n;
t = dt:dt:2*pi;
x = cos(t); y = sin(t);

A = A_template;
% A(A>0) = rand([1,nnz(A)])';
% A(sign(A)<0) = 0;

links_withnoise = randsample(index_template,num_damage);
A(links_withnoise)= 0;

B = zeros(size(A));
B(links_withnoise) = rand(size(links_withnoise));

% % nodes_added:
% whole_vec = 1:n^2;
% whole_vec(index_template)=[];%unconnected links in the original network
% links_added = randsample(whole_vec,num_damage);%adding links from the unconnected links
% % nodes_removed:
% links_removed = randsample(index_template,num_damage);
% 
% A(links_removed)= 0;
% 
% B = zeros(size(A));
% B(links_added) = rand(size(links_added));
% B(1:n+1:end)=0;% no loops

% figure;
G_A = digraph(A);
plot(G_A,'NodeColor','k','MarkerSize',0.5,'LineWidth',2,'ArrowSize',1,'XData',x,'YData',y,'NodeLabel',{})%1.1
hold on
G_B = digraph(B);
plot(G_B,'NodeColor','k','MarkerSize',0.5,'LineWidth',2.5,'ArrowSize',1,'XData',x,'YData',y,'NodeLabel',{})%1.6
set(gca,'visible','off')

end

