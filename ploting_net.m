function [] = ploting_net(A_template, percent_damage)

n = size(A_template,1);
index_template = find(A_template>0);
num_damage = round(percent_damage*length(index_template));

dt = 2*pi/n;
t = dt:dt:2*pi;
x = cos(t); y = sin(t);

A = A_template;
% A(A>0) = rand([1,nnz(A)])';
% A(sign(A)<0) = 0;

% nodes_added:
whole_vec = 1:n^2;
whole_vec(index_template)=[];%unconnected links in the original network
links_added = randsample(whole_vec,num_damage);%adding links from the unconnected links
% nodes_removed:
links_removed = randsample(index_template,num_damage);

A(links_removed)= 0;

B = zeros(size(A));
B(links_added) = rand(size(links_added));
B(1:n+1:end)=0;% no loops

% figure;
G_A = digraph(A);
plot(G_A,'NodeColor','k','MarkerSize',0.5,'LineWidth',1.1,'ArrowSize',1,'XData',x,'YData',y,'NodeLabel',{})
hold on
G_B = digraph(B);
plot(G_B,'NodeColor','k','MarkerSize',0.5,'LineWidth',1.6,'ArrowSize',1,'XData',x,'YData',y,'NodeLabel',{})
set(gca,'visible','off')

end

