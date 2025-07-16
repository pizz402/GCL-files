% n = size(X,2);%genes
% m = size(X,1);%cells

%% mod() for Z2

x1 = [0 0 0];
x2 = [0 1 1];
x3 = [1 0 1];
x4 = [1 1 0];

X = [x1;x2;x3;x4];
X = repmat(X,[4,1]);
X = 0.01 + X.*(1+0.01*rand(size(X)));

corr(X)
new_bcdistcorr(X)


%% mod() for  Z3

r1 = [zeros(9,1);ones(9,1);2*ones(9,1)];

r2tag = [zeros(3,1);ones(3,1);2*ones(3,1)];
r2 = repmat(r2tag,[3,1]);

r3tag = [0;1;2];
r3 = repmat(r3tag,[9,1]);

r4 = mod(r1+r2+r3,3);

X = [r1,r2,r3,r4];

corr(X)
new_bcdistcorr(X)

%% sum

n = 5;%genes
m = 100;%cells

X = rand(m,n-1)-0.5;
xn = sum(X.^2,2);

X = [X,xn];

% X = X./sum(X,2);

new_bcdistcorr_itr(X, 20)

%% multiple

n = 6;%genes
m = 10;%cells

X = rand(m,n-1);
xn = prod(X,2)./mean(X,2);

X = [X,xn];

X = X./sum(X,2);

imagesc(corr(X))
new_bcdistcorr_itr(X, 20)

