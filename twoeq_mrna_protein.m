function [zdot] = twoeq_mrna_protein(t,z,A,holdnodes)%n_hold

% dx/dt = mi*fi(y)-lambda_rna*xi
% dy/dt = ri*xi-lambda_prot*yi

num_genes = length(z)/2;
x = z(1:num_genes);
y = z(num_genes+1:2*num_genes);

lambda_rna = 1;
lambda_prot = 1;
r = 1;

% activation:
f = @(y) y./(1+y);

% self degradation: 
degx = @(x) - lambda_rna.*x;
degy = @(y) - lambda_prot.*y;

xdot = A*f(y)+degx(x);

switch nargin
    
    case 3 
       
    case 4 
        xdot(holdnodes) = 0;
        
end

ydot = r.*x+degy(y);

zdot = [xdot;ydot];

end