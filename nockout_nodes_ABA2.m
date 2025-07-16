function xdot = nockout_nodes_ABA2(t,x,A,B,A2,i)
% new = with C matrix

c_f = 1;
c_h = 1;

% self dynamics: 
f = @(x) - x.^c_f;
% activation:
g1 = @(x) (x.^c_h)./(1+(x.^c_h));
% inhibition:
g2 = @(x) 1./(1+(x.^c_h));


% ode
xdot = f(x)+A*g1(x)+B*g2(x)+(A2.genes_with_activationDim2).*(g1(x(A2.firstvec))).*(g1(x(A2.secondvec)));

xdot(i) = 0;

end