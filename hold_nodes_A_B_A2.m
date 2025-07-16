function xdot = hold_nodes_A_B_A2(t, x, A, B, A2, A2_twovec, holdnodes, c_hill)%
% new = with C matrix

c_f = 1;
% c_hill = 1;

% self dynamics: 
f = @(x) - x.^c_f;
% activation:
g1 = @(x) (x.^c_hill)./(1+(x.^c_hill));
% inhibition:
g2 = @(x) 1./(1+(x.^c_hill));

% % activation + inhibition
% xdot = f(x)+A*g1(x)+B*g2(x);

% % activation + inhibition + A2
xdot = f(x)+A*g1(x)+B*g2(x)+ A2.*g1(x(A2_twovec(:,1))).*g1(x(A2_twovec(:,2)));

xdot(holdnodes) = 0;

end

