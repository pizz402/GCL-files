function xdot = hold_nodes_A(t,x,A,holdnodes)
% new = with C matrix

c_f = 1;
c_h = 1;

% self dynamics: 
f = @(x) - x.^c_f;
% activation:
g1 = @(x) (x.^c_h)./(1+(x.^c_h));
% % inhibition:
% g2 = @(x) 1./(1+(x.^c_h));

xdot = f(x)+A*g1(x);

switch nargin
    
    case 3 % activation

    case 4 % activation + n_hold
        xdot(holdnodes) = 0;

end