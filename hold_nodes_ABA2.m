function xdot = hold_nodes_ABA2(t,x,A,Bvec,n_hold_vec)%,A2,n_hold)
% new = with C matrix

c_f = 1;
c_h = 1;

% self dynamics: 
f = @(x) - x.^c_f;
% activation:
g1 = @(x) (x.^c_h)./(1+(x.^c_h));
% % inhibition:
% g2 = @(x) 1./(1+(x.^c_h));

% ode
switch nargin
    
    case 3 % activation with noise on weights
        xdot = f(x)+A*g1(x);
    
    case 4 % activation with noise on weights + noise on degradation
        xdot = Bvec'.*f(x)+A*g1(x);
    
    case 5 % activation with noise on weights + noise on degradation + n_hold
        xdot = Bvec'.*f(x)+A*g1(x);
        xdot(n_hold_vec) = 0;
        
%     case 6 % activation + inhibition with noise on weights % B here is a matrix!!!!!!!!
%         xdot = f(x)+A*g1(x)+B*g2(x)+(A2.genes_with_activationDim2).*(g1(x(A2.firstvec))).*(g1(x(A2.secondvec)));
%         xdot(1:n_hold) = 0;

end

end