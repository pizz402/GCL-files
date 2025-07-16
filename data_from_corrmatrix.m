function M = data_from_corrmatrix(n,std,num_cells)
% std - only natural numbers
W = randn(n,std);
S = W*W' + diag(rand(1,n));
S = diag(1./sqrt(diag(S))) * S * diag(1./sqrt(diag(S)));

S(1:n+1:end) = 1;

M = copularnd('gaussian',S,num_cells);

end