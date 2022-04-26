% Randomized Pass-Efficient HOSVD
% input:
%   X: a data tensor of size I_1 x I_2 x ... x I_N
%   R: a vector containg multilinear ranks of the core tensor, note that
%      R_n <= I_n for all n=1, 2, ... , N.
%   k, s: sketching parameters, s > k. (both of size N x 1)
% ouput:
%   S: a core tensor of size R_1 x R_2 x ... x R_N
%   Q: a cell containing N factor matrices, factor Q{n} has size I_n x R_n
function [S,Q] = rpet(X, R, k, s)
N = ndims(X);
X_size = size(X);
Q = cell(N,1);
omega = cell(N,1);
H = X;
for n = 1:N
    [Q{n},~] = qr(double(tenmat(X, n))*randn(prod(X_size)/X_size(n) , k(n)),0);
    omega{n} = randn(s(n),X_size(n));
    H = ttm(H,omega{n},n);
end
% truncate thefactor matrices
for n = 1:N
    Q{n} = Q{n}(:,1:R(n));
end
% recover core tensor
S = H;
for n = 1:N
    S_size = size(S);
    S = lsqminnorm(omega{n}*Q{n},double(tenmat(S, n)));
    S = reshape(tensor(S),[S_size(1:n-1), R(n), S_size(n+1:end)]);
end