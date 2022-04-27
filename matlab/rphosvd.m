% Random Projection HOSVD
% Author: Wei-Chien Liao (liaowei@msu.edu)
% Tested by Shihab Shahriar Khan (khanmd@msu.edu)
% input:
%   X: a data tensor of size I_1 x I_2 x ... x I_N
%   R: a vector containg multilinear ranks of the core tensor, note that
%      R_n <= I_n for all n=1, 2, ... , N.
% ouput:
%   S: a core tensor of size R_1 x R_2 x ... x R_N
%   Q: a cell containing N factor matrices, factor Q{n} has size I_n x R_n
function [S,Q] = rphosvd(X,R)
Z = X;
N = ndims(X);
X_size = size(X);
Q = cell(N,1);
S = X;
for n = 1:N
    W = double(tenmat(Z , n))*randn(prod(X_size)/X_size(n) , R(n));
    [Q{n},~] = qr(W,0);
    S = ttm(S,Q{n}.',n);
end