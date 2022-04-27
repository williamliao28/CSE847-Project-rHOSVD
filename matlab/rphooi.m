% Random Projection HOOI
% Author: Wei-Chien Liao (liaowei@msu.edu)
% Tested by Shihab Shahriar Khan (khanmd@msu.edu)
% input:
%   X: a data tensor of size I_1 x I_2 x ... x I_N
%   R: a vector containg multilinear ranks of the core tensor, note that
%      R_n <= I_n for all n=1, 2, ... , N.
% ouput:
%   S: a core tensor of size R_1 x R_2 x ... x R_N
%   Q: a cell containing N factor matrices, factor Q{n} has size I_n x R_n
function [S,Q] = rphooi(X,R,maxiter)
Z = X;
N = ndims(X);
X_size = size(X);
Q = cell(N,1);
for n = 1:N
    Q{n} = randn(X_size(n) , R(n));
end

iter = 0;
while iter < maxiter
    for n = 1:N
        S = Z;
        for m = 1:N
            if m ~= n
                S = ttm(S,Q{m}.',m);
            end
        end
        W = double(tenmat(S , n))*randn(prod(X_size)/X_size(n) , R(n));
        [Q{n},~] = qr(W,0);
    end
    S = ttm(S,Q{N}.',N);
    iter = iter +1;
end