% Sequentially Truncated HOSVD
% Author: Wei-Chien Liao (liaowei@msu.edu)
% input:
%   X:       a data tensor of size I_1 x I_2 x ... x I_N
%   epsilon: tolerance
% ouput:
%   S: a core tensor of size R_1 x R_2 x ... x R_N
%   Q: a cell containing N factor matrices, factor Q{n} has size I_n x R_n
%   R: a vector containg multilinear ranks of the core tensor, note that
%      R_n <= I_n for all n=1, 2, ... , N.
function [S,Q,R] = sthosvd(X,epsilon)
S = X;
N = ndims(X);
X_size = size(X);
Q = cell(N,1);
for n = 1:N
    [Q{n},~,~] = svdsketch(double(tenmat(S , n)),epsilon/N);
    S = ttm(S,Q{n}.',n);
end
R = size(S);