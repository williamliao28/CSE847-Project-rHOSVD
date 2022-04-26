% Randomized Sequentially Truncated HOSVD
% input:
%   X: a data tensor of size I_1 x I_2 x ... x I_N
%   R: a vector containg multilinear ranks of the core tensor, note that
%      R_n <= I_n for all n=1, 2, ... , N.
% ouput:
%   S: a core tensor of size R_1 x R_2 x ... x R_N
%   Q: a cell containing N factor matrices, factor Q{n} has size I_n x R_n
function [S,Q] = rsthosvd(X,R)
S = X;
N = ndims(X);
Q = cell(N,1);
for n = 1:N
    [Q{n}, ~, ~] = rsvd(double(tenmat(S, n)), R(n), 10, 2);
    S = ttm(S,Q{n}.',n);
end