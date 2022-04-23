% Randomized Subspace Iteration
% input:
%   X: a matrix of size I x J
%   L: an integer indicating the number of vectors to sample
%   q: power iteration parameter
% output:
%   Q: I x L orthonormal matrix whose range approximates the range of X
function Q = rsi(X, L, q)
[~ , J] = size(X);
omega = randn(J, L);
Y = X*omega;
[Q , ~] = qr(Y,0);
for n = 1:q
    Y = X.'*Q;
    [Q , ~] = qr(Y,0);
    Y = X*Q;
    [Q , ~] = qr(Y,0);
end