% Basic Randomized SVD
% Author: Wei-Chien Liao (liaowei@msu.edu)
% input:
%   X: a data matrix of size I x J
%   R: target rank
%   p: oversampling parameter
%   q: power iteration parameter
% output:
%   SVD factor matrices: U (size I x R), S (size R x R), V (size R x J)
function [U, S, V] = rsvd(X,R,p,q)
Q = rsi(X, p+R, q);
B = Q.'*X;
[U1, S1, V1] = svds(B,R);
U1 = Q*U1;
U = U1(:,1:R);
S = S1(1:R,1:R);
V = V1(:,1:R);