% Author: Wei-Chien Liao (liaowei@msu.edu)
% tests for various randomized HOSVD algorithms
clear;
clc;
%% Test 1: Random Low Multilinear Rank Tensors
% randomly generate Gaussian core tensor
S = tensor(randn(20,40,30));
Q = cell(3,1);
Q{1} = randn(1000,20);
Q{2} = randn(1000,40);
Q{3} = randn(1000,30);
X = ttm(S, Q, 1:3);
[S1,Q1,R] = sthosvd(X,0.1);