clear; close all;

fileID = fopen('matrix_vals', 'r');
A = fscanf(fileID, '%f');

n = sqrt(length(A));
A = reshape(A,n,n);
[~,p] = chol(-A)

eigs = eig(-A);