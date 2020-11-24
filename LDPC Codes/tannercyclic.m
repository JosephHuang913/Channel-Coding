% Name: tannercyclic.m
% Description:
% Display the Tanner graph of a binary cyclic code

% Copyright (c) 2006. Robert Morelos-Zaragoza. All rights reserved.
clear all

% Length, dimension and generator polynomial as a vector (least first)
n=15; k=5; g = [1  0  1  0  0  1  1  0  1  1  1];

[parmat,genmat,k] = cyclgen(n,g,'nonsys');
h = parmat(1,1:k+1);        % Coefficients of PC polunomial (least first)

% Construct the extended parity-check matrix
H = [];
row = [ h zeros(1,n-k-1) ]; H = [ H; row ];
for i=1:n-1
    row_new(2:n) = row(1:n-1); row_new(1) = row(n);
    H = [H; row_new];
    row = row_new;
end

% Adjacency matrix as defined in Matlab
A = [ zeros(n,n) H'; H zeros(n,n) ];

xy = [];
for i=1:n
    xy = [ xy; i 2 ];                 % Locations of code nodes
end
for i=1:n
    xy = [ xy; i 1 ];                 % Locations of check nodes
end

gplot(A,xy)
axis([0 n+1 0.9 2.1])
title('Tanner graph of a binary cyclic code')