function [C] = generateCmatrix(n,N)
%generateCmatrix Generate C Matrix 
%  n is the length of signal
%  N is the windows of filter
%  output
%  C is a cell, C{i} denotes the ith C Matrix 
C = cell(N,1);
D = eye(n);
K = zeros(n,N);
for i = 1:1:N
    C{i} = zeros(n+N-1,n);
    C{i}(i:i+n-1,:) = D;
end
end

