function [Xh] = filtering(X,h,ab)
%filtering  x by filter h
%   X original signal, one row represents one sample
%   h filter, column vector
%  ab addictive backgroud noise
%  output
%  filtered signal Xh
m = size(X,1);
n =  size(X,2);
Xh = zeros(m,n);
N = length(h)-1;
C = cell(N,1);
D = eye(n);
 K = zeros(n,N-1);
c1 = X(:,1);
Xe = [repmat(c1,1,N-1) X];
for i = 1:1:N
    C{i} = zeros(n+N-1,n);
    C{i}(i:i+n-1,:) = D;
end

for j = 1:1:m
    x = Xe(j,:);   
    for i = 1:1:N
        K(:,i) = x*C{i};
    end
    Xh(j,:) = ([ones(n,1) K]*h)'+ab;
    %Xh(j,1:N) = mean(Xh(j,1:N));
end

end

