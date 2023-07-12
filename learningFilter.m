function [h,delabg] = learningFilter(xs,xt,N)
%learningFilter  estimate filter h
%   xs sample from source domain
%   xt  sample from target domain
%   N window size of filter
%   output
%   filter h
%   residual additive background noise
n = length(xs);
C = cell(N,1);
D = eye(n);
K = zeros(n,N);
xe = [zeros(1,N-1) xt];
for i = 1:1:N
    C{i} = zeros(n+N-1,n);
    C{i}(i:i+n-1,:) = D;
    K(:,i) = xe*C{i};
end
K = [ones(n,1) K];
K = prepocessing(K);
h  = pinv(K'*K)*K'*xs';
xts = filtering(xt,h,0);
delabg = xs-xts;
end

