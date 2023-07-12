function [Xp] = prepocessing(X)
%prepocessing preprocessing spectrum
%Keep 99%percent variance of X
[m,n] = size(X);
mx    = mean(X);
mcx   = (X-mx(ones(m,1),:));
[U,S,V] = svd(mcx);
s = diag(S);
var = cumsum(s)/sum(s);
for i =1:1:length(s)
    if(var(i)>0.99)
       s(i+1:end) = 0;
       break;
    end
end
S(1:length(s),1:length(s)) = diag(s);
Xp = U*S*V'+mx(ones(m,1),:);
end

