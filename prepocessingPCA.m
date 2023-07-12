function [Xp] = prepocessingPCA(X,maxpc)
%prepocessing preprocessing spectrum
%Keep 99%percent variance of X
[m,n] = size(X);
mx    = mean(X);
mcx   = (X-mx(ones(m,1),:));
[U,S,V] = svd(mcx);
s = diag(S);
s(maxpc+1:end) = 0;
S(1:length(s),1:length(s)) = diag(s);
Xp = U*S*V'+mx(ones(m,1),:);
end

