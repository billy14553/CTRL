function [xs] = generateVirtualsample(Xs,Ys,yt,eps)
%generateVirtualsample 
%  Xs  source domain data         -- matrix
%  Ys  source domain response --matrix or vector
%  yt   target domain response --scalar
[m,n] = size(Xs);
cov = 1/(m-1)*Ys'*Ys;
incov = inv(cov);
%if size(Ys,2)>1
   % Ys = Ys(:,1);
%end
%y = Ys;
%stdy = std(y)*eps;
a = zeros(m,1);
for i = 1:1:m
% a(i) = exp(-(yt-y(i))^2/(2*stdy*stdy)); 
a(i) = exp(-(yt-Ys(i,:))*incov*(yt-Ys(i,:))'/(2*eps*eps)); 
end
a = a/sum(a);
xs = (Xs)'*a;
xs = xs';
end

