fold = 5;
indices = crossvalind('Kfold', m, fold); 
RMSECV_PDS = zeros(length(ws),3);
for i = 1:1:length(ws)
    RMSEV = [0,0,0];
    for j = 1:1:fold
        test = (indices==j);
        train = ~test;
        F_TrainXs =  TrainXs(train,:);
        F_TrainXt =  TrainXt(train,:);
        F_TrainYs =  TrainYs(train,:);
        F_TrainYt =  TrainYt(train,:);         
        F_TestXs  =  TrainXs(test,:);
        F_TestXt  =  TrainXt(test,:);
        F_TestYs =  TrainYs(test,:);
        F_TestYt =  TrainYt(test,:); 
         
        [stdmp,stdvp] = stdgen(F_TrainXs,F_TrainXt,ws(i));
        cspec2p = stdize(F_TestXt,stdmp,stdvp);
        m1 = size(cspec2p,1);
        Yp = [ones(m1,1) cspec2p]*BETA;
        RMSEV(1) = RMSEV(1) + RMSE(Yp(:,1)-F_TestYt(:,1));
        RMSEV(2) = RMSEV(2) + RMSE(Yp(:,2)-F_TestYt(:,2));
        RMSEV(3) = RMSEV(3)+  RMSE(Yp(:,3)-F_TestYt(:,3));
    end
    RMSECV_PDS(i,1) = RMSEV(1)/fold;
    RMSECV_PDS(i,2) = RMSEV(2)/fold;
    RMSECV_PDS(i,3) = RMSEV(3)/fold;
end


% clear;
% res = [];
% a = 1.0;
% b = -5.1;
% c = -4.0;
% d = 2.0;
% x = -101:1:100;
%  
% for i = 1:1:length(x)
%     y(i) = a*x(i)^3+b*x(i)^2+c*x(i)+d;
%     if(i<2) 
%         continue;
%     end
%     if(y(i)==0)
%         res = [res x(i)];
%         continue;
%     end
%     if(y(i)*y(i-1)<0)
%         
%             yp   = y(i);
%             s = x(i-1);
%             e = x(i);
%             j = x(i);
%             while abs(yp)>10^(-9)
%                  j = (s+e)/2;
%                  yp = a*j^3+b*j^2+c*j+d;
%                  if yp*y(i)>0
%                      e = j;
%                  else
%                      s = j;
%                  end
%             end
%             res = [res j];
%         end
% end
%  
%  disp(res);
%  plot(x,y)