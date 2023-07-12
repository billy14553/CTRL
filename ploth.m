clear;close all;

% str = 'corn'
% load corn;
% X1 = m5spec.data;
% X2 = mp5spec.data;
% 
% y = propvals.data(:,1);

load nir_data;

str = 'gasline'
X1 = spec1(:,1:300);
X2 = spec2(:,1:300);
y = conc(:,3);

[m,n] = size(X1);
%X1 = imresize(X1,[m,floor(0.3*n)],'nearest');
%X2 = imresize(X2,[m,floor(0.3*n)],'nearest');


[XL,YL,XS,YS,BETA] = plsregress(X2,y,7);

%pid = 35;
%x1 = X1(pid,:);
%s1 = X2(pid,:);
%n = length(x1);
%N = 9;
% meanx = mean(X2);
% sspec = stdfir(X1,meanx,N,1);
% figure;
% subplot(1,2,1);
% plot(X2','g');
% hold on;
% 
% plot(X1','b')
% plot(sspec','r');
% subplot(1,2,2);
% y1 = [ones(80,1) X1]*BETA;
% y2 = [ones(80,1) sspec]*BETA;
% plot(y,'g');
% hold on;
% plot(y2,'r');
% plot(y1,'b');
% legend('真实y','FIR处理的m5样本预测','mp5样本预测');





%return;
%[stdmp,stdvp] = stdgen(s1,x1,1);
%disp(norm(stdmp));
%figure;
%plot(s1,'g')
%hold on;
%plot(stdvp,'r')
%return;
%figure;

%[B,FitInfo] = lasso(X2,y);
%BETA = [FitInfo.Intercept(59);B(:,59)];
%predict = [ones(80,1) X2]*BETA;
%MSE = mean((predict-y).^2);
%disp(MSE);
%plot(predict,'r'); hold on;
%plot(y,'g')
 
%figure;

pid = 6;
%x1 = mean(X1);
%s1 = mean(X2);
x1 =  X1(pid,:);
s1 = X2(pid,:);
n = length(x1);
N = 9;
C = [];
D = eye(n);
K = zeros(n,N);
xe = [zeros(1,N-1) x1];
size(xe);
for i = 1:1:N
    C{i} = zeros(n+N-1,n);
    C{i}(i:i+n-1,:) = D;
    K(:,i) = xe*C{i};
end
 
K = [ones(n,1) K];
h1  = pinv(K'*K)*K'*s1';
%plot(h);
figure;
subplot(1,2,1);
plot(K*h1,'r');
hold on;
plot(s1,'g');
plot(x1,'b');
%plot(s1'-K*h,'k');
%title('利用sample 1(m5->mp5)估计卷积核');
legend('mp5样本卷积后的光谱','m5样本','mp5样本');
%return;
diff1 = s1'-K*h1;
abg = dot([0;diff1],BETA);

%figure;
h2 = learningFilter(s1,x1,N);
xts1 = filtering(x1,h2,0);
%xts1(1:N) = mean(xts1(1:N)) ;
adb = s1-xts1;
subplot(1,2,2);
plot(xts1,'r');
hold on;
plot(s1,'g');
plot(x1,'b');
%plot(K*h,'k');

%return;
 %return;



figure;
for j = 1:1:6
 idx = (j+1)*3;   
x2 = X1(idx,:);
s2 = X2(idx,:);

xts2 = filtering(x2,h2,adb);
figure;
plot(xts2,'r');
hold on;
plot(s2,'g');
plot(x2,'b');
title(strcat('样本 ', num2str(idx)));
legend('mp5样本卷积后的光谱','m5样本','mp5样本','Location','southeast')
end
len = length(y);
for j = 1:1:len
 idx = j;   
x2 = X1(idx,:);
x2e = [zeros(1,N-1) x2];
K = [];
for i = 1:1:N
    C{i} = zeros(n+N-1,n);
    C{i}(i:i+n-1,:) = D;
    K(:,i) = (x2e*C{i})';
end
K = [ones(n,1) K]; 
y1(j) = dot([1 X1(j,:)],BETA);
y2(j) = dot([1;K*h1],BETA)+abg;

end
figure;
plot(y,'g');
hold on;
plot(y2,'r');
plot(y1,'b');
legend('真实y','转移后预测','转移前预测');

res = y' - y2;


rmse =sqrt(mean(res.^2));
RPD = std(y)/rmse;
std(y)
rmse
disp(RPD)

SST = sum((y-mean(y)).^2);
SSR = sum(res.^2)

R2 = 1-SSR/SST;

disp(sqrt(1/(1-R2)));

txt = sprintf('data:%s,RPD:%.2f',str,RPD);
title(txt);



len = length(y);
[XL,YL,XS,YS,BETA,PCTVAR,MSE] = plsregress(X2,conc,7);
y1 = [];
y2 = [];
XSh = X2;
Y = conc;
XTh = X1;
testratios = 0.3;
LV = 7;
[TrainIdx,TestIdx] = duplex(XSh,floor(m*testratios));
TrainXs = XSh(TrainIdx,:);
TrainYs = Y(TrainIdx,:);
TestXs = XSh(TestIdx,:);
TestYs = Y(TestIdx,:);
TrainXt = XTh(TrainIdx,:);
TrainYt = Y(TrainIdx,:);
TestXt = XTh(TestIdx,:);
TestYt = Y(TestIdx,:);
[XL,YL,XS,YS,BETA,PCTVAR,MSE] = plsregress(TrainXs,TrainYs,LV);
figure;
for j = 1:1:len
     idx = j;   
     x2 = X1(idx,:);
     xts2 = filtering(x2,h,adb);
     plot(x2,'b');
     hold on;
     plot(X2(idx,:),'g');
     plot(xts2,'r');
     y1(j,:) = [1 x2]*BETA;
     y2(j,:) = [1 xts2]*BETA ;
end
return;
for i = 1:1:size(conc,2)
    figure;
    plot(y1(:,i),'b');
    hold on;
    plot(y2(:,i),'r');
    plot(conc(:,i),'g');
end
 

