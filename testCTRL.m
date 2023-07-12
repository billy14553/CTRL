clear;close all;
load nir_data.mat;    
XSh = spec2;
XTh = spec1;
XSh = savgol(XSh,9,2,2);
XTh =  savgol(XTh,9,2,2);
Y = conc;
initParameter;
wavelength= lamda;
waveLabel = "WaveLength (nm)";
task = "corn_spec2_to_spec1";
tasknum = abs(char(task));
disp(task);
LV = 7;

p = size(Y,2);
m = size(XSh,1);
%Splite the data
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

sel = 3;
xt = TrainXt(sel,:);
yt = TrainYt(sel,:);
xs = TrainXs(sel,:);
h = learningFilter(xs,xt,9);
xts = filtering(xt,h,0);
delabg = xs - xts;
y1 = [];
y2 = [];
len = size(TestXt,1);
for j = 1:1:len
     idx = j;   
     x2 = TestXt(idx,:);
     xts2 = filtering(x2,h,delabg);
     y1(j,:) = [1 x2]*BETA;
     y2(j,:) = [1 xts2]*BETA ;
end
Y = TestYt;
RMSEP1 = [];
for i = 1:1:size(Y,2)
    figure;
    plot(y1(:,i),'b');
    hold on;
    plot(y2(:,i),'r');
    plot(Y(:,i),'g');
    str = sprintf("NOCT RMSEP:%f,CT RMSEP:%f",RMSE(y1(:,i)-Y(:,i)),RMSE(y2(:,i)-Y(:,i)));
    title(str);
    RMSEP1(i) = RMSE(y2(:,i)-Y(:,i));
    
end
disp(RMSEP1);
return;
xt1 = TrainXt(1,:);
xs1 = TrainXs(1,:);
 xts1 = filtering(xt1,h,delabg);
 plot(xs1);
 hold on;
 plot(xt1);
 plot(xts1);
 yt
 [1 xs1]*BETA
 [1 xt1]*BETA
[1 xts1]*BETA
[1 xts1]*BETA1
 