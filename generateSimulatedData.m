


colortable = obtainCorlor();

rng(2023);
%gernerate Simulated data

m = 80;  %number of samples
x = -5:0.1:24.9;
p = length(x); %numer of wavelength
lw = 1.5;
x1 = 1:p;
s1 = normpdf(x, 0, 1);
s2 = normpdf(x, 5, 2);
s3 = normpdf(x, 15, 4);
wavelength= 1:1:p;
waveLabel = "WaveLength (a.u.)";
figure;
plot(wavelength,s1,'linewidth',lw);
hold on;
plot(wavelength,s2,'linewidth',lw);
plot(wavelength,s3,'linewidth',lw);
xlabel(waveLabel);
ylabel("Intensity (a.u.)");
SCIPlot
updateLabel(wavelength,0);
S = [s1;s2;s3];
Y = normrnd(5,1,m,3);
X = Y*S;
%add noise (background noise(1%) and white noise(0.1%))
abg1 = generateAddictivebackgroud(x);
abg2 = generateAddictivebackgroud(x);
wn1 =  normrnd(0,1,m,length(x))/1000;
wn2 =  normrnd(0,1,m,length(x))/1000;
%Source Domain and Target Domain
XSh = filtering(X,[0;ones(9,1)/9],0);
XTh = filtering(X,[0;rand(9,1)/9],0);
XSh = XSh+ repmat(abg1/100,80,1);
XTh = XTh+ repmat(abg2/100,80,1);
XSh = XSh + wn1;
XTh = XTh + wn2;
%add noise to response matrix (1%)
Y = Y + normrnd(0,1,m,3)/100;
figure;
plot(XSh','r');
hold on;
plot(XTh','g');
