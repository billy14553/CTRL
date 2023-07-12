%default parameters

testratios = 0.3;
fold = 5;
LV = 3;
compArray = 1:1:19;
espArray = [10^3,10^2,10^1,10^0,10^-1,10^-2,10^-3,10^-4];
espArrayLabel = {"10^3","10^2","10^1","10^0","10^{-1}","10^{-2}","10^{-3}","10^{-4}"};
rhoArray = [-10^9,-10^8,-10^7,-10^6,-10^5,-10^4,-10^3,-10^2,-10^1,10^0,10^1,10^2,10^3,10^4,10^5,10^6,10^7,10^8,10^9];
eps = 0.01;
task = "simlated_data";
tasknum = abs(char(task));
visual = 'on';
ws = 5:2:65;
p = size(XSh,2);
wavelength= 1:1:p;
waveLabel = "WaveLength (a.u.)";
myesp = 0.01;
