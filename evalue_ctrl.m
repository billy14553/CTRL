%draw data
%pause(600);
colortable = obtainCorlor();
[m,n] = size(XSh);
% if n>50
%     XSh = imresize(XSh,[m,floor(0.3*n)],'nearest');
%     XTh = imresize(XTh,[m,floor(0.3*n)],'nearest');
% end
[m,n] = size(XSh);
f = figure('name',strcat(task,"_source_domain_spectra"),'visible',visual);
plot(XSh(1:2:m,:)','linewidth',1.5);
xlabel(waveLabel);
updateLabel(wavelength,0);
ylabel("Intensity (a.u.)");
SCIPlot;
MySaveFig(f,f.Name);
f = figure('name',strcat(task,"_target_domain_spectra"),'visible',visual);
plot(XTh(1:2:m,:)','linewidth',1.5);
xlabel(waveLabel);
updateLabel(wavelength,0);
ylabel("Intensity (a.u.)");
SCIPlot;
MySaveFig(f,f.Name);

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
[m,n] = size(TrainXs); %update m
indices = crossvalind('Kfold', m, fold);
save('dataset');
ms = mean(TrainXs);
mt = mean(TrainXt);
h = learningFilter(ms,mt,9);
mts = filtering(mt,h,0);
delabg = ms - mts;
f = figure('name',strcat(task,"_mean_spectra_transfer"),'visible',visual);
plot(ms,'linewidth',1.5);
hold on;
plot(mt,'linewidth',1.5);
plot(mts,'--','linewidth',1.5);
xlabel(waveLabel);
updateLabel(wavelength,0);
ylabel("Intensity (a.u.)");
SCIPlot;
legend(["$\bar{\mathbf{x}}_s$","$\bar{\mathbf{x}}_t$","$\hat{\bar{\mathbf{x}}}_s$"],'Interpreter','latex');
MySaveFig(f,f.Name);
%draw tranfered spectra
XTs = filtering(TrainXt,h,delabg);
f = figure('name',strcat(task,"_transfered_spectra_onesample"),'visible',visual);
disp(f.Name);
plot(TrainXt(floor(m/2),:),'linewidth',1.5);
hold on;
plot(TrainXs(floor(m/2),:),'linewidth',1.5,'Color','k');
plot(XTs(floor(m/2),:),'--','linewidth',1.5);
legend(["${\mathbf{x}}_t$","${\mathbf{x}}_s$","${\hat{\mathbf{x}}}_s$"],'Interpreter','latex');

xlabel(waveLabel);
updateLabel(wavelength,0);
ylabel("Intensity (a.u.)");
SCIPlot;
MySaveFig(f,f.Name);

f = figure('name',strcat(task,"_transfered_spectra"),'visible',visual);
plot(XTs(1:2:m,:)','linewidth',1.5);
xlabel(waveLabel);
updateLabel(wavelength,0);
ylabel("Intensity (a.u.)");
SCIPlot;
MySaveFig(f,f.Name);
%test by the filter obtain by mean spectrum
m1 = size(TestXt,1);
XTs = filtering(TestXt,h,delabg);
Yp = [ones(m1,1) XTs]*BETA;
mfRMSEP = zeros(1,p);
mfRDP = zeros(1,p);
for i = 1:1:p
    mfRMSEP(i) = RMSE(Yp(:,i)-TestYt(:,i));
    mfRDP(i) = std(TestYt(:,i))/mfRMSEP(i);
end

f = figure('name',strcat(task,"_mean_spectra_rmsep"),'visible',visual);
disp("mfRMSEP");
disp(mfRMSEP);
scatter(TestYt(:,1),Yp(:,1),'o','filled');
xlabel("y");
ylabel("yhat");
SCIPlot;
MySaveFig(f,f.Name);
disp(mean(mfRMSEP));
xt = TrainXt(6,:);
yt = TrainYt(6,:);
norm_dbg = zeros(length(ws),length(espArray));
espRMSEP = zeros(length(ws),length(espArray));
for k = 1:1:length(ws)
    for i = 1:1:length(espArray)
        eps = espArray(i);
        xs = generateVirtualsample(TrainXs,TrainYs,yt,eps);
        h = learningFilter(xs,xt,ws(k));
        xts = filtering(xt,h,0);
        delabg = xs - xts;
        norm_dbg(k,i) = norm(delabg);
        XTs = filtering(TestXt,h,delabg);
        Yp = [ones(m1,1) XTs]*BETA;
        for j = 1:1:p
            espRMSEP(k,i) = espRMSEP(k,i) + RMSE(Yp(:,j)-TestYt(:,j));
        end
    end
end

%ticklabel
yla = cell(1,length(ws));
for i = 1:1:length(ws)
    yla{i} = num2str(ws(i));
end


f = figure('name',strcat(task,"_esp_select"),'visible',visual);
if length(ws)==1
    plot(norm_dbg,'linewidth',1.5);
    ylabel("$||\Delta \mathbf{b}_g||$",'Interpreter','latex');
else
    imagesc(norm_dbg);
    ytk = yticks();
    ytk1 = cell(1,length(ytk));
    for i = 1:1:length(ytk)
       ytk1{i} = yla{ytk(i)};
    end
    ylabel("Windows Size");
    yticklabels(ytk1);
    colorbar;
end
xlabel("\epsilon");
xticklabels({"10^3","10^2","10^1","10^0","10^{-1}","10^{-2}","10^{-3}","10^{-4}"});
SCIPlot;
MySaveFig(f,f.Name);
if length(ws)>1
    [Xm,Ym] =  meshgrid(1:1:length(espArray),ws);
    f = figure('name',strcat(task,"_esp_select_mesh"),'visible',visual);
    mesh(Xm,Ym,norm_dbg,'linewidth',1);
    xlabel("\epsilon");
    xticklabels({"","10^2","10^0","10^{-2}","10^{-4}"});
    ylabel("Window Size");
    zlabel("$||\Delta \mathbf{b}_g||$",'Interpreter','latex');
    SCIPlot;
    grid on;
    view(54,32);
    MySaveFig(f,f.Name);
end

f = figure('name',strcat(task,"_esp_rmsep_fir"),'visible',visual);
if length(ws)==1
    plot(espRMSEP,'-o','linewidth',1.5,'MarkerFaceColor','auto','MarkerSize',2);
    ylabel("RMSEP");
else
    imagesc(espRMSEP);
    ylabel("Window Size");
    yticklabels(ytk1);
    colorbar;
end
xlabel("\epsilon");
xticklabels({"10^3","10^2","10^1","10^0","10^{-1}","10^{-2}","10^{-3}","10^{-4}"});
SCIPlot;
MySaveFig(f,f.Name);

tmp = min(norm_dbg(:));
[a,b] = find(norm_dbg==tmp);
RMSEP_CTLF = zeros(m,p);
RMSEP_NOCT = zeros(1,p);
RMSEP_SELF = zeros(1,p);
RDP_CTLF = zeros(m,p);
RDP_NOCT = zeros(1,p);
RDP_SELF = zeros(1,p);
for i = 1:1:m
    xt = TrainXt(i,:);
    yt = TrainYt(i,:);
    xs = generateVirtualsample(TrainXs,TrainYs,yt,myesp);
    h = learningFilter(xs,xt,floor(mean(a)));
    xts = filtering(xt,h,0);
    delabg = xs - xts;
    XTs = filtering(TestXt,h,delabg);
    Yp = [ones(m1,1) XTs]*BETA;
    for j = 1:1:p
        RMSEP_CTLF(i,j) = RMSE(Yp(:,j)-TestYt(:,j));
        RDP_CTLF(i,j) = std(TestYt(:,j))/RMSE(Yp(:,j)-TestYt(:,j));
    end
end
 

%figure;
TrainXTs = filtering(TrainXt,h,delabg);
f = figure('name',strcat(task,"_transfer_ctrl"),'visible',visual);
disp(f.Name);
plot(nan(),'-','Color',colortable(1,:),'LineWidth',1);
hold on
plot(nan(),'-','Color',colortable(5,:),'LineWidth',1);
plot(nan(),'--','Color',colortable(2,:),'LineWidth',1);
legend(["spectra from target domain","spectra from source domain","tranferred spectra"],'AutoUpdate','off');
 
plot(TrainXt(2:2:end,:)','Color',colortable(1,:),'linewidth',1);
 hold on;
 plot(TrainXs(2:2:end,:)','Color',colortable(5,:),'linewidth',1);
 plot(TrainXTs(2:2:end,:)','--','Color',colortable(2,:),'linewidth',1);
ylabel("Intensity (a.u.)");
 SCIPlot;
xlabel(waveLabel);
updateLabel(wavelength,0);
MySaveFig(f,f.Name);

Yp = [ones(m1,1) TestXt]*BETA;
for j = 1:1:p
    RMSEP_NOCT(1,j) = RMSE(Yp(:,j)-TestYt(:,j));
    RDP_NOCT(1,j) = std(TestYt(:,j))/RMSE(Yp(:,j)-TestYt(:,j));
end
Yp = [ones(m1,1) TestXs]*BETA;
for j = 1:1:p
    RMSEP_SELF(1,j) = RMSE(Yp(:,j)-TestYs(:,j));
    RDP_SELF(1,j) = std(TestYt(:,j))/RMSE(Yp(:,j)-TestYs(:,j));
end
disp("---------------RMSEP_CTLF_BEGIN-------------------");
disp(median(RMSEP_CTLF));
disp("---------------RMSEP_CTLF_END-------------------");
%return;
% RMSEPResult ="";
% RMSEPResult = RMSEPResult + sprintf("RMSEP_SELF:%.4f,%.4f,%.4f \n",RMSEP_SELF(1),RMSEP_SELF(2),RMSEP_SELF(3));
% RMSEPResult = RMSEPResult + sprintf("RMSEP_NOCT:%.4f,%.4f,%.4f \n",RMSEP_NOCT(1),RMSEP_NOCT(2),RMSEP_NOCT(3));
% mrmsep = mean(RMSEP_CTLF);
% RMSEPResult = RMSEPResult + sprintf("RMSEP_CTLF:%.4f,%.4f,%.4f \n",mrmsep(1),mrmsep(2),mrmsep(3));
% disp(RMSEPResult);

if exist('stdfir')==1
    f = errordlg('Install PLS TOOLBOX first','File Error');
    return;
end
%evaluate FIR on

RMSEP_FIR = zeros(length(ws),3);
RDP_FIR = zeros(length(ws),3);
for i = 1:1:length(ws)
    sspec = stdfir(TestXt,mean(TrainXs),ws(i));
    Yp = [ones(m1,1) sspec]*BETA;
    for j = 1:1:p
        RMSEP_FIR(i,j) = RMSE(Yp(:,j)-TestYt(:,j));
        RDP_FIR(i,j) = std(TestYt(:,j))/RMSE(Yp(:,j)-TestYt(:,j));
    end
end
f = figure('name',strcat(task,"_ws_rmsep_fir"),'visible',visual);
disp(f.Name);
plot(ws,RMSEP_FIR,'-o','linewidth',1.5,'MarkerFaceColor','auto','MarkerSize',2);
ylabel("RMSEP");
xlabel("Window Size");
SCIPlot;
MySaveFig(f,f.Name);

[a,b] = min(mean(RMSEP_FIR,2))
disp(sprintf("RMSEP_FIR: %f, best window:%d \n",a,ws(b)));

f = figure('name',strcat(task,"_transfered_spectra_fir"),'visible',visual);
disp(f.Name);
plot(nan(),'-','Color',colortable(1,:),'LineWidth',1);
hold on
plot(nan(),'-g','LineWidth',2);
plot(nan(),'--','Color',colortable(2,:),'LineWidth',1.0);

legend(["target spectra","reference spectrum","tranferred target spectra"],'AutoUpdate','off');
 
sspec = stdfir(TrainXt,mean(TrainXs),ws(b));
plot(TrainXt(2:2:m,:)','-','Color',colortable(1,:),'LineWidth',1);
hold on;
plot(mean(TrainXs),'-g','LineWidth',2);
plot(sspec(2:2:m,:)','--','Color',colortable(2,:),'LineWidth',1);

xlabel(waveLabel);
updateLabel(wavelength,0);
ylabel("Intensity (a.u.)");
SCIPlot;
MySaveFig(f,f.Name);

RMSECV_PDS = zeros(length(ws),p);
for i = 1:1:length(ws)
    RMSEV = zeros(1,p);
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
        for k = 1:1:p
            RMSEV(k) = RMSEV(k) + RMSE(Yp(:,k)-F_TestYt(:,k));
        end
    end
    for k = 1:1:p
        RMSECV_PDS(i,k) = RMSEV(k)/fold;
    end
end
f = figure('name',strcat(task,"_RMSECV_PDS"),'visible',visual);
disp(f.Name);
mrmsecv = mean(RMSECV_PDS,2);
h = plot(ws,mrmsecv,'-o','linewidth',1.5,'MarkerFaceColor','auto','MarkerSize',2);
ylabel("RMSECV");
xlabel("Window Size");
SCIPlot;
MySaveFig(f,f.Name);

%PDS can't  transfer by one reference
% for i = 2:5:size(TrainXs,1)
[~,i] = min(mrmsecv);
RMSEP_PDS = zeros(1,p);
RDP_PDS = zeros(1,p);
[stdmp,stdvp] = stdgen(TrainXs,TrainXt,ws(i));
if size(stdvp,1)>1
    stdvp = mean(stdvp);
end
m1 = size(TestXt,1);
cspec2p = stdize(TestXt,stdmp,stdvp);
Yp = [ones(m1,1) cspec2p]*BETA;
for i = 1:1:p
    RMSEP_PDS(i) =   RMSE(Yp(:,i)-TestYt(:,i));
    RDP_PDS(i) = std(TestYt(:,i))/RMSEP_PDS(i) ;
end
%disp("RMSEP_PDS");
%disp(RMSEP_PDS);
disp(task);
 %[status, cmdout] = system('E:\Users\xzh\anaconda3\python.exe pycalt.py')
% 
 load(strcat(task,"_result.mat"));
 disp("RMSEP_FIR & RDP");
 disp(RMSEP_FIR(LV,:));
 disp(RDP_FIR(LV,:));

 disp("RMSEP_DOP & RDP");
 disp(RMSEP_DOP);
 disp(RDP_DOP);
 disp("RMSEP_uDOP & RDP");
 disp(RMSEP_uDOP(LV,:));
 disp(RDP_uDOP(LV,:));
 disp("RMSEP_SST & RDP");
 disp(RMSEP_SST);
 disp(RDP_SST);
 disp("RMSEP_PDS & RDP");
 disp(RMSEP_PDS);
 disp(RDP_PDS);
 disp("RMSEP_NOCT & RDP");
 disp(RMSEP_NOCT);
  disp(RDP_NOCT);
 disp("RMSEP_SELF & RDP");
 disp(RMSEP_SELF);
 disp(RDP_SELF);
  disp("RMSEP_DIPLS & RDP");
 disp(RMSEP_DIPLS);
 disp(RDP_DIPLS);
 disp("RMSEP_CTLF & RDP");
 disp(median(RMSEP_CTLF));
 disp(median(RDP_CTLF));
 disp("mfRMSEP");
disp(mfRMSEP);
disp(mfRDP);
save(strcat(task,"_all_result"));

