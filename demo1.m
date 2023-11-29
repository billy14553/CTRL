clear;
generateSimulatedData 
ms = mean(XSh);
mt = mean(XTh);
obj = uCTRL(ms,mt,15);
mts = filtering(mt,obj.h,obj.delabg);

%plot
f = figure('name', "uCTRL_demo",'Color','w');
plot(ms,'linewidth',1.5);
hold on;
plot(mt,'linewidth',1.5);
plot(mts,'--','linewidth',1.5);
xlabel("Index");
ylabel("Intensity (a.u.)");
SCIPlot;
legend(["$\bar{\mathbf{x}}_s$","$\bar{\mathbf{x}}_t$","$\hat{\bar{\mathbf{x}}}_s$"],'Interpreter','latex');
