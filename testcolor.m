% x = 1:1:100;
% figure;
% for i = 1:1:7
%     plot(x+i*10,'Color',colortable(i,:),'LineWidth',1);
%     hold on;
% end
% legend on;

load corn;

X1 = mp5spec.data;
X2 = mp6spec.data;

for i = 5:30:255
 X3 = stdfir(X1,mean(X2),i);
 figure;
 plot(X1','b');
 hold on;
 plot(X2','g');
 plot(X3','r');
end
 
 