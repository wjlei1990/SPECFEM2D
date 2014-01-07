clear all;close all;

for l=200:400:6000

figure;
temp=load(['./OUTPUT_FILES/snapshot_' num2str(l,'%4.4i')]);
x=temp(:,1);
y=temp(:,2);
k=temp(:,5);
if l==200; 
max_k=max(abs(k)); 
end;

NEX=500;
NEY=fix(NEX*(max(y)-min(y))/(max(x)-min(x)));
NEZ=30;
NGLLX=5;
NGLLY=5;

[X,Y]=meshgrid([min(x):(max(x)-min(x))/NEX:max(x)],[min(y):(max(y)-min(y))/NEY:max(y)]);

K=griddata(x,y,k,X,Y,'linear',{'QJ'});

pcolor(X,Y,K);
shading interp;
axis equal;
colorbar off;
caxis([-1 1]*max_k*0.3);
axis off;

hold on;
plot([0 1],[0.2 0.2],'k');
saveas(gcf,['Results/PSV/output3/SV_',num2str(l,'%4.4i') '.eps'],'psc2');
end



%% plot seismogram
figure
temp1=load('./OUTPUT_FILES/seismograms_01_1');
subplot(3,1,1)
plot(temp1(:,1),temp1(:,2));
ylabel('amplitude');
%xlabel('time')
title('P');
text(18,0.05,'P');text(23,0.05,'SV');text(26,0.05,'Rayleigh')

temp1=load('./OUTPUT_FILES/seismograms_01_2');
subplot(3,1,2)
plot(temp1(:,1),temp1(:,2));
ylabel('amplitude'); 
%xlabel('time')
title('SH ');
%text(26,0.07,'S');text(32,0.07,'Love');

temp1=load('./OUTPUT_FILES/seismograms_01_3');
subplot(3,1,3)
plot(temp1(:,1),temp1(:,2));
ylabel('amplitude'); 
xlabel('time')
title('SV');
text(18,0.016,'P');text(23,0.016,'SV');text(26,0.016,'Rayleigh')
saveas(gcf,'Results/PSV/output3/seismograms.eps','psc2');



