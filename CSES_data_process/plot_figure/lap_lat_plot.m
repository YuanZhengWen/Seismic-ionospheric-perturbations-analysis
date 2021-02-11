%%
[orbit_0713D_data,value_Ne_0713D_day,value_Te_0713D_day]=dailydata(7,13,0,orbit_Select_data);
[orbit_0713A_data,value__Ne_0713A_day,value_Te_0713A_day]=dailydata(7,13,1,orbit_Select_data);
[orbit_0712D_data,value_Ne_0712D_day,value_Te_0712D_day]=dailydata(7,12,0,orbit_Select_data);
[orbit_0712A_data,value_Ne_0712A_day,value_Te_0712A_day]=dailydata(7,12,1,orbit_Select_data);
[orbit_0711D_data,value_Ne_0711D_day,value_Te_0711D_day]=dailydata(7,11,0,orbit_Select_data);
[orbit_0711A_data,value_Ne_0711A_day,value_Te_0711A_day]=dailydata(7,11,1,orbit_Select_data);
[orbit_0710D_data,value_Ne_0710D_day,value_Te_0710D_day]=dailydata10(7,10,0,orbit_Select_data);
[orbit_0710A_data,value_Ne_0710A_day,value_Te_0710A_day]=dailydata10(7,10,1,orbit_Select_data);

[orbit_0717D_data,value_Ne_0717D_day,value_Te_0717D_day]=dailydata(7,17,0,orbit_Select_data);
[orbit_0707D_data,value_Ne_0707D_day,value_Te_0707D_day]=dailydata(7,7,0,orbit_Select_data);
[orbit_0702D_data,value_Ne_0702D_day,value_Te_0702D_day]=dailydata(7,2,0,orbit_Select_data);
[orbit_0630D_data,value_Ne_0630D_day,value_Te_0630D_day]=dailydata10(6,30,0,orbit_Select_data);
[orbit_0627D_data,value_Ne_0627D_day,value_Te_0627D_day]=dailydata10(6,27,0,orbit_Select_data);
% [orbit_0622D_data,value_Ne_0622D_day,value_Te_0627D_day]=dailydata(6,22,0,orbit_Select_data);

[orbit_0720D_data,value_Ne_0720D_day,value_Te_0720D_day]=dailydata10(7,20,0,orbit_Select_data);
[orbit_0722D_data,value_Ne_0722D_day,value_Te_0722D_day]=dailydata(7,22,0,orbit_Select_data);

[orbit_0715D_data,value_Ne_0715D_day,value_Te_0715D_day]=dailydata(7,15,0,orbit_Select_data);
[orbit_0705D_data,value_Ne_0705D_day,value_Te_0705D_day]=dailydata(7,5,0,orbit_Select_data);

%%
orbit_0713D_data=sortrows(orbit_0713D_data,7);
orbit_0713A_data=sortrows(orbit_0713A_data,7);
orbit_0712D_data=sortrows(orbit_0712D_data,7);
orbit_0712A_data=sortrows(orbit_0712A_data,7);
orbit_0711D_data=sortrows(orbit_0711D_data,7);
orbit_0711A_data=sortrows(orbit_0711A_data,7);
orbit_0710D_data=sortrows(orbit_0710D_data,7);
orbit_0710A_data=sortrows(orbit_0710A_data,7);
orbit_0627D_data=sortrows(orbit_0627D_data,7);
orbit_0717D_data=sortrows(orbit_0717D_data,7);
orbit_0707D_data=sortrows(orbit_0707D_data,7);
orbit_0702D_data1=sortrows(orbit_0702D_data,7);
a=length(orbit_0702D_data1);
orbit_0702D_data=[];
 for i=1:a 
     if orbit_0702D_data1(i,1)==78280
         orbit_0702D_data=[orbit_0702D_data;orbit_0702D_data1(i,:)];
     end
 end
 
orbit_0720D_data=sortrows(orbit_0720D_data,7);
orbit_0722D_data=sortrows(orbit_0722D_data,7);
% orbit_0627D_data=sortrows(orbit_0627D_data,7);
% orbit_0622D_data=sortrows(orbit_0622D_data,7);
orbit_0715D_data=sortrows(orbit_0715D_data,7);

orbit_0705D_data=sortrows(orbit_0705D_data,7);
orbit_0630D_data=sortrows(orbit_0630D_data,7);
% orbit_0625D_data=sortrows(orbit_0625D_data,7);
% orbit_0620D_data=sortrows(orbit_0620D_data,7);




% Delete wrong data from 0711A and 0712A
orbit_0712A_data(42,:)=[];
orbit_0712A_data(48,:)=[];
orbit_0711A_data(30,:)=[];


%%
figure('Position',[40,40,1000,750],'Visible','on');
set(gcf,'color','white')
subplot(221)
% plot(orbit_0713D_data(:,7),orbit_0713D_data(:,8),'linewidth',1.5)
plot(orbit_0712D_data(:,7),orbit_0712D_data(:,8),'linewidth',1.5)
hold on
plot(orbit_0707D_data(:,7),orbit_0707D_data(:,8),'linewidth',1.5)
plot(orbit_0627D_data(:,7),orbit_0627D_data(:,8),'linewidth',1.5)
plot(orbit_0717D_data(:,7),orbit_0717D_data(:,8),'linewidth',1.5)
% plot(orbit_0722D_data(:,7),orbit_0722D_data(:,8),'linewidth',1.5)
% plot(orbit_0622D_data(:,7),orbit_0622D_data(:,8),'linewidth',1.5)

% plot(orbit_0711D_data(:,7),orbit_0711D_data(:,8),'linewidth',1.5)
% plot(orbit_0710D_data(:,7),orbit_0710D_data(:,8),'linewidth',1.5)

set(gca,'xlim',[-60,60],'linewidth',1.5)
x=Eq_lat;
y=0;
plot(x,y,'rp','markersize',10,'markerfacecolor','red')
plot([7.6 7.6],[0,10E10],'k--','linewidth',1.8)
legend({'07/12','07/07','06/27','07/17','EQ','Magnetic Equator'},'fontweight','bold',...
    'location','northwest')
title('Electron Density','Fontname','Times New roman','fontsize',12,'fontweight','bold')
xlabel('Latitude/{\circ}','fontsize',12,'fontweight','bold','fontname','times new roman')
ylabel('Ne(/m^{3})','fontsize',12,'fontweight','bold','fontname','times new roman')
x=Eq_lat;
set(gca,'ylim',[0,10e10],'linewidth',1.5);
grid on
legend('boxoff')
set(gca,'xtick',-60:20:60,'xTickLabel',{'60^{\circ}S','40^{\circ}S','20^{\circ}S',...
    '0^{\circ}','20^{\circ}N','40^{\circ}N','60^{\circ}N'})
set(gca,'ytick',[0,2e10,4e10,6e10,8e10,10e10],'yTicklabel',{'0','2E10','4E10','6E10','8E10','10E10'});
ax=gca;
ax.FontName = 'Arial';
ax.FontSize = 11;
%%
subplot(222)
plot(orbit_0710D_data(:,7),orbit_0710D_data(:,8),'linewidth',1.5)
hold on
% plot(orbit_0712A_data(:,7),orbit_0712A_data(:,8),'linewidth',1.5)
% plot(orbit_0711A_data(:,7),orbit_0711A_data(:,8),'linewidth',1.5)
% plot(orbit_0710A_data(:,7),orbit_0710A_data(:,8),'linewidth',1.5)
plot(orbit_0715D_data(:,7),orbit_0715D_data(:,8),'linewidth',1.5)
plot(orbit_0705D_data(:,7),orbit_0705D_data(:,8),'linewidth',1.5)
plot(orbit_0630D_data(:,7),orbit_0630D_data(:,8),'linewidth',1.5)
plot(orbit_0720D_data(:,7),orbit_0720D_data(:,8),'linewidth',1.5)
x=Eq_lat;
y=0;
plot(x,y,'rp','markersize',10,'markerfacecolor','red')
plot([7.6 7.6],[0,10E10],'k--','linewidth',1.5)
set(gca,'xlim',[-60,60],'linewidth',1.5)
set(gca,'ytick',[0,2e10,4e10,6e10,8e10,10e10],'yTicklabel',{'0','2E10','4E10','6E10','8E10','10E10'});
% set(gca,'ylim',[0,3.5E10])

legend({'07/10','07/15','07/05','06/30','07/20','EQ','Magnetic Equator'},'fontweight','bold',...
    'location','northwest')
legend('boxoff')
xlabel('Latitude/{\circ}','fontsize',12,'fontweight','bold','fontname','times new roman')
ylabel('Ne(/m^{3})','fontsize',12,'fontweight','bold','fontname','times new roman')
title('Electron Density','Fontname','Times New roman','fontsize',12,'fontweight','bold')
grid on
set(gca,'xtick',-60:20:60,'xTickLabel',{'60^{\circ}S','40^{\circ}S','20^{\circ}S',...
    '0^{\circ}','20^{\circ}N','40^{\circ}N','60^{\circ}N'})
ax=gca;
ax.FontName = 'Arial';
ax.FontSize = 11;
%%
subplot(223)
plot(orbit_0712D_data(:,7),orbit_0712D_data(:,9),'linewidth',1.5)
hold on
plot(orbit_0707D_data(:,7),orbit_0707D_data(:,9),'linewidth',1.5)
plot(orbit_0627D_data(:,7),(orbit_0627D_data(:,9)),'linewidth',1.5)
plot(orbit_0717D_data(:,7),orbit_0717D_data(:,9),'linewidth',1.5)
% plot(orbit_0722D_data(:,7),orbit_0722D_data(:,9),'linewidth',1.5)
x=Eq_lat;
y=1000;
plot(x,y,'rp','markersize',10,'markerfacecolor','red')
plot([7.6 7.6],[1000,5000],'k--','linewidth',1.8)

set(gca,'xlim',[-60,60],'linewidth',1.5)
set(gca,'ylim',[1000,5000]);


legend({'07/12','07/07','06/27','07/17','EQ','Magnetic Equator'},'fontweight','bold',...
    'location','northwest')
title('Electron Temperature','Fontname','Times New roman','fontsize',12,'fontweight','bold')
grid on
xlabel('Latitude/{\circ}','fontsize',12,'fontweight','bold','fontname','times new roman')
ylabel('Temperature(K)','fontsize',12,'fontweight','bold','fontname','times new roman')
legend('boxoff')
set(gca,'xtick',-60:20:60,'xTickLabel',{'60^{\circ}S','40^{\circ}S','20^{\circ}S',...
    '0^{\circ}','20^{\circ}N','40^{\circ}N','60^{\circ}N'})
ax=gca;
ax.FontName = 'Arial';
ax.FontSize = 11;
%%
subplot(224)
plot(orbit_0710D_data(:,7),orbit_0710D_data(:,9),'linewidth',1.5)
hold on
plot(orbit_0715D_data(:,7),orbit_0715D_data(:,9),'linewidth',1.5)
plot(orbit_0705D_data(:,7),orbit_0705D_data(:,9),'linewidth',1.5)
plot(orbit_0630D_data(:,7),smooth(orbit_0630D_data(:,9),5),'linewidth',1.5)
plot(orbit_0720D_data(:,7),orbit_0720D_data(:,9),'linewidth',1.5)
plot(x,y,'rp','markersize',10,'markerfacecolor','red')
plot([7.6 7.6],[1000,5000],'k--','linewidth',1.8)
x=Eq_lat;
y=1000;

set(gca,'xlim',[-60,60],'linewidth',1.5)
% set(gca,'ylim',[1000,5000]);
legend({'07/10','07/15','07/05','06/30','07/20','EQ','Magnetic Equator'},'fontweight','bold',...
    'location','northwest')
title('Electron Temperature','Fontname','Times New roman','fontsize',12,'fontweight','bold')
xlabel('Latitude/{\circ}','fontsize',12,'fontweight','bold','fontname','times new roman')
ylabel('Temperature(K)','fontsize',12,'fontweight','bold','fontname','times new roman')
grid on
legend('boxoff')
set(gca,'xtick',-60:20:60,'xTickLabel',{'60^{\circ}S','40^{\circ}S','20^{\circ}S',...
    '0^{\circ}','20^{\circ}N','40^{\circ}N','60^{\circ}N'})
ax=gca;
ax.FontName = 'Arial';
ax.FontSize = 11;