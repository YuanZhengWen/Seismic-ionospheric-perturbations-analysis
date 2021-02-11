%%
close all
%%figure('Position',[40,40,1250,1100],'Visible','on');
% figure('units','normalized','position',[ 0.0609    0.0593    0.5844    0.8463]);
figure('Position',[40,40,1250,1250],'Visible','on');
set(gcf,'color','white')
subplot(231);
set(gca,'units','normalized','position',[0.0125 0.58 0.35 0.35])

colormap jet
ax = usamap([nStartLat,nEndLat+0.52],[nStartLon,nEndLon]);
setm(ax,'PLabelMeridian','west','MLabelParallel','south');%ºá×Ý×ø±êÖáÎ»ÖÃ
geoshow('landareas.shp', 'FaceColor',[220,220,220]/255,'edgecolor',[211,211,211]/255);
hold on

scatterm(D710_d(:,7),D710_d(:,6),25,(D710_d(:,8)-Ne710_d')./Ne710_d','filled')
scatterm(D711_d(:,7),D711_d(:,6),25,(D711_d(:,8)-Ne711_d')./Ne711_d','filled')
scatterm(D712_d(:,7),D712_d(:,6),25,(D712_d(:,8)-Ne712_d')./Ne712_d','filled')
scatterm(D713_d(:,7),D713_d(:,6),25,(D713_d(:,8)-Ne713_d')./Ne713_d','filled')
scatterm(D714_d(:,7),D714_d(:,6),25,(D714_d(:,8)-Ne714_d')./Ne714_d','filled')

% Add text
% textm(-10.5,115,'0710','fontweight','bold','fontsize',10,'fontname','times new roman','color','black');
%    textm(14,124,'06:50','rotation',80,'fontweight','bold','fontsize',9,...
%       'color','black','fontname','times new roman')
%    textm(-20,117.5,'07:01','rotation',80,'fontweight','bold','fontsize',9,...
%        'color','black','fontname','times new roman')
%     textm(14,129,'05:31','rotation',80,'fontweight','bold','fontsize',9,...
%         'color','black','fontname','times new roman')
%   textm(-20,122.5,'05:42','rotation',80,'fontweight','bold','fontsize',9,...
%       'color','black','fontname','times new roman');
%  textm(-10.5,120,'0711','fontweight','bold','fontsize',10,'fontname','times new roman','color','black');
%     textm(14,134,'05:12','rotation',80,'fontweight','bold','fontsize',9,...
%         'color','black','fontname','times new roman');
%     textm(-20,127.5,'05:23','rotation',80,'fontweight','bold','fontsize',9,...
%        'color','black','fontname','times new roman');
%  textm(-10.5,125,'0712','fontweight','bold','fontsize',10,'fontname','times new roman','color','black');
%   textm(-10.5,130,'0713','fontweight','bold','fontsize',10,'fontname','times new roman','color','black');
%    textm(6.8,134,'July 13th 05:53','rotation',80,'fontweight','bold','fontsize',10,...
%        'color','black','fontname','times new roman')
%     textm(-20,128,'06:04','rotation',80,'fontweight','bold','fontsize',10,...
%        'color','black','fontname','times new roman')
% textm(18.5,118,'0710','fontweight','bold','fontsize',10,'fontname','times new roman','color','black');
%  textm(-19.5,108.5,'1450-1501','fontweight','bold','fontsize',8.5,'fontname','times new roman','color','black');
%  textm(18.5,124,'0711','fontweight','bold','fontsize',10,'fontname','times new roman','color','black');
%  textm(-19.5,116,'1331-1342','fontweight','bold','fontsize',8.5,'fontname','times new roman','color','black');
%  textm(18.5,130,'0712','fontweight','bold','fontsize',10,'fontname','times new roman','color','black');
%  textm(18.5,136,'0713','fontweight','bold','fontsize',10,'fontname','times new roman','color','black');

 textm(8,119,'07/10 06:50 UT','rotation',80,'fontweight','bold','fontsize',9,...
     'color','black','fontname','times new roman')
 textm(-20,113.5,'07:01 UT','rotation',80,'fontweight','bold','fontsize',9,...
     'color','black','fontname','times new roman')

 textm(8,124,'07/11 05:31 UT','rotation',80,'fontweight','bold','fontsize',9,...
     'color','black','fontname','times new roman')
 textm(-20,118.5,'05:42 UT','rotation',80,'fontweight','bold','fontsize',9,...
     'color','black','fontname','times new roman')

 textm(8,129,'07/12 05:12 UT','rotation',80,'fontweight','bold','fontsize',9,...
     'color','black','fontname','times new roman')
 textm(-20,123.5,'05:23 UT','rotation',80,'fontweight','bold','fontsize',9,...
     'color','black','fontname','times new roman')

 textm(8,134,'07/13 05:53 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')
 textm(-20,128,'06:04 UT','rotation',80,'fontweight','bold','fontsize',9,...
     'color','black','fontname','times new roman')
 
 textm(8,139,'07/09 04:34 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')
textm(-20,133.5,'04:45 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')

textm(17,110,'a','fontweight','bold','fontsize',12)

plotm(Eq_lat,Eq_lon,'rp','MarkerFaceColor','r','LineWidth',2 )
set(gca,'Clim',[-1.5,1.5])
cb = colorbar('vertical','Ticks',[-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5],'ticklabels',{'-200%','-150%','-100%','-50%','0','50%','100%','150%','200%','250%'});
% cb=colorbar('vertical');
circlem(Eq_lat,Eq_lon,1247.5,'edgecolor','red','linestyle','--','linewidth',1)
title('Deviation of Electron Density','fontname','times new roman','fontsize',12)
% title('July 10-13, 2019');
cb.Title.String='dNe';
ax=gca;
ax.FontName='Times New Roman';
ax.FontSize=11;
ax.FontWeight='bold';
%%
subplot(232)
set(gca,'units','normalized','position',[0.33 0.58 0.35 0.35])
% set(gca,'units','normalized','position',[0.13 0.11 0.38 0.38]);

ax = usamap([nStartLat,nEndLat+0.52],[nStartLon,nEndLon]);
setm(ax,'PLabelMeridian','west','MLabelParallel','south');%ºá×Ý×ø±êÖáÎ»ÖÃ
geoshow('landareas.shp', 'FaceColor',[220,220,220]/255,'edgecolor',[211,211,211]/255);
hold on
scatterm(D710_d(:,7),D710_d(:,6),25,(D710_d(:,9)-Te710_d')./Te710_d','filled')
scatterm(D711_d(:,7),D711_d(:,6),25,(D711_d(:,9)-Te711_d')./Te711_d','filled')
scatterm(D712_d(:,7),D712_d(:,6),25,(D712_d(:,9)-Te712_d')./Te712_d','filled')
scatterm(D713_d(:,7),D713_d(:,6),25,(D713_d(:,9)-Te713_d')./Te713_d','filled')
scatterm(D714_d(:,7),D714_d(:,6),25,(D714_d(:,9)-Te714_d')./Te714_d','filled')

%Add text
 textm(8,119,'07/10 06:50 UT','rotation',80,'fontweight','bold','fontsize',9,...
     'color','black','fontname','times new roman')
 textm(-20,113.5,'07:01 UT','rotation',80,'fontweight','bold','fontsize',9,...
     'color','black','fontname','times new roman')

 textm(8,124,'07/11 05:31 UT','rotation',80,'fontweight','bold','fontsize',9,...
     'color','black','fontname','times new roman')
 textm(-20,118.5,'05:42 UT','rotation',80,'fontweight','bold','fontsize',9,...
     'color','black','fontname','times new roman')

 textm(8,129,'07/12 05:12 UT','rotation',80,'fontweight','bold','fontsize',9,...
     'color','black','fontname','times new roman')
 textm(-20,123.5,'05:23 UT','rotation',80,'fontweight','bold','fontsize',9,...
     'color','black','fontname','times new roman')

 textm(8,134,'07/13 05:53 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')
 textm(-20,128,'06:04 UT','rotation',80,'fontweight','bold','fontsize',9,...
     'color','black','fontname','times new roman')
 
 textm(8,139,'07/09 04:34 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')
textm(-20,133.5,'04:45 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')

textm(17,110,'b','fontweight','bold','fontsize',12)

plotm(Eq_lat,Eq_lon,'rp','MarkerFaceColor','r','LineWidth',2 )
set(gca,'Clim',[-1.5,1.5])
cb = colorbar('vertical','Ticks',[-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3],'ticklabels',{'-200%','-150%','-100%','-50%','0','50%','100%','150%','200%','250%','300%'});
% cb=colorbar('vertical');
circlem(Eq_lat,Eq_lon,1247.5,'edgecolor','red','linestyle','--','linewidth',1)
title('Deviation of Electron Temperature','fontname','times new roman','fontsize',12)
cb.Title.String='dTe'
ax=gca;
ax.FontName='Times New Roman';
ax.FontSize=11;
ax.FontWeight='bold';
%%
subplot(233);
% set(gca,'units','normalized','position',[0.13 0.583837209302326 0.38 0.38])
set(gca,'units','normalized','position',[0.65 0.58 0.35 0.35])

colormap jet
ax = usamap([nStartLat,nEndLat+0.52],[nStartLon,nEndLon]);
setm(ax,'PLabelMeridian','west','MLabelParallel','south');%ºá×Ý×ø±êÖáÎ»ÖÃ
geoshow('landareas.shp', 'FaceColor',[220,220,220]/255,'edgecolor',[211,211,211]/255);
hold on

scatterm(pap_D710_d(:,11),pap_D710_d(:,10),25,(pap_D710_d(:,14)-O710_d')./O710_d','filled')
scatterm(pap_D711_d(:,11),pap_D711_d(:,10),25,(pap_D711_d(:,14)-O711_d')./O711_d','filled')
scatterm(pap_D712_d(:,11),pap_D712_d(:,10),25,(pap_D712_d(:,14)-O712_d')./O712_d','filled')
scatterm(pap_D713_d(:,11),pap_D713_d(:,10),25,(pap_D713_d(:,14)-O713_d')./O713_d','filled')
scatterm(pap_D714_d(:,11),pap_D714_d(:,10),25,(pap_D714_d(:,14)-O714_d')./O714_d','filled')
% Add text
 textm(8,119,'07/10 06:50 UT','rotation',80,'fontweight','bold','fontsize',9,...
     'color','black','fontname','times new roman')
 textm(-20,113.5,'07:01 UT','rotation',80,'fontweight','bold','fontsize',9,...
     'color','black','fontname','times new roman')

 textm(8,124,'07/11 05:31 UT','rotation',80,'fontweight','bold','fontsize',9,...
     'color','black','fontname','times new roman')
 textm(-20,118.5,'05:42 UT','rotation',80,'fontweight','bold','fontsize',9,...
     'color','black','fontname','times new roman')

 textm(8,129,'07/12 05:12 UT','rotation',80,'fontweight','bold','fontsize',9,...
     'color','black','fontname','times new roman')
 textm(-20,123.5,'05:23 UT','rotation',80,'fontweight','bold','fontsize',9,...
     'color','black','fontname','times new roman')

 textm(8,134,'07/13 05:53 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')
 textm(-20,128,'06:04 UT','rotation',80,'fontweight','bold','fontsize',9,...
     'color','black','fontname','times new roman')
 
 textm(8,139,'07/09 04:34 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')
textm(-20,133.5,'04:45 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')

textm(17,110,'c','fontweight','bold','fontsize',12)

plotm(Eq_lat,Eq_lon,'rp','MarkerFaceColor','r','LineWidth',2)
set(gca,'Clim',[-1.5,1.5])
cb = colorbar('vertical','Ticks',[-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5],'ticklabels',{'-200%','-150%','-100%','-50%','0','50%','100%','150%','200%','250%','300%','350%'});
% cb=colorbar('vertical');
circlem(Eq_lat,Eq_lon,1247.5,'edgecolor','red','linestyle','--','linewidth',1)
title('Deviation of O^{+} Density ','fontname','times new roman','fontsize',12)
% c.Label.String='Deviation of Electron Density';
cb.Title.String='dN_{O^{+}}'
ax=gca;
ax.FontName='Times New Roman';
ax.FontSize=11;
ax.FontWeight='bold';
%%
subplot(234)
set(gca,'units','normalized','position',[0.15 0.14 0.35 0.35]);

ax = usamap([nStartLat,nEndLat+0.52],[nStartLon,nEndLon]);
setm(ax,'PLabelMeridian','west','MLabelParallel','south');%ºá×Ý×ø±êÖáÎ»ÖÃ
geoshow('landareas.shp', 'FaceColor',[220,220,220]/255,'edgecolor',[211,211,211]/255);
hold on

scatterm(pap_D710_d(:,11),pap_D710_d(:,10),25,(pap_D710_d(:,15)-Ti710_d')./Ti710_d','filled')
scatterm(pap_D711_d(:,11),pap_D711_d(:,10),25,(pap_D711_d(:,15)-Ti711_d')./Ti711_d','filled')
scatterm(pap_D712_d(:,11),pap_D712_d(:,10),25,(pap_D712_d(:,15)-Ti712_d')./Ti712_d','filled')
scatterm(pap_D713_d(:,11),pap_D713_d(:,10),25,(pap_D713_d(:,15)-Ti713_d')./Ti713_d','filled')
scatterm(pap_D714_d(:,11),pap_D714_d(:,10),25,(pap_D714_d(:,15)-Ti714_d')./Ti714_d','filled')

% Add text
 textm(8,119,'07/10 06:50 UT','rotation',80,'fontweight','bold','fontsize',9,...
     'color','black','fontname','times new roman')
 textm(-20,113.5,'07:01 UT','rotation',80,'fontweight','bold','fontsize',9,...
     'color','black','fontname','times new roman')

 textm(8,124,'07/11 05:31 UT','rotation',80,'fontweight','bold','fontsize',9,...
     'color','black','fontname','times new roman')
 textm(-20,118.5,'05:42 UT','rotation',80,'fontweight','bold','fontsize',9,...
     'color','black','fontname','times new roman')

 textm(8,129,'07/12 05:12 UT','rotation',80,'fontweight','bold','fontsize',9,...
     'color','black','fontname','times new roman')
 textm(-20,123.5,'05:23 UT','rotation',80,'fontweight','bold','fontsize',9,...
     'color','black','fontname','times new roman')

 textm(8,134,'07/13 05:53 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')
 textm(-20,128,'06:04 UT','rotation',80,'fontweight','bold','fontsize',9,...
     'color','black','fontname','times new roman')
 
 textm(8,139,'07/09 04:34 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')
textm(-20,133.5,'04:45 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')

textm(17,110,'d','fontweight','bold','fontsize',12)

plotm(Eq_lat,Eq_lon,'rp','MarkerFaceColor','r','LineWidth',2 )
set(gca,'Clim',[-1.5,1.5])
cb = colorbar('vertical','Ticks',[-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3],'ticklabels',{'-200%','-150%','-100%','-50%','0','50%','100%','150%','200%','250%','300%'});
%cb=colorbar('vertical');
circlem(Eq_lat,Eq_lon,1247.5,'edgecolor','red','linestyle','--','linewidth',1)
title('Deviation of Ion Temperature','fontname','times new roman','fontsize',12)
ax=gca;
ax.FontName='Times New Roman';
cb.Title.String='dTi'
ax=gca;
ax.FontName='Times New Roman';
ax.FontSize=11;
ax.FontWeight='bold';
%%
subplot(236);
set(gca,'units','normalized','position',[0.5 0.14 0.35 0.35])
set(gcf,'color','white')
colormap jet
ax = usamap([nStartLat,nEndLat+0.52],[nStartLon,nEndLon]);
setm(ax,'PLabelMeridian','west','MLabelParallel','south');%ºá×Ý×ø±êÖáÎ»ÖÃ
geoshow('landareas.shp', 'FaceColor',[220,220,220]/255,'edgecolor',[211,211,211]/255);
hold on

scatterm(He_0710D(:,11),He_0710D(:,10),25,(He_0710D(:,13)-He710_d')./He710_d','filled')
scatterm(He_0711D(:,11),He_0711D(:,10),25,(He_0711D(:,13)-He711_d')./He711_d','filled')
scatterm(He_0712D(:,11),He_0712D(:,10),25,(He_0712D(:,13)-He712_d')./He712_d','filled')
scatterm(He_0713D(:,11),He_0713D(:,10),25,(He_0713D(:,13)-He713_d')./He713_d','filled')
scatterm(He_0714D(:,11),He_0714D(:,10),25,(He_0714D(:,13)-He714_d')./He714_d','filled')
% Add text
 textm(8,119,'07/10 06:50 UT','rotation',80,'fontweight','bold','fontsize',9,...
     'color','black','fontname','times new roman')
 textm(-20,113.5,'07:01 UT','rotation',80,'fontweight','bold','fontsize',9,...
     'color','black','fontname','times new roman')

 textm(8,124,'07/11 05:31 UT','rotation',80,'fontweight','bold','fontsize',9,...
     'color','black','fontname','times new roman')
 textm(-20,118.5,'05:42 UT','rotation',80,'fontweight','bold','fontsize',9,...
     'color','black','fontname','times new roman')

 textm(8,129,'07/12 05:12 UT','rotation',80,'fontweight','bold','fontsize',9,...
     'color','black','fontname','times new roman')
 textm(-20,123.5,'05:23 UT','rotation',80,'fontweight','bold','fontsize',9,...
     'color','black','fontname','times new roman')

 textm(8,134,'07/13 05:53 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')
 textm(-20,128,'06:04 UT','rotation',80,'fontweight','bold','fontsize',9,...
     'color','black','fontname','times new roman')
 
 textm(8,139,'07/09 04:34 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')
textm(-20,133.5,'04:45 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')


textm(17,110,'e','fontweight','bold','fontsize',12)


plotm(Eq_lat,Eq_lon,'rp','MarkerFaceColor','r','LineWidth',2 )
 set(gca,'Clim',[-1.5,1.5])
cb = colorbar('vertical','Ticks',[-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3],'ticklabels',{'-200%','-150%','-100%','-50%','0','50%','100%','150%','200%','250%','300%'});
% cb=colorbar('vertical');
circlem(Eq_lat,Eq_lon,1247.5,'edgecolor','red','linestyle','--','linewidth',1)
title('Deviation of He^{+} Density','fontname','times new roman','fontsize',12)
cb.Title.String='dN_{He^{+}}';
ax=gca;
ax.FontName='Times New Roman';
ax.FontSize=11;
ax.FontWeight='bold';