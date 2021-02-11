%%
clc
close all
figure('Position',[40,40,1250,1250],'Visible','on');
set(gcf,'color','white')
subplot(231);
set(gca,'units','normalized','position',[0.0125 0.58 0.35 0.35])

colormap jet
ax = usamap([nStartLat,nEndLat+0.52],[nStartLon,nEndLon]);
setm(ax,'PLabelMeridian','west','MLabelParallel','south');%ºá×Ý×ø±êÖáÎ»ÖÃ
geoshow('landareas.shp', 'FaceColor',[220,220,220]/255,'edgecolor',[211,211,211]/255);
hold on
% scatterm(D703_d(:,7),D703_d(:,6),25,(D703_d(:,8)-Ne703_d')./Ne703_d','filled')
% scatterm(D702_d(:,7),D702_d(:,6),25,(D702_d(:,8)-Ne702_d')./Ne702_d','filled')
% scatterm(D701_d(:,7),D701_d(:,6),25,(D701_d(:,8)-Ne701_d')./Ne701_d','filled')
% scatterm(D630_d(:,7),D630_d(:,6),25,(D630_d(:,8)-Ne630_d')./Ne630_d','filled')
% scatterm(D629_d(:,7),D629_d(:,6),25,(D629_d(:,8)-Ne629_d')./Ne629_d','filled')

scatterm(D707_d(:,7),D707_d(:,6),25,(D707_d(:,8)-Ne707_d')./Ne707_d','filled')
scatterm(D704_d(:,7),D704_d(:,6),25,(D704_d(:,8)-Ne704_d')./Ne704_d','filled')
scatterm(D705_d(:,7),D705_d(:,6),25,(D705_d(:,8)-Ne705_d')./Ne705_d','filled')
scatterm(D701_d(:,7),D701_d(:,6),25,(D701_d(:,8)-Ne701_d')./Ne701_d','filled')
scatterm(D708_d(:,7),D708_d(:,6),25,(D708_d(:,8)-Ne708_d')./Ne708_d','filled')

% scatterm(D715_d(:,7),D715_d(:,6),25,(D715_d(:,8)-Ne715_d')./Ne715_d','filled')
% scatterm(D716_d(:,7),D716_d(:,6),25,(D716_d(:,8)-Ne716_d')./Ne716_d','filled')
% scatterm(D717_d(:,7),D717_d(:,6),25,(D717_d(:,8)-Ne717_d')./Ne717_d','filled')

% scatterm(D713_d(:,7),D713_d(:,6),25,(D713_d(:,8)-Ne713_d')./Ne713_d','filled')
% scatterm(D719_d(:,7),D719_d(:,6),25,(D719_d(:,8)-Ne719_d')./Ne719_d','filled')

% Add text
textm(8,119,'07/05 06:50 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')
textm(-20,113.5,'07:00 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')


textm(8,124,'07/06 05:31 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')
textm(-20,118.5,'05:41 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')


textm(8,129,'07/07 05:23 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')
textm(-20,123.5,'05:12 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')

textm(8,139,'07/04 04:34 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')
textm(-20,133.5,'04:45 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')

 textm(8,134,'07/08 05:53 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')
 textm(-20,128,'06:04 UT','rotation',80,'fontweight','bold','fontsize',9,...
     'color','black','fontname','times new roman')

textm(17,110,'a','fontweight','bold','fontsize',12)

plotm(Eq_lat,Eq_lon,'rp','MarkerFaceColor','r','LineWidth',2 )
set(gca,'Clim',[-1.5,1.5])
cb = colorbar('vertical','Ticks',[-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5],'ticklabels',{'-200%','-150%','-100%','-50%','0','50%','100%','150%','200%','250%'});
% cb=colorbar('vertical');
circlem(Eq_lat,Eq_lon,1247.5,'edgecolor','red','linestyle','--','linewidth',1)
title('Deviation of Electron Density','fontname','times new roman','fontsize',12)
cb.Title.String='dNe'
ax=gca;
ax.FontName='Times New Roman';
ax.FontSize=11;
ax.FontWeight='bold';
%%
subplot(232)
set(gca,'units','normalized','position',[0.33 0.58 0.35 0.35])
ax = usamap([nStartLat,nEndLat+0.52],[nStartLon,nEndLon]);
setm(ax,'PLabelMeridian','west','MLabelParallel','south');%ºá×Ý×ø±êÖáÎ»ÖÃ
geoshow('landareas.shp', 'FaceColor',[220,220,220]/255,'edgecolor',[211,211,211]/255);
hold on
% scatterm(D703_d(:,7),D703_d(:,6),25,(D703_d(:,9)-Te703_d')./Te703_d','filled')
% scatterm(D702_d(:,7),D702_d(:,6),25,(D702_d(:,9)-Te702_d')./Te702_d','filled')
% scatterm(D701_d(:,7),D701_d(:,6),25,(D701_d(:,9)-Te701_d')./Te701_d','filled')
% scatterm(D630_d(:,7),D630_d(:,6),25,(D630_d(:,9)-Te630_d')./Te630_d','filled')
% scatterm(D629_d(:,7),D629_d(:,6),25,(D629_d(:,9)-Te629_d')./Te629_d','filled')

scatterm(D707_d(:,7),D707_d(:,6),25,(D707_d(:,9)-Te707_d')./Te707_d','filled')
scatterm(D706_d(:,7),D706_d(:,6),25,(D706_d(:,9)-Te706_d')./Te706_d','filled')
scatterm(D705_d(:,7),D705_d(:,6),25,(D705_d(:,9)-Te705_d')./Te705_d','filled')
scatterm(D704_d(:,7),D704_d(:,6),25,(D704_d(:,9)-Te704_d')./Te704_d','filled')
scatterm(D708_d(:,7),D708_d(:,6),25,(D708_d(:,9)-Te708_d')./Te708_d','filled')

% scatterm(D717_d(:,7),D717_d(:,6),25,(D717_d(:,8)-Ne717_d')./Ne717_d','filled')
% scatterm(D715_d(:,7),D715_d(:,6),25,(D715_d(:,9)-Te715_d')./Te715_d','filled')
% scatterm(D716_d(:,7),D716_d(:,6),25,(D716_d(:,9)-Te716_d')./Te716_d','filled')
% scatterm(D717_d(:,7),D717_d(:,6),25,(D717_d(:,9)-Te717_d')./Te717_d','filled')
% scatterm(D718_d(:,7),D718_d(:,6),25,(D718_d(:,9)-Te718_d')./Te718_d','filled')
% scatterm(D719_d(:,7),D719_d(:,6),25,(D719_d(:,9)-Te719_d')./Te719_d','filled')
% Add text
textm(8,119,'07/05 06:50 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')
textm(-20,113.5,'07:00 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')


textm(8,124,'07/06 05:31 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')
textm(-20,118.5,'05:41 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')


textm(8,129,'07/07 05:23 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')
textm(-20,123.5,'05:12 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')

textm(8,139,'07/04 04:34 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')
textm(-20,133.5,'04:45 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')

 textm(8,134,'07/08 05:53 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')
 textm(-20,128,'06:04 UT','rotation',80,'fontweight','bold','fontsize',9,...
     'color','black','fontname','times new roman')

textm(17,110,'b','fontweight','bold','fontsize',12)

plotm(Eq_lat,Eq_lon,'rp','MarkerFaceColor','r','LineWidth',2 )
set(gca,'Clim',[-1.5,1.5])
cb = colorbar('vertical','Ticks',[-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5],'ticklabels',{'-150%','-100%','-50%','0','50%','100%','150%','200%','250%','300%','350%'});
%  cb=colorbar('vertical');
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
set(gca,'units','normalized','position',[0.65 0.58 0.35 0.35]);

colormap jet
ax = usamap([nStartLat,nEndLat+0.52],[nStartLon,nEndLon]);
setm(ax,'PLabelMeridian','west','MLabelParallel','south');%ºá×Ý×ø±êÖáÎ»ÖÃ
geoshow('landareas.shp', 'FaceColor',[220,220,220]/255,'edgecolor',[211,211,211]/255);
hold on
% scatterm(pap_D701_d(:,11),pap_D701_d(:,10),25,(pap_D701_d(:,14)-O701_d')./O701_d','filled')
% scatterm(pap_D702_d(:,11),pap_D702_d(:,10),25,(pap_D702_d(:,14)-O702_d')./O702_d','filled')
% scatterm(pap_D703_d(:,11),pap_D703_d(:,10),25,(pap_D703_d(:,14)-O703_d')./O703_d','filled')
% scatterm(pap_D630_d(:,11),pap_D630_d(:,10),25,(pap_D630_d(:,14)-O630_d')./O630_d','filled')
% scatterm(pap_D629_d(:,11),pap_D629_d(:,10),25,(pap_D629_d(:,14)-O629_d')./O629_d','filled')

scatterm(pap_D704_d(:,11),pap_D704_d(:,10),25,(pap_D704_d(:,14)-O704_d')./O704_d','filled')
scatterm(pap_D705_d(:,11),pap_D705_d(:,10),25,(pap_D705_d(:,14)-O705_d')./O705_d','filled')
scatterm(pap_D706_d(:,11),pap_D706_d(:,10),25,(pap_D706_d(:,14)-O706_d')./O706_d','filled')
scatterm(pap_D707_d(:,11),pap_D707_d(:,10),25,(pap_D707_d(:,14)-O707_d')./O707_d','filled')
scatterm(pap_D708_d(:,11),pap_D708_d(:,10),25,(pap_D708_d(:,14)-O708_d')./O708_d','filled')

% scatterm(pap_D713_d(:,11),pap_D713_d(:,10),25,(pap_D713_d(:,14)-O713_d')./O713_d','filled')
% scatterm(pap_D715_d(:,11),pap_D715_d(:,10),25,(pap_D715_d(:,14)-O715_d')./O715_d','filled')
% scatterm(pap_D716_d(:,11),pap_D716_d(:,10),25,(pap_D716_d(:,14)-O716_d')./O716_d','filled')
% scatterm(pap_D717_d(:,11),pap_D717_d(:,10),25,(pap_D717_d(:,14)-O717_d')./O717_d','filled')
% scatterm(pap_D713_d(:,11),pap_D713_d(:,10),25,(pap_D713_d(:,14)-O713_d')./O713_d','filled')
% scatterm(pap_D719_d(:,11),pap_D719_d(:,10),25,(pap_D719_d(:,14)-O719_d')./O719_d','filled')

plotm(Eq_lat,Eq_lon,'rp','MarkerFaceColor','r','LineWidth',2)
set(gca,'Clim',[-1.5,1.5])
cb = colorbar('vertical','Ticks',[-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5],'ticklabels',{'-150%','-100%','-50%','0','50%','100%','150%','200%','250%','300%','350%'});
% cb=colorbar('vertical');
circlem(Eq_lat,Eq_lon,1247.5,'edgecolor','red','linestyle','--','linewidth',1)
title('Deviation of O^{+} Density ','fontname','times new roman','fontsize',12)
% Add text
textm(8,119,'07/05 06:50 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')
textm(-20,113.5,'07:00 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')


textm(8,124,'07/06 05:31 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')
textm(-20,118.5,'05:41 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')


textm(8,129,'07/07 05:23 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')
textm(-20,123.5,'05:12 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')

textm(8,139,'07/04 04:34 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')
textm(-20,133.5,'04:45 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')

 textm(8,134,'07/08 05:53 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')
 textm(-20,128,'06:04 UT','rotation',80,'fontweight','bold','fontsize',9,...
     'color','black','fontname','times new roman')
textm(17,110,'c','fontweight','bold','fontsize',12)
% c.Label.String='Deviation of Electron Density';
cb.Title.String='dN_{O^{+}}';
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


% scatterm(pap_D703_d(:,11),pap_D703_d(:,10),25,(pap_D703_d(:,15)-Ti703_d')./Ti703_d','filled')
% scatterm(pap_D702_d(:,11),pap_D702_d(:,10),25,(pap_D702_d(:,15)-Ti702_d')./Ti702_d','filled')
% scatterm(pap_D701_d(:,11),pap_D701_d(:,10),25,(pap_D701_d(:,15)-Ti701_d')./Ti701_d','filled')
% scatterm(pap_D630_d(:,11),pap_D630_d(:,10),25,(pap_D630_d(:,15)-Ti630_d')./Ti630_d','filled')
% scatterm(pap_D629_d(:,11),pap_D629_d(:,10),25,(pap_D629_d(:,15)-Ti629_d')./Ti629_d','filled')

scatterm(pap_D704_d(:,11),pap_D704_d(:,10),25,(pap_D704_d(:,15)-Ti704_d')./Ti704_d','filled')
scatterm(pap_D705_d(:,11),pap_D705_d(:,10),25,(pap_D705_d(:,15)-Ti705_d')./Ti705_d','filled')
scatterm(pap_D706_d(:,11),pap_D706_d(:,10),25,(pap_D706_d(:,15)-Ti706_d')./Ti706_d','filled')
scatterm(pap_D707_d(:,11),pap_D707_d(:,10),25,(pap_D707_d(:,15)-Ti707_d')./Ti707_d','filled')
scatterm(pap_D708_d(:,11),pap_D708_d(:,10),25,(pap_D708_d(:,15)-Ti708_d')./Ti708_d','filled')

% scatterm(pap_D715_d(:,11),pap_D715_d(:,10),25,(pap_D715_d(:,15)-Ti715_d')./Ti715_d','filled')
% scatterm(pap_D716_d(:,11),pap_D716_d(:,10),25,(pap_D716_d(:,15)-Ti716_d')./Ti716_d','filled')
% scatterm(pap_D717_d(:,11),pap_D717_d(:,10),25,(pap_D717_d(:,15)-Ti717_d')./Ti717_d','filled')
% scatterm(pap_D718_d(:,11),pap_D718_d(:,10),25,(pap_D718_d(:,15)-Ti713_d')./Ti718_d','filled')
% scatterm(pap_D719_d(:,11),pap_D719_d(:,10),25,(pap_D719_d(:,15)-Ti719_d')./Ti719_d','filled')

% Add text
textm(8,119,'07/05 06:50 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')
textm(-20,113.5,'07:00 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')


textm(8,124,'07/06 05:31 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')
textm(-20,118.5,'05:41 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')


textm(8,129,'07/07 05:23 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')
textm(-20,123.5,'05:12 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')

textm(8,139,'07/04 04:34 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')
textm(-20,133.5,'04:45 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')

 textm(8,134,'07/08 05:53 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')
 textm(-20,128,'06:04 UT','rotation',80,'fontweight','bold','fontsize',9,...
     'color','black','fontname','times new roman')
 
textm(17,110,'d','fontweight','bold','fontsize',12)
plotm(Eq_lat,Eq_lon,'rp','MarkerFaceColor','r','LineWidth',2 )

set(gca,'Clim',[-1.5,1.5])
cb = colorbar('vertical','Ticks',[-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5],'ticklabels',{'-150%','-100%','-50%','0','50%','100%','150%','200%','250%','300%','350%'});
circlem(Eq_lat,Eq_lon,1247.5,'edgecolor','red','linestyle','--','linewidth',1)
title('Deviation of Ion Temperature','fontname','times new roman','fontsize',12)
cb.Title.String='dTi';
ax=gca;
ax.FontName='Times New Roman';
ax.FontSize=11;
ax.FontWeight='bold';
%%
subplot(236);
set(gca,'units','normalized','position',[0.5 0.14 0.35 0.35])
colormap jet
ax = usamap([nStartLat,nEndLat+0.52],[nStartLon,nEndLon]);
setm(ax,'PLabelMeridian','west','MLabelParallel','south');%ºá×Ý×ø±êÖáÎ»ÖÃ
geoshow('landareas.shp', 'FaceColor',[220,220,220]/255,'edgecolor',[211,211,211]/255);
hold on
% scatterm(He_0701D(:,11),He_0701D(:,10),25,(He_0701D(:,13)-He701_d')./He701_d','filled')
% scatterm(He_0702D(:,11),He_0702D(:,10),25,(He_0702D(:,13)-He702_d')./He702_d','filled')
% scatterm(He_0703D(:,11),He_0703D(:,10),25,(He_0703D(:,13)-He703_d')./He703_d','filled')
% scatterm(He_0630D(:,11),He_0630D(:,10),25,(He_0630D(:,13)-He630_d')./He630_d','filled')
% scatterm(He_0629D(:,11),He_0629D(:,10),25,(He_0629D(:,13)-He629_d')./He629_d','filled')

 scatterm(He_0706D(:,11),He_0706D(:,10),25,(He_0706D(:,13)-He706_d')./He706_d','filled')
 scatterm(He_0705D(:,11),He_0705D(:,10),25,(He_0705D(:,13)-He705_d')./He705_d','filled')
 scatterm(He_0704D(:,11),He_0704D(:,10),25,(He_0704D(:,13)-He704_d')./He704_d','filled')
scatterm(He_0707D(:,11),He_0707D(:,10),25,(He_0707D(:,13)-He707_d')./He707_d','filled')
scatterm(He_0708D(:,11),He_0708D(:,10),25,(He_0708D(:,13)-He708_d')./He708_d','filled')
% scatterm(pap_D715_d(:,11),pap_D715_d(:,10),25,(pap_D715_d(:,13)-He715_d')./He715_d','filled')

% scatterm(He_0716D(:,11),He_0716D(:,10),25,(He_0716D(:,13)-He716_d')./He716_d','filled')
% scatterm(He_0715D(:,11),He_0715D(:,10),25,(He_0715D(:,13)-He715_d')./He715_d','filled')
% scatterm(He_0717D(:,11),He_0717D(:,10),25,(He_0717D(:,13)-He717_d')./He717_d','filled')
% scatterm(He_0713D(:,11),He_0713D(:,10),25,(He_0713D(:,13)-He713_d')./He713_d','filled')
% scatterm(He_0719D(:,11),He_0719D(:,10),25,(He_0719D(:,13)-He719_d')./He719_d','filled')
% Add text
textm(8,119,'07/05 06:50 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')
textm(-20,113.5,'07:00 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')


textm(8,124,'07/06 05:31 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')
textm(-20,118.5,'05:41 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')


textm(8,129,'07/07 05:23 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')
textm(-20,123.5,'05:12 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')

textm(8,139,'07/04 04:34 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')
textm(-20,133.5,'04:45 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')

 textm(8,134,'07/08 05:53 UT','rotation',80,'fontweight','bold','fontsize',9,...
    'color','black','fontname','times new roman')
 textm(-20,128,'06:04 UT','rotation',80,'fontweight','bold','fontsize',9,...
     'color','black','fontname','times new roman')

textm(17,110,'e','fontweight','bold','fontsize',12)

plotm(Eq_lat,Eq_lon,'rp','MarkerFaceColor','r','LineWidth',2 )
 set(gca,'Clim',[-1.5,1.5])
cb = colorbar('vertical','Ticks',[-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5],'ticklabels',{'-150%','-100%','-50%','0','50%','100%','150%','200%','250%','300%','350%'});
%  cb=colorbar('vertical');
circlem(Eq_lat,Eq_lon,1247.5,'edgecolor','red','linestyle','--','linewidth',1)
title('Deviation of He^{+} Density','fontname','times new roman','fontsize',12)
cb.Title.String='dN_{He^{+}}';
ax=gca;
ax.FontName='Times New Roman';
ax.FontSize=11;
ax.FontWeight='bold';