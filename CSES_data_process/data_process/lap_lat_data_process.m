%%  Read the data(.txt)
% filepath = 'F:\张衡一号LAP数据\LAP_txt';
% filepath = 'E:\MATLAB\Mycode\Reserch\NNE_of_Laiwui_EQ\LAP_TXT';

%filepath = 'F:\CSES\send\lapwithtrap';
filepath = 'F:\CSES\send\LAPwithtime\LAPwithtime';
data_list = dir([filepath,'\*.txt']); % READ the file(.txt)
num_file = length(data_list);
data = [];

for k = 1: num_file
    %location = h5read([filepath filename(k).name],'/LonLat');
    data1 = importdata(data_list(k).name);
% Data selection
Eq_lat = -0.52;
Eq_lon =128.17 ;

    data = [data;data1]; 
end
%%
nStartLat = Eq_lat-60;
nEndLat = Eq_lat+60;
nStartLon = Eq_lon-18;
nEndLon = Eq_lon+15;


%% Selecting the data of JUNE 14th-30th
orbit_Select_data1 = []; % JUNE 14th-30th
orbit_selection=1; % 0 for descending; 1 for ascending

% Data Selection Process
for j = 1:length(data)
    
if ( (data(j,6) >= nStartLon) && (data(j,6) <= nEndLon)) % data(j,2)==0)

      if data(j,4) == 6 && data(j,5) >= 25
            orbit_Select_data1 = [orbit_Select_data1;data(j,:)];
      end
     
end

end
%% Selecting the data of JULY 1st- 14th
orbit_Select_data2 = []; % JULY 1st- 14th
for j = 1:length(data)
    
if ((data(j,6) >= nStartLon) && (data(j,6) <= nEndLon))% && data(j,2)==0)
   
    
      if data(j,4) == 7 && data(j,5) <= 22
            orbit_Select_data2=[orbit_Select_data2;data(j,:)];
      end
end

end
%% The data 30 days before the EQ
% Adding Select_data1 and Select_data2 5
orbit_Select1_data = [orbit_Select_data1;orbit_Select_data2];
%%  Delete the data with error (-99999)
[c,~] = size(orbit_Select1_data);
orbit_Select_data = []; % Selected Data with Descending orbits
% Select_data_a=[]; % Selected Data with Ascending orbits

for i = 1:cs
    
    if orbit_Select1_data(i,8) ~= -99999
        orbit_Select_data = [orbit_Select_data;orbit_Select1_data(i,:)];      
          
    end  
    
end
%% orbit_overview
%f1 = figure(1)
figure('Position',[40,40,900,600],'Visible','on');
colormap Jet
ax = usamap([nStartLat,nEndLat],[nStartLon,nEndLon]);
setm(ax,'PLabelMeridian','west','MLabelParallel','south');%横纵坐标轴位置
geoshow('landareas.shp', 'FaceColor', 'none');
 hold on
 plotm(Eq_lat,Eq_lon,'rp','MarkerFaceColor','r','LineWidth',2 )
 plotm( orbit_Select_data(:,7), orbit_Select_data(:,6),'.b')
 % plotm( Select_data_a(:,7), Select_data_a(:,6),'.b')
 
 title('Descending orbits from 2019-06-14 to 2019-07-14','fontname','times new roman','fontsize',15) 
 % title('Ascending orbits from 2019-06-14 to 2019-07-14','fontname','times new roman','fontsize',15)
 