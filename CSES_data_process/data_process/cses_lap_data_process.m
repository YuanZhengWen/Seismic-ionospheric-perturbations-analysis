%  Read the data(.txt)
filepath = 'F:\CSES\send\LAPwithtime\LAPwithtime';
data_list = dir([filepath,'\*.txt']); % READ the file(.txt)
num_file = length(data_list);
data = [];

for k = 1: num_file
    %location = h5read([filepath filename(k).name],'/LonLat');
    data1 = importdata(data_list(k).name);
    % Data selection
    Eq_lat = -0.52;
    Eq_lon =128.17;
    
    data = [data;data1];
end
% Covrage of the earthquake region
nStartLat = Eq_lat-20;
nEndLat = Eq_lat+20;
nStartLon = Eq_lon-20;
nEndLon = Eq_lon+20;


% Selecting the data of JUNE 14th-30th
Select_data1 = []; % JUNE 14th-30th
orbit_selection=1; % 0 for descending; 1 for ascending

% Data Selection Process
for j = 1:length(data)
    
    if ( (data(j,7) >= nStartLat) && (data(j,7) <= nEndLat) && (data(j,6) >= nStartLon) && (data(j,6) <= nEndLon)) % data(j,2)==0)
        
        if data(j,4) == 6 && data(j,5) >= 14
            Select_data1 = [Select_data1;data(j,:)];
        end
        
    end
    
end
% Selecting the data of JULY 1st- 14th
Select_data2 = []; % JULY 1st- 14th
for j = 1:length(data)
    
    if ((data(j,7) >= nStartLat) && (data(j,7) <= nEndLat) && (data(j,6) >= nStartLon) && (data(j,6) <= nEndLon))% && data(j,2)==0)
        
        
        if data(j,4) == 7 && data(j,5) <=14
            Select_data2=[Select_data2;data(j,:)];
        end
    end
    
end

Select1_data = [Select_data1;Select_data2]; % 30 days' data before the EQ
%%  Delete the data with error (-99999)
[c,~] = size(Select1_data);
Select_data = []; % Selected Data with Descending orbits
% Select_data_a=[]; % Selected Data with Ascending orbits

for i = 1:c
    
    if Select1_data(i,8) ~= -99999
        Select_data = [Select_data;Select1_data(i,:)];
        
    end
    
end
%% Get the daily data

% orbit observation data and grid data
[D706_d,Ne706d,Te706d]=dailydata(7,6,0,Select_data);
[D706_a,Ne706a,Te706a]=dailydata(7,6,1,Select_data);

[D703_d,Ne703d,Te703d]=dailydata(7,3,0,Select_data);
[D703_a,Ne703a,Te703a]=dailydata(7,3,1,Select_data);

[D701_d,Ne701d,Te701d]=dailydata(6,26,0,Select_data);
[D701_a,Ne701a,Te701a]=dailydata(6,26,1,Select_data);

[D630_d,Ne630d,Te630d]=dailydata10(6,30,0,Select_data);
[D630_a,Ne630a,Te630a]=dailydata10(6,30,1,Select_data);

[D629_d,Ne629d,Te629d]=dailydata19(6,29,0,Select_data);
[D629_a,Ne629a,Te629a]=dailydata19(6,29,1,Select_data);

[D709_d,Ne709d,Te709d]=dailydata19(7,9,0,Select_data);
[D709_a,Ne709a,Te709a]=dailydata19(7,9,1,Select_data);

[D708_d,Ne708d,Te708d]=dailydata10(7,8,0,Select_data);
[D708_a,Ne708a,Te708a]=dailydata10(7,8,1,Select_data);

[D710_d,Ne710d,Te710d]=dailydata10(7,10,0,Select_data);
[D710_a,Ne710a,Te710a]=dailydata10(7,10,1,Select_data);

[D711_d,Ne711d,Te711d]=dailydata(7,11,0,Select_data);
[D711_a,Ne711a,Te711a]=dailydata(7,11,1,Select_data);

[D712_d,Ne712d,Te712d]=dailydata(7,12,0,Select_data);
[D712_a,Ne712a,Te712a]=dailydata(7,12,1,Select_data);

[D713_d,Ne713d,Te713d]=dailydata(7,13,0,Select_data);
[D713_a,Ne713a,Te713a]=dailydata(7,13,1,Select_data);

 [D715_d,Ne715d,Te715d]=dailydata10(7,15,0,Select_data);
 [D715_a,Ne715a,Te715a]=dailydata10(7,15,1,Select_data);

 [D716_d,Ne716d,Te716d]=dailydata(7,16,0,Select_data);
 [D716_a,Ne716a,Te716a]=dailydata(7,16,1,Select_data);

 [D717_d,Ne717d,Te717d]=dailydata(7,17,0,Select_data);
 [D717_a,Ne717a,Te717a]=dailydata(7,17,1,Select_data);

 [D718_d,Ne718d,Te718d]=dailydata(7,13,0,Select_data);
 [D718_a,Ne718a,Te718a]=dailydata(7,13,1,Select_data);

 [D719_d,Ne719d,Te719d]=dailydata19(7,19,0,Select_data);
 [D719_a,Ne719a,Te719a]=dailydata19(7,19,1,Select_data);

%%
% orbit mean data
[Ne710_d,Te710_d]=orbitdata(D710_d,Select_data,0);
[Ne710_a,Te710_a]=orbitdata(D710_a,Select_data,1);
[Ne711_d,Te711_d]=orbitdata(D711_d,Select_data,0);
[Ne711_a,Te711_a]=orbitdata(D711_a,Select_data,0);
[Ne712_d,Te712_d]=orbitdata(D712_d,Select_data,0);
[Ne712_a,Te712_a]=orbitdata(D712_a,Select_data,1);
[Ne713_d,Te713_d]=orbitdata(D713_d,Select_data,0);
[Ne713_a,Te713_a]=orbitdata(D713_a,Select_data,1);

[Ne708_d,Te708_d]=orbitdata(D708_d,Select_data,0);
[Ne708_a,Te708_a]=orbitdata(D708_a,Select_data,1);

[Ne709_d,Te709_d]=orbitdata(D709_d,Select_data,0);
[Ne709_a,Te709_a]=orbitdata(D709_a,Select_data,1);

[Ne701_d,Te701_d]=orbitdata(D701_d,Select_data,0);
[Ne701_a,Te701_a]=orbitdata(D701_a,Select_data,1);

[Ne703_d,Te703_d]=orbitdata(D703_d,Select_data,0);
[Ne703_a,Te703_a]=orbitdata(D703_a,Select_data,1);

[Ne630_d,Te630_d]=orbitdata(D630_d,Select_data,0);
[Ne630_a,Te630_a]=orbitdata(D630_a,Select_data,1);

[Ne629_d,Te629_d]=orbitdata(D629_d,Select_data,0);
[Ne629_a,Te629_a]=orbitdata(D629_a,Select_data,1);

 [Ne715_d,Te715_d]=orbitdata(D715_d,Select_data,0);
 [Ne715_a,Te715_a]=orbitdata(D715_a,Select_data,1);

 [Ne716_d,Te716_d]=orbitdata(D716_d,Select_data,0);
 [Ne716_a,Te716_a]=orbitdata(D716_a,Select_data,1);

 [Ne717_d,Te717_d]=orbitdata(D717_d,Select_data,0);
 [Ne717_a,Te717_a]=orbitdata(D717_a,Select_data,1);

 [Ne718_d,Te718_d]=orbitdata(D718_d,Select_data,0);
[Ne718_a,Te718_a]=orbitdata(D718_a,Select_data,1);

 [Ne719_d,Te719_d]=orbitdata(D719_d,Select_data,0);
 [Ne719_a,Te719_a]=orbitdata(D719_a,Select_data,1);
 
 [D703_d,Ne703d,Te703d]=dailydata(7,3,0,Select_data);
[D703_a,Ne703a,Te703a]=dailydata(7,3,1,Select_data);

[D706_d,Ne706d,Te706d]=dailydata(7,6,0,Select_data);
[D706_a,Ne706a,Te706a]=dailydata(7,6,1,Select_data);

[D714_d,Ne714d,Te714d]=dailydata19(7,14,0,Select_data);
[D714_a,Ne714a,Te714a]=dailydata19(7,14,1,Select_data);

[Ne706_d,Te706_d]=orbitdata(D706_d,Select_data,0);
[Ne706_a,Te706_a]=orbitdata(D706_a,Select_data,1);


[D705_d,Ne705d,Te705d]=dailydata10(7,5,0,Select_data);
[D705_a,Ne705a,Te705a]=dailydata10(7,5,1,Select_data);

[Ne705_d,Te705_d]=orbitdata(D705_d,Select_data,0);
[Ne705_a,Te705_a]=orbitdata(D705_a,Select_data,1);

[D704_d,Ne704d,Te704d]=dailydata19(7,4,0,Select_data);
[D704_a,Ne704a,Te704a]=dailydata19(7,4,1,Select_data);

[Ne704_d,Te704_d]=orbitdata(D704_d,Select_data,0);
[Ne704_a,Te704_a]=orbitdata(D704_a,Select_data,1);

[D707_d,Ne707d,Te707d]=dailydata(7,7,0,Select_data);
[D707_a,Ne707a,Te707a]=dailydata(7,7,1,Select_data);

[Ne707_d,Te707_d]=orbitdata(D707_d,Select_data,0);
[Ne707_a,Te707_a]=orbitdata(D707_a,Select_data,1);

[D708_d,Ne708d,Te708d]=dailydata(7,8,0,Select_data);
[D708_a,Ne708a,Te708a]=dailydata(7,8,1,Select_data);

[Ne708_d,Te708_d]=orbitdata(D708_d,Select_data,0);
[Ne708_a,Te708_a]=orbitdata(D708_a,Select_data,1);

[D702_d,Ne702d,Te702d]=dailydata(7,2,0,Select_data);
[D702_a,Ne702a,Te702a]=dailydata(7,2,1,Select_data);

[Ne702_d,Te702_d]=orbitdata(D702_d,Select_data,0);
[Ne702_a,Te702_a]=orbitdata(D702_a,Select_data,1);

[Ne703_d,Te703_d]=orbitdata(D703_d,Select_data,0);
[Ne703_a,Te703_a]=orbitdata(D703_a,Select_data,1);

[Ne714_d,Te714_d]=orbitdata(D714_d,Select_data,0);
[Ne714_a,Te714_a]=orbitdata(D714_a,Select_data,1);

