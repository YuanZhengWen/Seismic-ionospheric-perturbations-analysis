%%  Read the data(.txt)
pap_filepath = 'F:\CSES\send\PAP\PAP';
pap_data_list = dir([pap_filepath,'\*.txt']); % READ the file(.txt)
pap_num_file = length(pap_data_list);
pap_data = [];

for k = 1: pap_num_file
    %location = h5read([filepath filename(k).name],'/LonLat');
    pap_data1 = importdata(pap_data_list(k).name);
    % Data selection
    Eq_lat = -0.52;
    Eq_lon = 128.17;
    
    pap_data = [pap_data;pap_data1];
end
%%
nStartLat = Eq_lat-20;
nEndLat = Eq_lat+20;
nStartLon = Eq_lon-20;
nEndLon = Eq_lon+20;


%% Selecting the data of JUNE 14th-30th
pap_Select_data1 = []; % JUNE 14th-30th
pap_orbit_selection=1; % 0 for descending; 1 for ascending

% Data Selection Process
for j = 1:length(pap_data)
    
    if ( (pap_data(j,11) >= nStartLat) && (pap_data(j,11) <= nEndLat) && (pap_data(j,10) >= nStartLon) && (pap_data(j,10) <= nEndLon)) % data(j,2)==0)
        
        if pap_data(j,4) == 6 && pap_data(j,5) >= 14
            pap_Select_data1 = [pap_Select_data1;pap_data(j,:)];
        end
        
    end
    
end
%% Selecting the data of JULY 1st- 14th
pap_Select_data2 = []; % JULY 1st- 14th
for j = 1:length(pap_data)
    
    if ((pap_data(j,11) >= nStartLat) && (pap_data(j,11) <= nEndLat) && (pap_data(j,10) >= nStartLon) && (pap_data(j,10) <= nEndLon))
        
        
        if pap_data(j,4) == 7 && pap_data(j,5) <= 17
            pap_Select_data2=[pap_Select_data2;pap_data(j,:)];
        end
    end
    
end
%% The data 30 days before the EQ
% Adding Select_data1 and Select_data2
pap_Select_data = [pap_Select_data1;pap_Select_data2];
%% PAP H and He density error data selection
for j=1:length(pap_Select_data)
    if pap_Select_data(j,12)<1
        pap_Select_data(j,12) = 0;
    end
    
    if pap_Select_data(j,13)<1
        pap_Select_data(j,13) = 0;
    end
end
%% Get the daily data of pap

% orbit observation data and grid data
[pap_D710_d,H710d,He710d,O710d,Ti710d]=papdailydata10(7,10,0,pap_Select_data);
[pap_D710_a,H710a,He710a,O710a,Ti710a]=papdailydata10(7,10,1,pap_Select_data);
%%
[pap_D701_d,H701d,He701d,O701d,Ti701d]=papdailydata(6,26,0,pap_Select_data);
[pap_D701_a,H701a,He701a,O701a,Ti701a]=papdailydata(6,26,1,pap_Select_data);

[pap_D630_d,H630d,He630d,O630d,Ti630d]=papdailydata10(6,30,0,pap_Select_data);
[pap_D630_a,H630a,He630a,O630a,Ti630a]=papdailydata10(6,30,1,pap_Select_data);

[pap_D629_d,H629d,He629d,O629d,Ti629d]=papdailydata19(6,29,0,pap_Select_data);
[pap_D629_a,H629a,He629a,O629a,Ti629a]=papdailydata19(6,29,1,pap_Select_data);

[pap_D703_d,H703d,He703d,O703d,Ti703d]=papdailydata(7,3,0,pap_Select_data);
[pap_D703_a,H703a,He703a,O703a,Ti703a]=papdailydata(7,3,1,pap_Select_data);

[pap_D711_d,H711d,He711d,O711d,Ti711d]=papdailydata(7,11,0,pap_Select_data);
[pap_D711_a,H711a,He711a,O711a,Ti711a]=papdailydata(7,11,1,pap_Select_data);

[pap_D712_d,H712d,He712d,O712d,Ti712d]=papdailydata(7,12,0,pap_Select_data);
[pap_D712_a,H712a,He712a,O712a,Ti712a]=papdailydata(7,12,1,pap_Select_data);

[pap_D713_d,H713d,He713d,O713d,Ti713d]=papdailydata(7,13,0,pap_Select_data);
[pap_D713_a,H713a,He713a,O713a,Ti713a]=papdailydata(7,13,1,pap_Select_data);

[pap_D706_d,H706d,He706d,O706d,Ti706d]=papdailydata10(7,6,0,pap_Select_data);
[pap_D706_a,H706a,He706a,O706a,Ti706a]=papdailydata10(7,6,1,pap_Select_data);

[pap_D705_d,H705d,He705d,O705d,Ti705d]=papdailydata10(7,5,0,pap_Select_data);
[pap_D705_a,H705a,He705a,O705a,Ti705a]=papdailydata10(7,5,1,pap_Select_data);

[pap_D704_d,H704d,He704d,O704d,Ti704d]=papdailydata4(7,4,0,pap_Select_data);
[pap_D704_a,H704a,He704a,O704a,Ti704a]=papdailydata10(7,4,1,pap_Select_data);

[pap_D702_d,H702d,He702d,O702d,Ti702d]=papdailydata10(7,2,0,pap_Select_data);
[pap_D702_a,H702a,He702a,O702a,Ti702a]=papdailydata10(7,2,1,pap_Select_data);

[pap_D707_d,H707d,He707d,O707d,Ti707d]=papdailydata10(7,7,0,pap_Select_data);
[pap_D707_a,H707a,He707a,O707a,Ti707a]=papdailydata10(7,7,1,pap_Select_data);

[pap_D708_d,H708d,He708d,O708d,Ti708d]=papdailydata(7,8,0,pap_Select_data);
[pap_D708_a,H708a,He708a,O708a,Ti708a]=papdailydata(7,8,1,pap_Select_data);

[pap_D709_d,H709d,He709d,O707d,Ti709d]=papdailydata19(7,9,0,pap_Select_data);
[pap_D709_a,H709a,He709a,O709a,Ti709a]=papdailydata19(7,9,1,pap_Select_data);

[pap_D715_d,H715d,He715d,O715d,Ti715d]=papdailydata10(7,15,0,pap_Select_data);
[pap_D715_a,H715a,He715a,O715a,Ti715a]=papdailydata10(7,15,1,pap_Select_data);

[pap_D716_d,H716d,He716d,O716d,Ti716d]=papdailydata(7,16,0,pap_Select_data);
[pap_D716_a,H716a,He716a,O716a,Ti716a]=papdailydata(7,16,1,pap_Select_data);

[pap_D717_d,H717d,He717d,O717d,Ti717d]=papdailydata(7,17,0,pap_Select_data);
[pap_D717_a,H717a,He717a,O717a,Ti717a]=papdailydata(7,17,1,pap_Select_data);

[pap_D718_d,H718d,He718d,O718d,Ti718d]=papdailydata(7,13,0,pap_Select_data);
[pap_D718_a,H718a,He718a,O718a,Ti718a]=papdailydata(7,13,1,pap_Select_data);

[pap_D714_d,H714d,He714d,O714d,Ti714d]=papdailydata19(7,14,0,pap_Select_data);
[pap_D714_a,H714a,He714a,O714a,Ti714a]=papdailydata19(7,14,1,pap_Select_data);
%%
[H_0701D,He_0701D]=dataselect(pap_D701_d);
[H_0701A,He_0701A]=dataselect(pap_D701_a);

[H_0701D,He_0701D]=dataselect(pap_D701_d);
[H_0701A,He_0701A]=dataselect(pap_D701_a);

[H_0703D,He_0703D]=dataselect(pap_D703_d);
[H_0703A,He_0703A]=dataselect(pap_D703_a);

[H_0630D,He_0630D]=dataselect(pap_D630_d);
[H_0630A,He_0630A]=dataselect(pap_D630_a);

[H_0629D,He_0629D]=dataselect(pap_D629_d);
[H_0629A,He_0629A]=dataselect(pap_D629_a);

[H_0710D,He_0710D]=dataselect(pap_D710_d);
[H_0710A,He_0710A]=dataselect(pap_D710_a);

[H_0710D,He_0710D]=dataselect(pap_D710_d);
[H_0710A,He_0710A]=dataselect(pap_D710_a);
[H_0711D,He_0711D]=dataselect(pap_D711_d);
[H_0711A,He_0711A]=dataselect(pap_D711_a);
[H_0712D,He_0712D]=dataselect(pap_D712_d);
[H_0712A,He_0712A]=dataselect(pap_D712_a);
[H_0713D,He_0713D]=dataselect(pap_D713_d);
[H_0713A,He_0713A]=dataselect(pap_D713_a);

[H_0715D,He_0715D]=dataselect(pap_D715_d);
[H_0715A,He_0715A]=dataselect(pap_D715_a);
[H_0716D,He_0716D]=dataselect(pap_D716_d);
[H_0716A,He_0716A]=dataselect(pap_D716_a);
[H_0717D,He_0717D]=dataselect(pap_D717_d);
[H_0717A,He_0717A]=dataselect(pap_D717_a);
[H_0718D,He_0718D]=dataselect(pap_D718_d);
[H_0718A,He_0718A]=dataselect(pap_D718_a);
[H_0714D,He_0714D]=dataselect(pap_D714_d);
[H_0714A,He_0714A]=dataselect(pap_D714_a);


[H_0709D,He_0709D]=dataselect(pap_D709_d);
[H_0709A,He_0709A]=dataselect(pap_D709_a);
[H_0708D,He_0708D]=dataselect(pap_D708_d);
[H_0708A,He_0708A]=dataselect(pap_D708_a);

[H_0706D,He_0706D]=dataselect(pap_D706_d);
[H_0706A,He_0706A]=dataselect(pap_D706_a);
[H_0705D,He_0705D]=dataselect(pap_D705_d);
[H_0705A,He_0705A]=dataselect(pap_D705_a);
[H_0704D,He_0704D]=dataselect(pap_D704_d);
[H_0704A,He_0704A]=dataselect(pap_D704_a);
[H_0707D,He_0707D]=dataselect(pap_D707_d);
[H_0707A,He_0707A]=dataselect(pap_D707_a);
[H_0702D,He_0702D]=dataselect(pap_D702_d);
[H_0702A,He_0702A]=dataselect(pap_D702_a);


%% He Density Orbit Mean Data
[He701_d]=paporbitdataHe(He_0701D,pap_Select_data,0);
[He701_a]=paporbitdataHe(He_0701A,pap_Select_data,1);

[He703_d]=paporbitdataHe(He_0703D,pap_Select_data,0);
[He703_a]=paporbitdataHe(He_0703A,pap_Select_data,1);

[He630_d]=paporbitdataHe(He_0630D,pap_Select_data,0);
[He630_a]=paporbitdataHe(He_0630A,pap_Select_data,1);

[He629_d]=paporbitdataHe(He_0629D,pap_Select_data,0);
[He629_a]=paporbitdataHe(He_0629A,pap_Select_data,1);


[He710_d]=paporbitdataHe(He_0710D,pap_Select_data,0);
[He710_a]=paporbitdataHe(He_0710A,pap_Select_data,1);
[He711_d]=paporbitdataHe(He_0711D,pap_Select_data,0);
[He711_a]=paporbitdataHe(He_0711A,pap_Select_data,1);
[He712_d]=paporbitdataHe(He_0712D,pap_Select_data,0);
[He712_a]=paporbitdataHe(He_0712A,pap_Select_data,1);
[He713_d]=paporbitdataHe(He_0713D,pap_Select_data,0);
[He713_a]=paporbitdataHe(He_0713A,pap_Select_data,1);

[He715_d]=paporbitdataHe(He_0715D,pap_Select_data,0);
[He715_a]=paporbitdataHe(He_0715A,pap_Select_data,1);
[He716_d]=paporbitdataHe(He_0716D,pap_Select_data,0);
[He716_a]=paporbitdataHe(He_0716A,pap_Select_data,1);
[He717_d]=paporbitdataHe(He_0717D,pap_Select_data,0);
[He717_a]=paporbitdataHe(He_0717A,pap_Select_data,1);
[He718_d]=paporbitdataHe(He_0718D,pap_Select_data,0);
[He718_a]=paporbitdataHe(He_0718A,pap_Select_data,1);
[He714_d]=paporbitdataHe(He_0714D,pap_Select_data,0);
[He714_a]=paporbitdataHe(He_0714A,pap_Select_data,1);

[He706_d]=paporbitdataHe(He_0706D,pap_Select_data,0);
[He706_a]=paporbitdataHe(He_0706A,pap_Select_data,1);
[He705_d]=paporbitdataHe(He_0705D,pap_Select_data,0);
[He705_a]=paporbitdataHe(He_0705A,pap_Select_data,1);
[He704_d]=paporbitdataHe(He_0704D,pap_Select_data,0);
[He704_a]=paporbitdataHe(He_0704A,pap_Select_data,1);
[He707_d]=paporbitdataHe(He_0707D,pap_Select_data,0);
[He707_a]=paporbitdataHe(He_0707A,pap_Select_data,1);
[He702_d]=paporbitdataHe(He_0702D,pap_Select_data,0);
[He702_a]=paporbitdataHe(He_0702A,pap_Select_data,1);
[He708_d]=paporbitdataHe(He_0708D,pap_Select_data,0);
[He708_a]=paporbitdataHe(He_0708A,pap_Select_data,1);
[He709_d]=paporbitdataHe(He_0709D,pap_Select_data,0);
[He709_a]=paporbitdataHe(He_0709A,pap_Select_data,1);
%% H Density Orbit Mean Data
[H710_d]=paporbitdataH(H_0710D,pap_Select_data);
[H710_a]=paporbitdataH(H_0710A,pap_Select_data);
[H711_d]=paporbitdataH(H_0711D,pap_Select_data);
[H711_a]=paporbitdataH(H_0711A,pap_Select_data);
[H712_d]=paporbitdataH(H_0712D,pap_Select_data);
[H712_a]=paporbitdataH(H_0712A,pap_Select_data);
[H713_d]=paporbitdataH(H_0713D,pap_Select_data);
[H713_a]=paporbitdataH(H_0713A,pap_Select_data);
%%
% orbit mean data
[~,~,O701_d,Ti701_d]=paporbitdata(pap_D701_d,pap_Select_data,0);
[~,~,O701_a,Ti701_a]=paporbitdata(pap_D701_a,pap_Select_data,1);

[~,~,O703_d,Ti703_d]=paporbitdata(pap_D703_d,pap_Select_data,0);
[~,~,O703_a,Ti703_a]=paporbitdata(pap_D703_a,pap_Select_data,1);

[~,~,O630_d,Ti630_d]=paporbitdata(pap_D630_d,pap_Select_data,0);
[~,~,O630_a,Ti630_a]=paporbitdata(pap_D630_a,pap_Select_data,1);

[~,~,O629_d,Ti629_d]=paporbitdata(pap_D629_d,pap_Select_data,0);
[~,~,O629_a,Ti629_a]=paporbitdata(pap_D629_a,pap_Select_data,1);

[~,~,O715_d,Ti715_d]=paporbitdata(pap_D715_d,pap_Select_data,0);
[~,~,O715_a,Ti715_a]=paporbitdata(pap_D715_a,pap_Select_data,1);
[~,~,O716_d,Ti716_d]=paporbitdata(pap_D716_d,pap_Select_data,0);
[~,~,O716_a,Ti716_a]=paporbitdata(pap_D716_a,pap_Select_data,1);
[~,~,O717_d,Ti717_d]=paporbitdata(pap_D717_d,pap_Select_data,0);
[~,~,O717_a,Ti717_a]=paporbitdata(pap_D717_a,pap_Select_data,1);
[~,~,O718_d,Ti718_d]=paporbitdata(pap_D713_d,pap_Select_data,0);
[~,~,O718_a,Ti718_a]=paporbitdata(pap_D713_a,pap_Select_data,1);
[~,~,O714_d,Ti714_d]=paporbitdata(pap_D714_d,pap_Select_data,0);
[~,~,O714_a,Ti714_a]=paporbitdata(pap_D714_a,pap_Select_data,1);



[~,~,O710_d,Ti710_d]=paporbitdata(pap_D710_d,pap_Select_data,0);
[~,~,O710_a,Ti710_a]=paporbitdata(pap_D710_a,pap_Select_data,1);
[~,~,O711_d,Ti711_d]=paporbitdata(pap_D711_d,pap_Select_data,0);
[~,~,O711_a,Ti711_a]=paporbitdata(pap_D711_a,pap_Select_data,1);
[~,~,O712_d,Ti712_d]=paporbitdata(pap_D712_d,pap_Select_data,0);
[~,~,O712_a,Ti712_a]=paporbitdata(pap_D712_a,pap_Select_data,1);
[~,~,O713_d,Ti713_d]=paporbitdata(pap_D713_d,pap_Select_data,0);
[~,~,O713_a,Ti713_a]=paporbitdata(pap_D713_a,pap_Select_data,1);

[~,~,O706_d,Ti706_d]=paporbitdata(pap_D706_d,pap_Select_data,0);
[~,~,O706_a,Ti706_a]=paporbitdata(pap_D706_a,pap_Select_data,1);
[~,~,O705_d,Ti705_d]=paporbitdata(pap_D705_d,pap_Select_data,0);
[~,~,O705_a,Ti705_a]=paporbitdata(pap_D705_a,pap_Select_data,1);
[~,~,O704_d,Ti704_d]=paporbitdata(pap_D704_d,pap_Select_data,0);
[~,~,O704_a,Ti704_a]=paporbitdata(pap_D704_a,pap_Select_data,1);

[~,~,O702_d,Ti702_d]=paporbitdata(pap_D702_d,pap_Select_data,0);
[~,~,O702_a,Ti702_a]=paporbitdata(pap_D702_a,pap_Select_data,1);

[~,~,O707_d,Ti707_d]=paporbitdata(pap_D707_d,pap_Select_data,0);
[~,~,O707_a,Ti707_a]=paporbitdata(pap_D707_a,pap_Select_data,1);

[~,~,O708_d,Ti708_d]=paporbitdata(pap_D708_d,pap_Select_data,0);
[~,~,O708_a,Ti708_a]=paporbitdata(pap_D708_a,pap_Select_data,1);

[~,~,O709_d,Ti709_d]=paporbitdata(pap_D709_d,pap_Select_data,0);
[~,~,O709_a,Ti709_a]=paporbitdata(pap_D709_a,pap_Select_data,1);
