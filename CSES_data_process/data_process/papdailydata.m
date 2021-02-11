% Function for obtaining obeservation data
% Daily average data of electron density (Ne) and temperature (Te)
% Input: month, day, orbit (0 (D) or 1 (A)) 
% Output: Daily data; 
% [observation_data,value_Ne_day,value_Te_day]=Dayily_data(7,13,1,Select_data)
function [pap_observation_data,value_H_day,value_He_day,value_O_day,value_Ti_day]=papdailydata(month,day,orbit,pap_Select_data)
pap_observation_data = [];
 [q,~] = size(pap_Select_data);
 
 for g = 1:q
     
        if (pap_Select_data(g,4) == month) && (pap_Select_data(g,5) == day && pap_Select_data(g,2) == orbit && pap_Select_data(g,10)>=118 && pap_Select_data(g,10)<=141 ) 
            pap_observation_data = [pap_observation_data;pap_Select_data(g,:)];
        end
        
 end
Eq_lat = -0.52;
Eq_lon = 128.17;
nStartLat = Eq_lat-20;
nEndLat = Eq_lat+20;
nStartLon = Eq_lon-20;
nEndLon = Eq_lon+20;

% nStartLat = Eq_lat-15;
% nEndLat = Eq_lat+15;
% nStartLon = Eq_lon-15;
% nEndLon = Eq_lon+15;

nSplitLat = 1;
nSplitLon = 2;
% nSplitLat = 2;
% nSplitLon = 4;

nLon1 = nStartLon :nSplitLon:nEndLon;

nLat1 = 20 :-nSplitLat:-20;
% important
[X1,Y1] = meshgrid(nLon1,nLat1');
[nY1Len,nX1Len] = size(X1);
Z1 = NaN(nY1Len,nX1Len);  

% Choosing the centain day:0713
% The same thounghts as Calc each value in every cell

[c,~] = size(pap_observation_data);
[e,f] = size(X1);

sum_He_day = zeros(1,e*f);
sum_H_day = zeros(1,e*f);
sum_O_day = zeros(1,e*f);
sum_Ti_day = zeros(1,e*f);
pap_anv_H = zeros(1,e*f);
pap_anv_He = zeros(1,e*f);
pap_anv_O= zeros(1,e*f);
pap_anv_Ti= zeros(1,e*f);
count_day = zeros(1,e*f);
count1 = 0;

for lon = nStartLon :nSplitLon:nEndLon   
    for lat = nEndLat :-nSplitLat: nStartLat
        count1 = count1 + 1;
        ct = 0;
        for i = 1: c          
            if  ( pap_observation_data(i,10) >= lon )  && (pap_observation_data(i,10) < lon + nSplitLon ) && (pap_observation_data(i,11) >= lat )  && (pap_observation_data(i,11) < lat + nSplitLat  )
                ct = ct + 1;
                sum_H_day(1,count1) = sum_H_day(1,count1) + pap_observation_data(i,12);
                sum_He_day(1,count1) = sum_He_day(1,count1) + pap_observation_data(i,13);
                sum_O_day(1,count1) = sum_O_day(1,count1) + pap_observation_data(i,14);
                sum_Ti_day(1,count1) = sum_Ti_day(1,count1) + pap_observation_data(i,15);
            end
            count(1,count1) = ct;
            pap_anv_H(1,count1) = sum_H_day(1,count1) / ct;
            pap_anv_He(1,count1) = sum_He_day(1,count1) / ct;
            pap_anv_O(1,count1) = sum_O_day(1,count1) / ct;
            pap_anv_Ti(1,count1) = sum_Ti_day(1,count1) / ct;
        end
    end  
end
% Changing the one-dimensional anv_Ne(Te)_day into two-dimensional value_Ne(Te)_day

value_H_day = zeros(e,f);
value_He_day = zeros(e,f);
value_O_day = zeros(e,f);
value_Ti_day = zeros(e,f);
[e,f] = size(X1);

for i = 1:e
    
    for j = 1:f
        value_H_day(i,j) = pap_anv_H((j-1)*e + i);
        value_He_day(i,j) = pap_anv_He((j-1)*e + i);
        value_O_day(i,j) = pap_anv_O((j-1)*e + i);
        value_Ti_day(i,j) = pap_anv_Ti((j-1)*e + i);
    
end
end