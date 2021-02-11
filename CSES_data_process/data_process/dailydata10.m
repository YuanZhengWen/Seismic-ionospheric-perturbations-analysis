% Function for obtaining obeservation data
% Daily average data of electron density (Ne) and temperature (Te)
% Input: month, day, orbit (0 (D) or 1 (A)) 
% Output: Daily data; 
% [observation_data,value_Ne_day,value_Te_day]=Dayily_data(7,13,1,Select_data)
function [observation_data,value_Ne_day,value_Te_day]=dailydata10(month,day,orbit,Select_data)
observation_data = [];
 [q,~] = size(Select_data);
 
 for g = 1:q
     
        if (Select_data(g,4) == month) && (Select_data(g,5) == day && Select_data(g,2) == orbit && Select_data(g,6)>=110 && Select_data(g,6)<=138 ) 
            observation_data = [observation_data;Select_data(g,:)];
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

[c,~] = size(observation_data);
[e,f] = size(X1);

sum_Ne_day = zeros(1,e*f);
sum_Te_day = zeros(1,e*f);
anv_Te_day = zeros(1,e*f);
anv_Ne_day = zeros(1,e*f);
count_day = zeros(1,e*f);
count1 = 0;

for lon = nStartLon :nSplitLon:nEndLon   
    for lat = nEndLat :-nSplitLat: nStartLat
        count1 = count1 + 1;
        ct = 0;
        for i = 1: c          
            if  ( observation_data(i,6) >= lon )  && (observation_data(i,6) < lon + nSplitLon ) && (observation_data(i,7) >= lat )  && (observation_data(i,7) < lat + nSplitLat  )
                ct = ct + 1;
                sum_Ne_day(1,count1) = sum_Ne_day(1,count1) + observation_data(i,8);
                sum_Te_day(1,count1) = sum_Te_day(1,count1) + observation_data(i,9);
            end
            count(1,count1) = ct;
            anv_Te_day(1,count1) = sum_Te_day(1,count1) / ct;
            anv_Ne_day(1,count1) = sum_Ne_day(1,count1) / ct;
        end
    end  
end
% Changing the one-dimensional anv_Ne(Te)_day into two-dimensional value_Ne(Te)_day

value_Ne_day = zeros(e,f);
value_Te_day = zeros(e,f);
[e,f] = size(X1);

for i = 1:e
    
    for j = 1:f
        value_Ne_day(i,j) = anv_Ne_day((j-1)*e + i);
        value_Te_day(i,j) = anv_Te_day((j-1)*e + i);
    end
    
end
end