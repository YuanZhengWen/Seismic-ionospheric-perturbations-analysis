function [pap_anv_H,pap_anv_He,pap_anv_O,pap_anv_Ti]=paporbitdata(Daydata,pap_Select_data,orbit)

% orbit_start_lat=D713_d(:,7)-1;
% orbit_end_lat=D713_d(:,7)+1;
% orbit_start_lon=D713_d(:,6)-2 ;
% orbit_end_lon=D713_d(:,6)+2;
[a,~] = size(Daydata);
pap_anv_O=[];
pap_anv_He=[];
pap_anv_H=[];
pap_anv_Ti=[];
pap_sum_O = 0;
pap_sum_He = 0;
pap_sum_H=0;
pap_sum_Ti=0;
count = 0;
ct = 0;
ct2 = 0;
ct3 = 0;   
for i=1:a 
       
    for j=1:length(pap_Select_data)          
        
        if (pap_Select_data(j,10)>= Daydata(i,10)-2) && (pap_Select_data(j,10)<= Daydata(i,10)+2) && (pap_Select_data(j,11)>=Daydata(i,11)-1) && (pap_Select_data(j,11)<=Daydata(i,11)+1) && (pap_Select_data(j,2)==orbit)
            ct = ct + 1;
            pap_sum_O=pap_sum_O+pap_Select_data(j,14);
            pap_sum_Ti=pap_sum_Ti+pap_Select_data(j,15);
%             pap_sum_H =  pap_sum_H  + pap_Select_data(j,12);
%             pap_sum_He =  pap_sum_He  + pap_Select_data(j,13);
             if pap_Select_data(j,12) ~= 0
                 ct2 = ct2 +1;
                pap_sum_H =  pap_sum_H  + pap_Select_data(j,12);
             end
            if pap_Select_data(j,13) ~= 0
                 ct3 = ct3 + 1;
             pap_sum_He =  pap_sum_He  + pap_Select_data(j,13);
            end
        end
        pap_anv_H(i) = pap_sum_H/ ct2;
        pap_anv_He(i) = pap_sum_He/ ct3;   
        pap_anv_O(i)=pap_sum_O/ct;
        pap_anv_Ti(i)=pap_sum_Ti/ct;
    end
end

end