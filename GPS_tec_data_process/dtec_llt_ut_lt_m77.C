//===============================================================
//Demeter Case Code "dtec_llt_ut_lt_m77.C" 
//TAO DAN 2017.01.04 updated
//Wang 2020.08.15 updated

//===============================================================
#include <math.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <TText.h>
#include <TFile.h>
#include <TChain.h>
#include <TH1.h>
#include <TH2.h>
#include <TMultiGraph.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TLatex.h>
#include <TMarker.h>
#include <TGaxis.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TStyle.h>
#include <fstream> 
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <iomanip>
#include <string.h>
#include <TColor.h>
#include <TPaletteAxis.h>
#include <TLegend.h>
#include <TBox.h>
#include <TArrow.h>
#include <TRandom.h>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>
#include <THistPainter.h>
#include <TEllipse.h>

using namespace std;

//g++ dtec_llt_ut_lt_m77.C -O2 -Wall -fPIC -pthread -m64 -I/home/wang/CERNroot/root/root/include/ -L/home/wang/CERNroot/root/root/lib/ -lGraf -lGraf3d -lCore -lCint -lTree -lRIO  -lHist -lGpad -o dtec_llt_ut_lt_m77

#define PI 3.141592653589793238463 

float vtec_30d_quartile_h02_on[5183][4];
float vtec_30d_quartile_h24_on[5183][4];
float vtec_30d_quartile_h46_on[5183][4];

int main(int argc, char** argv) 
{
float map_calc_quartile(int ngridx, int vind, int nq, float *qteccode);
float map_bilinear_interpolation(float *datapoint4); 

TH2F *hm = new TH2F("hm","Mercator",360,-180,180,175,-87.5,87.5);

ifstream in1;
in1.open("/home/wang/Projects/EQ/otherdata/earth.dat");
double maplon,maplat;

 while(in1.good())
{
 in1>>maplon>>maplat;

 hm->Fill(maplon,maplat,1);
} 
 in1.close();

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//read GPS-TEC-2HR data 30d-before and 10d-after
float teclat[71],teclon[73];
for(int i=0;i<71;i++){teclat[i]=87.5-2.5*i;}
for(int i=0;i<73;i++){teclon[i]=-180+5.0*i;}

//GPS-TEC: NWC ON; EQ7.7-(107.419,-9.284)-20.0km-south of Java-06/07/17-08:19:26UTC

int arr_h02_on[62],arr_h24_on[62],arr_h46_on[62];
for(int i=0;i<31;i++){arr_h02_on[2*i]=468+12*i+3;arr_h02_on[2*i+1]=468+12*i+4;}
//for(int i=0;i<31;i++){arr_h02_on[2*i]=528+12*i+4;arr_h02_on[2*i+1]=528+12*i+5;}

for(int i=0;i<31;i++){arr_h24_on[2*i]=528+12*i+4;arr_h24_on[2*i+1]=528+12*i+5;}
//for(int i=0;i<31;i++){arr_h24_on[2*i]=552+12*i+0;arr_h24_on[2*i+1]=552+12*i+1;}
//for(int i=0;i<31;i++){arr_h24_on[2*i]=552+12*i+3;arr_h24_on[2*i+1]=552+12*i+4;}

//for(int i=0;i<31;i++){arr_h46_on[2*i]=552+12*i+0;arr_h46_on[2*i+1]=552+12*i+1;}
for(int i=0;i<31;i++){arr_h46_on[2*i]=552+12*i+3;arr_h46_on[2*i+1]=552+12*i+4;}
//for(int i=0;i<31;i++){arr_h46_on[2*i]=636+12*i+8;arr_h46_on[2*i+1]=636+12*i+9;}

float vtec_max_h02_on[5183],vtec_max_h24_on[5183],vtec_max_h46_on[5183];
float vdtec_h02_on[5183],vdtec_h24_on[5183],vdtec_h46_on[5183];

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//loop for three TEC anomaly grid
for(int ngrid=0;ngrid<3;ngrid++){
ifstream in2;
in2.open("/home/wang/Projects/EQ/Yahama/GPS_TEC2HR_IGS_Yahama_20190430_0724.txt");
//EQ:06-07-17 08:19:26UTC; +75:05-03; +60:05-18; +45:06-02; +30:06-17; -10:07-27  
//dtec extreme maximum time ponit: 03:00,05:00,13:00 07-15; 30d before 07-17: 06-17
//06/17(540)-07/17(912)--31d
int nrecord2=62;   //01:00([0-2]UTC) in per 31 days
float tecyear_on[nrecord2],tecsec_on[nrecord2],teccode_on[nrecord2][5183];

for(int i=0;i<nrecord2;i++){
    tecyear_on[i]=0.;tecsec_on[i]=0.;
    for(int j=0;j<5183;j++){teccode_on[i][j]=0.;}
}

//char arr2[90000];
//int m2=1;
int nrd2=0;
int nrd2_temp0=0;
float vtecyeartemp_on=0., vtecsectemp_on=0., vtectemp_on=0.;

 /*while(m2<119)
{
  in2.getline(arr2,90000);
  m2++;
}
 in2.getline(arr2,90000); */

while(in2.good())
{
 int nrd2_temp1=0;

   if(ngrid==0){
     for(int k1=0;k1<62;k1++){
        if(nrd2==arr_h02_on[k1]){nrd2_temp1=nrd2;break;}
        else{nrd2_temp1=0;}
     }
   }
   if(ngrid==1){
     for(int k1=0;k1<62;k1++){
        if(nrd2==arr_h24_on[k1]){nrd2_temp1=nrd2;break;}
        else{nrd2_temp1=0;}
     }
   }
   if(ngrid==2){
     for(int k1=0;k1<62;k1++){
        if(nrd2==arr_h46_on[k1]){nrd2_temp1=nrd2;break;}
        else{nrd2_temp1=0;}
     }
   }

 if(nrd2_temp1!=0){
    //cout<<nrd2_temp1<<endl;
    in2>>tecyear_on[nrd2_temp0]>>tecsec_on[nrd2_temp0];
    for(int k2=0;k2<5183;k2++){in2>>teccode_on[nrd2_temp0][k2];}
    nrd2_temp0++;
 }
 else{
    in2>>vtecyeartemp_on>>vtecsectemp_on;
    for(int k3=0;k3<5183;k3++){in2>>vtectemp_on;}
 }

 nrd2++;
} 
in2.close();
//cout<<tecyear_on[61]<<";"<<tecsec_on[61]<<";"<<teccode_on[61][0]<<";"<<teccode_on[61][5182]<<";"<<endl;

//TEC ; the global map: grid 5040=72*70
float teccode_c2h_on[31][5183];  //372=31*12
for(int i=0;i<31;i++){
   for(int j=0;j<5183;j++){
      teccode_c2h_on[i][j]=0.;
      teccode_c2h_on[i][j]=0.5*(teccode_on[2*i][j]+teccode_on[2*i+1][j]); //mid-piont value of 2h time period
   }
}
//cout<<teccode_map_grid_on[30][71][69]<<endl;

float vtec_max_on[5183]; //vtec_max_on[*][*]=teccode_map_grid_on[30][*][*]
 
for(int l1=0;l1<5183;l1++){
                float tec_quartile_temp_on[30];
                	  if(ngrid==0){ vtec_max_on[l1]=teccode_c2h_on[28][l1];}
		          if(ngrid==1){ vtec_max_on[l1]=teccode_c2h_on[28][l1];}
		          if(ngrid==2){ vtec_max_on[l1]=teccode_c2h_on[28][l1];}

                for(int l2=0;l2<30;l2++){
                tec_quartile_temp_on[l2]=0.;
                tec_quartile_temp_on[l2]=teccode_c2h_on[l2][l1];
                }

   map_calc_quartile(ngrid,l1,30,tec_quartile_temp_on);
   //cout<<vtec_max_on[l1]<<";"<<vtec_30d_quartile_on[l1][1]<<";"<<vtec_30d_quartile_on[l1][3]<<";"<<endl;
   
   //cout<<vtec_max_on[0][0]<<endl;
 }

if(ngrid==0){
  for(int l2=0;l2<5183;l2++){
      vtec_max_h02_on[l2]=vtec_max_on[l2];
      if(vtec_30d_quartile_h02_on[l2][3]!=0.){vdtec_h02_on[l2]=(vtec_max_h02_on[l2]-vtec_30d_quartile_h02_on[l2][1])/vtec_30d_quartile_h02_on[l2][3];}
  }
}
if(ngrid==1){
  for(int l2=0;l2<5183;l2++){
      vtec_max_h24_on[l2]=vtec_max_on[l2];
      if(vtec_30d_quartile_h24_on[l2][3]!=0.){vdtec_h24_on[l2]=(vtec_max_h24_on[l2]-vtec_30d_quartile_h24_on[l2][1])/vtec_30d_quartile_h24_on[l2][3];}
  }
}
if(ngrid==2){
  for(int l2=0;l2<5183;l2++){
      vtec_max_h46_on[l2]=vtec_max_on[l2];
      if(vtec_30d_quartile_h46_on[l2][3]!=0.){vdtec_h46_on[l2]=(vtec_max_h46_on[l2]-vtec_30d_quartile_h46_on[l2][1])/vtec_30d_quartile_h46_on[l2][3];}
  }
}

}
//end loop+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

for(int l3=0;l3<5183;l3++){
   if(vtec_max_h02_on[l3]==0. || vtec_max_h24_on[l3]==0. || vtec_max_h46_on[l3]==0. ){
      //vtec_max_h02_on[l3]=1.0E-5;vtec_max_h24_on[l3]=1.0E-5;vtec_max_h46_on[l3]=1.0E-5;
   }
   for(int l4=0;l4<4;l4++){
     if(vtec_30d_quartile_h02_on[l3][l4]==0. || vtec_30d_quartile_h24_on[l3][l4]==0. || vtec_30d_quartile_h46_on[l3][l4]==0. ){ 
        //vtec_30d_quartile_h02_on[l3][l4]=1.0E-5;vtec_30d_quartile_h24_on[l3][l4]=1.0E-5;vtec_30d_quartile_h46_on[l3][l4]=1.0E-5;
     }
   }  
   if(vdtec_h02_on[l3]>-2. && vdtec_h02_on[l3]<2.){vdtec_h02_on[l3]=-5.5;} //for changing the background color of the plot
   if(vdtec_h24_on[l3]>-2. && vdtec_h24_on[l3]<2.){vdtec_h24_on[l3]=-5.5;}
   if(vdtec_h46_on[l3]>-2. && vdtec_h46_on[l3]<2.){vdtec_h46_on[l3]=-5.5;} 
}

TH2D *h00 = new TH2D("h00"," ",73,-182.5,182.5,71,-88.75,88.75);
TH2D *h01 = new TH2D("h01"," ",73,-182.5,182.5,71,-88.75,88.75);
TH2D *h02 = new TH2D("h02"," ",73,-182.5,182.5,71,-88.75,88.75);
TH2D *h03 = new TH2D("h03"," ",73,-182.5,182.5,71,-88.75,88.75);
TH2D *h10 = new TH2D("h10"," ",73,-182.5,182.5,71,-88.75,88.75);
TH2D *h11 = new TH2D("h11"," ",73,-182.5,182.5,71,-88.75,88.75);
TH2D *h12 = new TH2D("h12"," ",73,-182.5,182.5,71,-88.75,88.75);
TH2D *h13 = new TH2D("h13"," ",73,-182.5,182.5,71,-88.75,88.75);
TH2D *h20 = new TH2D("h20"," ",73,-182.5,182.5,71,-88.75,88.75);
TH2D *h21 = new TH2D("h21"," ",73,-182.5,182.5,71,-88.75,88.75);
TH2D *h22 = new TH2D("h22"," ",73,-182.5,182.5,71,-88.75,88.75);
TH2D *h23 = new TH2D("h23"," ",73,-182.5,182.5,71,-88.75,88.75);

float loneq_on[1],lateq_on[1],loneq_size_on[1],lateq_size_on[1];
loneq_on[0]=128.17;lateq_on[0]=-0.52;
loneq_size_on[0]=138.35;lateq_size_on[0]=-0.52;
float lon_zoom_on[5]={95.,170.,170.,95.,95.},lat_zoom_on[5]={-22.,-22.,18.,18.,-22.};

//draw the magnetic equator
float mag_eq_lon[73];
float mag_eq_lat[73];


ifstream in3;
in3.open("/home/wang/Projects/EQ/otherdata/Mag_Equator.txt");
int n_mag_eq=0;

 while(in3.good())
{
 in3>>mag_eq_lat[n_mag_eq]>>mag_eq_lon[n_mag_eq];
 n_mag_eq++;
} 
 in3.close();

TGraph *g_mag_eq = new TGraph(73,mag_eq_lon,mag_eq_lat);

TGraph *g00 = new TGraph(1,loneq_on,lateq_on);
TGraph *g01 = new TGraph(1,loneq_size_on,lateq_size_on);
TGraph *g02 = new TGraph(5,lon_zoom_on,lat_zoom_on);
 
for(int k2=0;k2<71;k2++){  
for(int k1=0;k1<73;k1++){
    h00->Fill(teclon[k1],teclat[k2],vtec_max_h02_on[k1+73*k2]);
    h10->Fill(teclon[k1],teclat[k2],vtec_max_h24_on[k1+73*k2]);
    h20->Fill(teclon[k1],teclat[k2],vtec_max_h46_on[k1+73*k2]);

    h01->Fill(teclon[k1],teclat[k2],vtec_30d_quartile_h02_on[k1+73*k2][1]);
    h11->Fill(teclon[k1],teclat[k2],vtec_30d_quartile_h24_on[k1+73*k2][1]);
    h21->Fill(teclon[k1],teclat[k2],vtec_30d_quartile_h46_on[k1+73*k2][1]);

       h02->Fill(teclon[k1],teclat[k2],vdtec_h02_on[k1+73*k2]);
       h03->Fill(teclon[k1],teclat[k2],vdtec_h02_on[k1+73*k2]);

       h12->Fill(teclon[k1],teclat[k2],vdtec_h24_on[k1+73*k2]);
       h13->Fill(teclon[k1],teclat[k2],vdtec_h24_on[k1+73*k2]);

       h22->Fill(teclon[k1],teclat[k2],vdtec_h46_on[k1+73*k2]);
       h23->Fill(teclon[k1],teclat[k2],vdtec_h46_on[k1+73*k2]);
}}   


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TH2D *h30 = new TH2D("h30"," ",73,-182.5,182.5,71,-88.75,88.75);
TH2D *h31 = new TH2D("h31"," ",73,-182.5,182.5,71,-88.75,88.75);
TH2D *h32 = new TH2D("h32"," ",73,-182.5,182.5,71,-88.75,88.75);
TH2D *h33 = new TH2D("h33"," ",73,-182.5,182.5,71,-88.75,88.75);
TH2D *h40 = new TH2D("h40"," ",73,-182.5,182.5,71,-88.75,88.75);
TH2D *h41 = new TH2D("h41"," ",73,-182.5,182.5,71,-88.75,88.75);
TH2D *h42 = new TH2D("h42"," ",73,-182.5,182.5,71,-88.75,88.75);
TH2D *h43 = new TH2D("h43"," ",73,-182.5,182.5,71,-88.75,88.75);
TH2D *h50 = new TH2D("h50"," ",73,-182.5,182.5,71,-88.75,88.75);
TH2D *h51 = new TH2D("h51"," ",73,-182.5,182.5,71,-88.75,88.75);
TH2D *h52 = new TH2D("h52"," ",73,-182.5,182.5,71,-88.75,88.75);
TH2D *h53 = new TH2D("h53"," ",73,-182.5,182.5,71,-88.75,88.75);

for(int nperiod=0;nperiod<3;nperiod++){ //start loop for three periods 

//Global fixed LT at 11:00-13:00LT(4:00-6:00UT) @107.419E  time zone (+7)=lt_time_bin[18]
float lt_lon_bin[25][2]; //two bin limits
float lt_time_bin[25];
float lt_lon_ind[25][3];
float arr_h1113lt_on[25][31][2];

for(int i=0;i<25;i++){
  lt_time_bin[i]=-12+i;
  if(i==0){lt_lon_bin[i][0]=-180.;lt_lon_bin[i][1]=-172.5;lt_lon_ind[i][0]=0;lt_lon_ind[i][1]=1;lt_lon_ind[i][2]=1;}
  if(i==24){lt_lon_bin[i][0]=172.5;lt_lon_bin[i][1]=180.;lt_lon_ind[i][0]=71;lt_lon_ind[i][1]=72;lt_lon_ind[i][2]=72;}
  if(i!=0 && i!=24){
      lt_lon_bin[i][0]=-172.5+15.*(i-1);lt_lon_bin[i][1]=-172.5+15.*i;
      lt_lon_ind[i][0]=3*(i-1)+2;lt_lon_ind[i][1]=3*(i-1)+3;lt_lon_ind[i][2]=3*(i-1)+4;
  }

  for(int j=0;j<31;j++){
	if(nperiod==0){arr_h1113lt_on[i][j][0]=468+12*(j+1)-0.5*(i+1)+2;arr_h1113lt_on[i][j][1]=468+12*(j+1)-0.5*(i+1)+3;} //LT for 02-04UT
	//if(nperiod==0){arr_h1113lt_on[i][j][0]=528+12*(j+1)-0.5*(i+1)+2;arr_h1113lt_on[i][j][1]=528+12*(j+1)-0.5*(i+1)+3;}
           

	if(nperiod==1){arr_h1113lt_on[i][j][0]=528+12*(j+1)-0.5*(i+1)+3;arr_h1113lt_on[i][j][1]=528+12*(j+1)-0.5*(i+1)+4;} //LT for 04-06UT
	//if(nperiod==1){arr_h1113lt_on[i][j][0]=552+12*(j+1)-0.5*(i+1)-1;arr_h1113lt_on[i][j][1]=552+12*(j+1)-0.5*(i+1);}
	//if(nperiod==1){arr_h1113lt_on[i][j][0]=552+12*(j+1)-0.5*(i+1)-1;arr_h1113lt_on[i][j][1]=552+12*(j+1)-0.5*(i+1);}
	
	if(nperiod==2){arr_h1113lt_on[i][j][0]=552+12*(j+1)-0.5*(i+1)+2;arr_h1113lt_on[i][j][1]=552+12*(j+1)-0.5*(i+1)+3;} //LT for 12-14UT
	//if(nperiod==2){arr_h1113lt_on[i][j][0]=552+12*(j+1)-0.5*(i+1)+7;arr_h1113lt_on[i][j][1]=552+12*(j+1)-0.5*(i+1)+8;}
	//if(nperiod==2){arr_h1113lt_on[i][j][0]=636+12*(j+1)-0.5*(i+1)+7+1;arr_h1113lt_on[i][j][1]=636+12*(j+1)-0.5*(i+1)+8+1;}
  }

  //cout<<arr_h1113lt_on[i][30][0]<<";"<<arr_h1113lt_on[i][30][1]<<endl;
}
 
float xvtec_max_on[5183],xvtec_q1_on[5183],xvtec_median_on[5183],xvtec_q2_on[5183],xvdtec_on[5183];
for(int i=0;i<5183;i++){xvtec_max_on[i]=0.;xvtec_q1_on[i]=0.;xvtec_median_on[i]=0.;xvtec_q2_on[i]=0.;xvdtec_on[i]=0.;}

//start the loop for 25 time zones
for(int nt=0;nt<25;nt++){  //25 time zones
//read TEC data again
//GPS-TEC: NWC ON; EQ7.7-(107.419,-9.284)-20.0km-south of Java-06/07/17-08:19:26UTC
ifstream in4;
in4.open("/home/wang/Projects/EQ/Yahama/GPS_TEC2HR_IGS_Yahama_20190430_0724.txt");
//EQ:06-07-17 08:19:26UTC; +75:05-03; +60:05-18; +45:06-02; +30:06-17; -10:07-27     
int xnrecord2=62; ////06/17(540)-07/17(912)--31d  
float xtecyear_on[xnrecord2],xtecsec_on[xnrecord2],xteccode_on[xnrecord2][213];  //LON:71*3 for each time zone
for(int i=0;i<xnrecord2;i++){
for(int j=0;j<213;j++){
   xteccode_on[i][j]=0.;
}}

int xarr_ltt_on[213];
for(int i=0;i<71;i++){
   xarr_ltt_on[3*i]=lt_lon_ind[nt][0]+73*i;xarr_ltt_on[3*i+1]=lt_lon_ind[nt][1]+73*i;xarr_ltt_on[3*i+2]=lt_lon_ind[nt][2]+73*i;
}
//cout<<xarr_ltt_on[0]<<";"<<xarr_ltt_on[2]<<endl;

//char xarr2[90000];
//int xm2=1;
int xnrd2=0;
int xnrd2_temp0=0;
float xvtecyeartemp_on=0., xvtecsectemp_on=0., xvtectemp_on=0.;

 /*while(xm2<119)
{
  in4.getline(xarr2,90000);
  xm2++;
}
 in4.getline(xarr2,90000);*/ 

while(in4.good())
{
 int xnrd2_temp1=0;

 for(int k1=0;k1<31;k1++){  //time limits for each time zone
   if(xnrd2>=arr_h1113lt_on[nt][k1][0] && xnrd2<=arr_h1113lt_on[nt][k1][1]){xnrd2_temp1=xnrd2;break;}
   else{xnrd2_temp1=0;}
 }

 if(xnrd2_temp1!=0){
    in4>>xtecyear_on[xnrd2_temp0]>>xtecsec_on[xnrd2_temp0];
    int xnrd2_temp2=0;

    int gap3=0;
    for(int k2=0;k2<5183;k2++){
       int xnrd2_temp_flag=0;

       for(int k3=0;k3<213;k3++){
          if(k2==xarr_ltt_on[k3]){
             xnrd2_temp_flag=1;
             if(nt==0 || nt==24){gap3++;}
             break;
          }
       } 

       if(xnrd2_temp_flag!=0){
          if(gap3==2 && (nt==0 || nt==24)){
             in4>>xteccode_on[xnrd2_temp0][xnrd2_temp2];xnrd2_temp2++;
             xteccode_on[xnrd2_temp0][xnrd2_temp2]=xteccode_on[xnrd2_temp0][xnrd2_temp2-1];xnrd2_temp2++;gap3=0;
          }
          else{in4>>xteccode_on[xnrd2_temp0][xnrd2_temp2];xnrd2_temp2++;}
       }
       else{in4>>xvtectemp_on;}
    }
    if(nt%2==0){xnrd2_temp0=xnrd2_temp0+2;}
    else{xnrd2_temp0++;}
 }
 else{
    in4>>xvtecyeartemp_on>>xvtecsectemp_on;
    for(int k4=0;k4<5183;k4++){in4>>xvtectemp_on;}
 }

 xnrd2++;
} 
in4.close();

//if(nt==0){cout<<xteccode_on[0][6]<<";"<<xteccode_on[0][7]<<";"<<xteccode_on[0][8]<<";"<<endl;}

float xteccode_c2h_on[31][213];
for(int i=0;i<31;i++){
   for(int j=0;j<213;j++){
      xteccode_c2h_on[i][j]=0.;
      if(nt%2==0){xteccode_c2h_on[i][j]=xteccode_on[2*i][j]+xteccode_on[2*i+1][j];} 
      else{xteccode_c2h_on[i][j]=0.5*(xteccode_on[2*i][j]+xteccode_on[2*i+1][j]);}
   }
}

for(int l1=0;l1<213;l1++){
   float xtec_quartile_temp_on[30];

   for(int l2=0;l2<30;l2++){
      xtec_quartile_temp_on[l2]=0.;
      xtec_quartile_temp_on[l2]=xteccode_c2h_on[l2][l1];
   }
   
   int mapind0=xarr_ltt_on[l1];
     	   if(nperiod==0){ xvtec_max_on[mapind0]=xteccode_c2h_on[28][l1];}
	   if(nperiod==1){ xvtec_max_on[mapind0]=xteccode_c2h_on[28][l1];}
	   if(nperiod==2){ xvtec_max_on[mapind0]=xteccode_c2h_on[28][l1];}
   xvtec_q1_on[mapind0]=map_calc_quartile(111,l1,30,xtec_quartile_temp_on);
   xvtec_median_on[mapind0]=map_calc_quartile(112,l1,30,xtec_quartile_temp_on); 
   xvtec_q2_on[mapind0]=map_calc_quartile(113,l1,30,xtec_quartile_temp_on);  

   if(nperiod==0||2){xvdtec_on[mapind0]=(xvtec_max_on[mapind0]-xvtec_median_on[mapind0])/(xvtec_q2_on[mapind0]-xvtec_q1_on[mapind0]);}
   if(nperiod==1){xvdtec_on[mapind0]=(xvtec_max_on[mapind0]-xvtec_median_on[mapind0])/(xvtec_q2_on[mapind0]-xvtec_q1_on[mapind0]);}
   //cout<<xvdtec_on[mapind0]<<endl;  
   if(xvdtec_on[mapind0]>-2. && xvdtec_on[mapind0]<2.){xvdtec_on[mapind0]=-7.;}
}

}//end the loop for 25 time zones

for(int k2=0;k2<71;k2++){  
for(int k1=0;k1<73;k1++){
   if(nperiod==0){
      h30->Fill(teclon[k1],teclat[k2],xvtec_max_on[k1+73*k2]);
      h31->Fill(teclon[k1],teclat[k2],xvtec_median_on[k1+73*k2]);
      h32->Fill(teclon[k1],teclat[k2],xvdtec_on[k1+73*k2]);
      h33->Fill(teclon[k1],teclat[k2],xvdtec_on[k1+73*k2]);
   //   if(teclon[k1]>95. && teclon[k1]<170. && teclat[k2]>-22. && teclat[k2]<18.){cout<<"1st "<<(xvdtec_on[k1+73*k2] - 2) /2 * 100<<"%"<<endl;}
   }
   if(nperiod==1){
      h40->Fill(teclon[k1],teclat[k2],xvtec_max_on[k1+73*k2]);
      h41->Fill(teclon[k1],teclat[k2],xvtec_median_on[k1+73*k2]);
      h42->Fill(teclon[k1],teclat[k2],xvdtec_on[k1+73*k2]);
      h43->Fill(teclon[k1],teclat[k2],xvdtec_on[k1+73*k2]);
    //  if(teclon[k1]>95. && teclon[k1]<170. && teclat[k2]>-22. && teclat[k2]<18.){cout<<"2nd "<<(xvdtec_on[k1+73*k2] - 2) /2 * 100<<"%"<<endl;}
   }
   if(nperiod==2){
      h50->Fill(teclon[k1],teclat[k2],xvtec_max_on[k1+73*k2]);
      h51->Fill(teclon[k1],teclat[k2],xvtec_median_on[k1+73*k2]);
      h52->Fill(teclon[k1],teclat[k2],xvdtec_on[k1+73*k2]);
      h53->Fill(teclon[k1],teclat[k2],xvdtec_on[k1+73*k2]);
    //  if(teclon[k1]>95. && teclon[k1]<170. && teclat[k2]>-22. && teclat[k2]<18.){cout<<"3rd "<<(xvdtec_on[k1+73*k2] - 2) /2 * 100<<"%"<<endl;}    
   }
}}
   
}//end loop for three periods
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


//plotting================================================================
gStyle->SetOptStat(000);
gStyle->SetPalette(55);
//gStyle->SetTitleFontSize(0.0);
gStyle->SetErrorX(0);

TApplication theApp("App", &argc, argv);

TCanvas *c0 = new TCanvas("c0","TEC",960,550);
c0->Connect("Closed()", "TApplication", &theApp, "Terminate()");

TPad *pad01 = new TPad("pad01","pad01",0.02,0.66,0.35,0.97);
pad01->Draw();
TPad *pad02 = new TPad("pad02","pad02",0.02,0.35,0.35,0.66);
pad02->Draw();
TPad *pad03 = new TPad("pad03","pad03",0.02,0.02,0.35,0.35);
pad03->Draw();
TPad *pad04 = new TPad("pad04","pad04",0.35,0.66,0.67,0.97);
pad04->Draw();
TPad *pad05 = new TPad("pad05","pad05",0.35,0.35,0.67,0.66);
pad05->Draw();
TPad *pad06 = new TPad("pad06","pad06",0.35,0.02,0.67,0.35);
pad06->Draw();
TPad *pad07 = new TPad("pad07","pad07",0.67,0.66,0.99,0.97);
pad07->Draw();
TPad *pad08 = new TPad("pad08","pad08",0.67,0.35,0.99,0.66);
pad08->Draw();
TPad *pad09 = new TPad("pad09","pad09",0.67,0.02,0.99,0.35);
pad09->Draw();


//++++++++++++++++
pad01->cd();
pad01->SetTopMargin(0.08);
pad01->SetBottomMargin(0.11);
pad01->SetLeftMargin(0.10);
pad01->SetRightMargin(0.09);
h30->Draw("colz ");
h30->GetXaxis()->SetNdivisions(210);
h30->GetYaxis()->SetNdivisions(410);
h30->GetXaxis()->SetTickLength(0.025);
h30->GetYaxis()->SetTickLength(0.015);
h30->GetXaxis()->SetLabelSize(0.05);
h30->GetYaxis()->SetLabelSize(0.05);
h30->GetXaxis()->SetLabelOffset(0.02);
h30->GetXaxis()->CenterTitle();
h30->GetYaxis()->CenterTitle();
h30->GetYaxis()->SetDecimals(kTRUE);
h30->GetXaxis()->SetRangeUser(-180., 180.);
h30->GetYaxis()->SetRangeUser(-81.25, 86.25);

hm->Draw("CONT2 SAME");
hm->SetLineColor(kGray+3);
hm->SetLineWidth(0.1);

g02->SetLineColor(kRed);
g02->Draw("SAME");

g00->SetMarkerStyle(29);
g00->SetMarkerSize(1.2);
g00->SetMarkerColor(kRed);
g00->Draw("P SAME");

gPad->Update();
TPaletteAxis *palette30 = (TPaletteAxis*)h30->GetListOfFunctions()->FindObject("palette");
if(palette30) {
   palette30->SetX1NDC(0.922); 
   palette30->SetX2NDC(0.935); 
   gPad->Modified(); 
}

h30->GetZaxis()->SetRangeUser(0,30);
h30->GetZaxis()->SetNdivisions(206);
h30->GetZaxis()->SetLabelSize(0.049);
h30->GetZaxis()->SetLabelOffset(0.008);
h30->GetZaxis()->SetTickLength(0.009);
h30->GetZaxis()->CenterTitle();

g_mag_eq->SetLineColor(kBlack);
g_mag_eq->Draw("C SAME");

TLatex *t300 = new TLatex();
   t300->SetTextSize(0.064);
   t300->DrawLatex(-172,90,"2019/07/06 global fixed 14:00-16:00 LT LLT map");

TLatex *t301 = new TLatex();
   t301->SetTextSize(0.07);
   t301->SetTextAngle(90);
   t301->DrawLatex(-213.,-35.,"Latitude (#circ)");

TLatex *t302 = new TLatex();
   t302->SetTextSize(0.05);
   t302->SetTextAngle(90);
   t302->DrawLatex(220.,-30.,"TEC (TECU)");

TLatex *t303 = new TLatex();
   t303->SetTextSize(0.07);
   t303->DrawLatex(-225.,90.,"(a)");

pad02->cd();
pad02->SetTopMargin(0.08);
pad02->SetBottomMargin(0.11);
pad02->SetLeftMargin(0.10);
pad02->SetRightMargin(0.09);
h31->Draw("colz ");
h31->GetXaxis()->SetNdivisions(210);
h31->GetYaxis()->SetNdivisions(410);
h31->GetXaxis()->SetTickLength(0.025);
h31->GetYaxis()->SetTickLength(0.015);
h31->GetXaxis()->SetLabelSize(0.05);
h31->GetYaxis()->SetLabelSize(0.05);
h31->GetXaxis()->SetLabelOffset(0.02);
h31->GetXaxis()->CenterTitle();
h31->GetYaxis()->CenterTitle();
h31->GetYaxis()->SetDecimals(kTRUE);
h31->GetXaxis()->SetRangeUser(-180., 180.);
h31->GetYaxis()->SetRangeUser(-81.25, 86.25);

hm->Draw("CONT2 SAME");
hm->SetLineColor(kGray+3);
hm->SetLineWidth(0.1);

g02->SetLineColor(kRed);
g02->Draw("SAME");

g00->SetMarkerStyle(29);
g00->SetMarkerSize(1.2);
g00->SetMarkerColor(kRed);
g00->Draw("P SAME");

gPad->Update();
TPaletteAxis *palette31 = (TPaletteAxis*)h31->GetListOfFunctions()->FindObject("palette");
if(palette31) {
   palette31->SetX1NDC(0.922); 
   palette31->SetX2NDC(0.935); 
   gPad->Modified(); 
}

h31->GetZaxis()->SetRangeUser(0,30);
h31->GetZaxis()->SetNdivisions(206);
h31->GetZaxis()->SetLabelSize(0.049);
h31->GetZaxis()->SetLabelOffset(0.008);
h31->GetZaxis()->SetTickLength(0.009);
h31->GetZaxis()->CenterTitle();

g_mag_eq->SetLineColor(kBlack);
g_mag_eq->Draw("C SAME");

TLatex *t310 = new TLatex();
   t310->SetTextSize(0.064);
   t310->DrawLatex(-190,90,"2019/06/07-07/07 14:00-16:00 LT 30-day median LLT map");

TLatex *t311 = new TLatex();
   t311->SetTextSize(0.07);
   t311->SetTextAngle(90);
   t311->DrawLatex(-213.,-35.,"Latitude (#circ)");

TLatex *t312 = new TLatex();
   t312->SetTextSize(0.05);
   t312->SetTextAngle(90);
   t312->DrawLatex(220.,-30.,"TEC (TECU)");

TLatex *t313 = new TLatex();
   t313->SetTextSize(0.07);
   t313->DrawLatex(-225.,90.,"(b)");

pad03->cd();
pad03->SetTopMargin(0.08);
pad03->SetBottomMargin(0.14);
pad03->SetLeftMargin(0.10);
pad03->SetRightMargin(0.09);
h32->Draw("colz ");
h32->GetXaxis()->SetNdivisions(210);
h32->GetYaxis()->SetNdivisions(410);
h32->GetXaxis()->SetTickLength(0.025);
h32->GetYaxis()->SetTickLength(0.015);
h32->GetXaxis()->SetLabelSize(0.05);
h32->GetYaxis()->SetLabelSize(0.05);
h32->GetXaxis()->SetLabelOffset(0.02);
h32->GetXaxis()->CenterTitle();
h32->GetYaxis()->CenterTitle();
h32->GetYaxis()->SetDecimals(kTRUE);
h32->GetXaxis()->SetRangeUser(95., 170.);
h32->GetYaxis()->SetRangeUser(-22., 18.);
//h32->GetXaxis()->SetRangeUser(-180., 180.);
//h32->GetYaxis()->SetRangeUser(-81.25, 86.25);


hm->Draw("CONT2 SAME");
hm->SetLineColor(kGray+3);
hm->SetLineWidth(0.1);

g00->SetMarkerStyle(29);
g00->SetMarkerSize(1.2);
g00->SetMarkerColor(kRed);
g00->Draw("P SAME");

g01->SetMarkerStyle(5);
g01->SetMarkerSize(1.2);
g01->SetMarkerColor(kRed);
//g01->Draw("P SAME");

//g02->SetLineColor(kRed);
//g02->Draw("SAME");

//g_mag_eq->SetLineColor(kBlack);
//g_mag_eq->Draw("C SAME");

TEllipse *el02 = new TEllipse(128.17,-0.52,12.,12.);
   el02->SetFillStyle(0);
   el02->SetLineStyle(2);
   el02->SetLineColor(kRed);
   el02->Draw();

gPad->Update();
TPaletteAxis *palette32 = (TPaletteAxis*)h32->GetListOfFunctions()->FindObject("palette");
if(palette32) {
   palette32->SetX1NDC(0.922); 
   palette32->SetX2NDC(0.935); 
   gPad->Modified(); 
}

h32->GetZaxis()->SetRangeUser(0.,3.);
//h32->GetZaxis()->SetRangeUser(-2.,6.);
//h32->GetZaxis()->SetRangeUser(-3.,7.);
h32->GetZaxis()->SetNdivisions(204);
h32->GetZaxis()->SetLabelSize(0.049);
h32->GetZaxis()->SetLabelOffset(0.008);
h32->GetZaxis()->SetTickLength(0.009);
h32->GetZaxis()->CenterTitle();

TLatex *t320 = new TLatex();
   t320->SetTextSize(0.064);
   t320->DrawLatex(90.,20.,"2019/07/06 14:00-16:00 LT difference LLT (|#DeltaTEC|>2.0) map");
//  t320->DrawLatex(-150,90.,"2019/07/06 14:00-16:00 LT difference LLT map");
   
TLatex *t321 = new TLatex();
   t321->SetTextSize(0.07);
   t321->SetTextAngle(90);
   t321->DrawLatex(86.5,-12.,"Latitude (#circ)");
 //  t321->DrawLatex(-213,-35.,"Latitude (#circ)");  

TLatex *t322 = new TLatex();
   t322->SetTextSize(0.066);
   t322->DrawLatex(122.5,-30.5,"Longitude (#circ)");
  //t322->DrawLatex(-50,-105,"Longitude (#circ)");

TLatex *t323 = new TLatex();
   t323->SetTextSize(0.05);
   t323->SetTextAngle(90);
   t323->DrawLatex(180.,-5.5,"#DeltaTEC");
 // t323->DrawLatex(220.,-30.,"TEC(TECU)");

TLatex *t324 = new TLatex();
   t324->SetTextSize(0.075);
   t324->DrawLatex(82.5,20.,"(d)");
   //t324->DrawLatex(-225.,90.,"(c)");

//++++++++++++++++
pad04->cd();
pad04->SetTopMargin(0.08);
pad04->SetBottomMargin(0.11);
pad04->SetLeftMargin(0.08);
pad04->SetRightMargin(0.09);
h40->Draw("colz ");
h40->GetXaxis()->SetNdivisions(210);
h40->GetYaxis()->SetNdivisions(410);
h40->GetXaxis()->SetTickLength(0.025);
h40->GetYaxis()->SetTickLength(0.015);
h40->GetXaxis()->SetLabelSize(0.05);
h40->GetYaxis()->SetLabelSize(0.05);
h40->GetXaxis()->SetLabelOffset(0.02);
h40->GetXaxis()->CenterTitle();
h40->GetYaxis()->CenterTitle();
h40->GetYaxis()->SetDecimals(kTRUE);
h40->GetXaxis()->SetRangeUser(-180., 180.);
h40->GetYaxis()->SetRangeUser(-81.25, 86.25);

hm->Draw("CONT2 SAME");
hm->SetLineColor(kGray+3);
hm->SetLineWidth(0.1);

g02->SetLineColor(kRed);
g02->Draw("SAME");

g00->SetMarkerStyle(29);
g00->SetMarkerSize(1.2);
g00->SetMarkerColor(kRed);
g00->Draw("P SAME");

gPad->Update();
TPaletteAxis *palette40 = (TPaletteAxis*)h40->GetListOfFunctions()->FindObject("palette");
if(palette40) {
   palette40->SetX1NDC(0.922); 
   palette40->SetX2NDC(0.935); 
   gPad->Modified(); 
}

h40->GetZaxis()->SetRangeUser(0,30);
h40->GetZaxis()->SetNdivisions(206);
h40->GetZaxis()->SetLabelSize(0.049);
h40->GetZaxis()->SetLabelOffset(0.008);
h40->GetZaxis()->SetTickLength(0.009);
h40->GetZaxis()->CenterTitle();

g_mag_eq->SetLineColor(kBlack);
g_mag_eq->Draw("C SAME");

TLatex *t400 = new TLatex();
   t400->SetTextSize(0.064);
   t400->DrawLatex(-173,90,"2019/07/11 global fixed 16:00-18:00 LT LLT map");

TLatex *t401 = new TLatex();
   t401->SetTextSize(0.05);
   t401->SetTextAngle(90);
   t401->DrawLatex(220.,-30.,"TEC (TECU)");

pad05->cd();
pad05->SetTopMargin(0.08);
pad05->SetBottomMargin(0.11);
pad05->SetLeftMargin(0.08);
pad05->SetRightMargin(0.09);
h41->Draw("colz ");
h41->GetXaxis()->SetNdivisions(210);
h41->GetYaxis()->SetNdivisions(410);
h41->GetXaxis()->SetTickLength(0.025);
h41->GetYaxis()->SetTickLength(0.015);
h41->GetXaxis()->SetLabelSize(0.05);
h41->GetYaxis()->SetLabelSize(0.05);
h41->GetXaxis()->SetLabelOffset(0.02);
h41->GetXaxis()->CenterTitle();
h41->GetYaxis()->CenterTitle();
h41->GetYaxis()->SetDecimals(kTRUE);
h41->GetXaxis()->SetRangeUser(-180., 180.);
h41->GetYaxis()->SetRangeUser(-81.25, 86.25);

hm->Draw("CONT2 SAME");
hm->SetLineColor(kGray+3);
hm->SetLineWidth(0.1);

g02->SetLineColor(kRed);
g02->Draw("SAME");

g00->SetMarkerStyle(29);
g00->SetMarkerSize(1.2);
g00->SetMarkerColor(kRed);
g00->Draw("P SAME");

gPad->Update();
TPaletteAxis *palette41 = (TPaletteAxis*)h41->GetListOfFunctions()->FindObject("palette");
if(palette41) {
   palette41->SetX1NDC(0.922); 
   palette41->SetX2NDC(0.935); 
   gPad->Modified(); 
}

h41->GetZaxis()->SetRangeUser(0,30);
h41->GetZaxis()->SetNdivisions(206);
h41->GetZaxis()->SetLabelSize(0.049);
h41->GetZaxis()->SetLabelOffset(0.008);
h41->GetZaxis()->SetTickLength(0.009);
h41->GetZaxis()->CenterTitle();

g_mag_eq->SetLineColor(kBlack);
g_mag_eq->Draw("C SAME");

TLatex *t410 = new TLatex();
   t410->SetTextSize(0.064);
   t410->DrawLatex(-191,90,"2019/06/12-07/12 16:00-18:00 LT 30-day median LLT map");

TLatex *t411 = new TLatex();
   t411->SetTextSize(0.05);
   t411->SetTextAngle(90);
   t411->DrawLatex(220.,-30.,"TEC (TECU)");

pad06->cd();
pad06->SetTopMargin(0.08);
pad06->SetBottomMargin(0.14);
pad06->SetLeftMargin(0.08);
pad06->SetRightMargin(0.09);
h42->Draw("colz ");
h42->GetXaxis()->SetNdivisions(210);
h42->GetYaxis()->SetNdivisions(410);
h42->GetXaxis()->SetTickLength(0.025);
h42->GetYaxis()->SetTickLength(0.015);
h42->GetXaxis()->SetLabelSize(0.05);
h42->GetYaxis()->SetLabelSize(0.05);
h42->GetXaxis()->SetLabelOffset(0.02);
h42->GetXaxis()->CenterTitle();
h42->GetYaxis()->CenterTitle();
h42->GetYaxis()->SetDecimals(kTRUE);
h42->GetXaxis()->SetRangeUser(95., 170.);
h42->GetYaxis()->SetRangeUser(-22., 18.);
//h32->GetXaxis()->SetRangeUser(-180., 180.);
//h32->GetYaxis()->SetRangeUser(-81.25, 86.25);

hm->Draw("CONT2 SAME");
hm->SetLineColor(kGray+3);
hm->SetLineWidth(0.1);

g00->SetMarkerStyle(29);
g00->SetMarkerSize(1.2);
g00->SetMarkerColor(kRed);
g00->Draw("P SAME");

//g02->SetLineColor(kRed);
//g02->Draw("SAME");

//g_mag_eq->SetLineColor(kBlack);
//g_mag_eq->Draw("C SAME");

TEllipse *el12 = new TEllipse(128.17,-0.52,12.,12.);
   el12->SetFillStyle(0);
   el12->SetLineStyle(2);
   el12->SetLineColor(kRed);
   el12->Draw();

gPad->Update();
TPaletteAxis *palette42 = (TPaletteAxis*)h42->GetListOfFunctions()->FindObject("palette");
if(palette42) {
   palette42->SetX1NDC(0.922); 
   palette42->SetX2NDC(0.935); 
   gPad->Modified(); 
}

h42->GetZaxis()->SetRangeUser(0.,3.);
//h42->GetZaxis()->SetRangeUser(-2.,6.);
//h42->GetZaxis()->SetRangeUser(-3.,7.);
h42->GetZaxis()->SetNdivisions(204);
h42->GetZaxis()->SetLabelSize(0.049);
h42->GetZaxis()->SetLabelOffset(0.008);
h42->GetZaxis()->SetTickLength(0.009);
h42->GetZaxis()->CenterTitle();

TLatex *t420 = new TLatex();
   t420->SetTextSize(0.064);
   t420->DrawLatex(90.,20.,"2019/07/11 16:00-18:00 LT difference LLT (|#DeltaTEC|>2.0) map");
//   t420->DrawLatex(-150,90.,"2019/07/11 16:00-18:00LT difference LLT map");

TLatex *t422 = new TLatex();
   t422->SetTextSize(0.066);
   t422->DrawLatex(122.5,-30.5,"Longitude (#circ)");
  // t422->DrawLatex(-50,-115,"Longitude (#circ)");

TLatex *t423 = new TLatex();
   t423->SetTextSize(0.05);
   t423->SetTextAngle(90);
   t423->DrawLatex(180.,-5.5,"#DeltaTEC");
  // t423->DrawLatex(220.,-30.,"TEC(TECU)");

//++++++++++++++++
pad07->cd();
pad07->SetTopMargin(0.08);
pad07->SetBottomMargin(0.11);
pad07->SetLeftMargin(0.08);
pad07->SetRightMargin(0.09);
h50->Draw("colz ");
h50->GetXaxis()->SetNdivisions(210);
h50->GetYaxis()->SetNdivisions(410);
h50->GetXaxis()->SetTickLength(0.025);
h50->GetYaxis()->SetTickLength(0.015);
h50->GetXaxis()->SetLabelSize(0.05);
h50->GetYaxis()->SetLabelSize(0.05);
h50->GetXaxis()->SetLabelOffset(0.02);
h50->GetXaxis()->CenterTitle();
h50->GetYaxis()->CenterTitle();
h50->GetYaxis()->SetDecimals(kTRUE);
h50->GetXaxis()->SetRangeUser(-180., 180.);
h50->GetYaxis()->SetRangeUser(-81.25, 86.25);

hm->Draw("CONT2 SAME");
hm->SetLineColor(kGray+3);
hm->SetLineWidth(0.1);

g02->SetLineColor(kRed);
g02->Draw("SAME");

g00->SetMarkerStyle(29);
g00->SetMarkerSize(1.2);
g00->SetMarkerColor(kRed);
g00->Draw("P SAME");

gPad->Update();
TPaletteAxis *palette50 = (TPaletteAxis*)h50->GetListOfFunctions()->FindObject("palette");
if(palette50) {
   palette50->SetX1NDC(0.922); 
   palette50->SetX2NDC(0.935); 
   gPad->Modified(); 
}

h50->GetZaxis()->SetRangeUser(0,30);
h50->GetZaxis()->SetNdivisions(206);
h50->GetZaxis()->SetLabelSize(0.049);
h50->GetZaxis()->SetLabelOffset(0.008);
h50->GetZaxis()->SetTickLength(0.009);
h50->GetZaxis()->CenterTitle();

g_mag_eq->SetLineColor(kBlack);
g_mag_eq->Draw("C SAME");

TLatex *t500 = new TLatex();
   t500->SetTextSize(0.064);
   t500->DrawLatex(-173,90,"2019/07/13 global fixed 14:00-16:00 LT LLT map");
//   t500->DrawLatex(-173,90,"2019/07/20 global fixed 00:00-02:00LT LLT map");

TLatex *t501 = new TLatex();
   t501->SetTextSize(0.05);
   t501->SetTextAngle(90);
   t501->DrawLatex(220.,-30.,"TEC (TECU)");

pad08->cd();
pad08->SetTopMargin(0.08);
pad08->SetBottomMargin(0.11);
pad08->SetLeftMargin(0.08);
pad08->SetRightMargin(0.09);
h51->Draw("colz ");
h51->GetXaxis()->SetNdivisions(210);
h51->GetYaxis()->SetNdivisions(410);
h51->GetXaxis()->SetTickLength(0.025);
h51->GetYaxis()->SetTickLength(0.015);
h51->GetXaxis()->SetLabelSize(0.05);
h51->GetYaxis()->SetLabelSize(0.05);
h51->GetXaxis()->SetLabelOffset(0.02);
h51->GetXaxis()->CenterTitle();
h51->GetYaxis()->CenterTitle();
h51->GetYaxis()->SetDecimals(kTRUE);
h51->GetXaxis()->SetRangeUser(-180., 180.);
h51->GetYaxis()->SetRangeUser(-81.25, 86.25);

hm->Draw("CONT2 SAME");
hm->SetLineColor(kGray+3);
hm->SetLineWidth(0.1);

g02->SetLineColor(kRed);
g02->Draw("SAME");

g00->SetMarkerStyle(29);
g00->SetMarkerSize(1.2);
g00->SetMarkerColor(kRed);
g00->Draw("P SAME");

gPad->Update();
TPaletteAxis *palette51 = (TPaletteAxis*)h51->GetListOfFunctions()->FindObject("palette");
if(palette51) {
   palette51->SetX1NDC(0.922); 
   palette51->SetX2NDC(0.935); 
   gPad->Modified(); 
}

h51->GetZaxis()->SetRangeUser(0,30);
h51->GetZaxis()->SetNdivisions(206);
h51->GetZaxis()->SetLabelSize(0.049);
h51->GetZaxis()->SetLabelOffset(0.008);
h51->GetZaxis()->SetTickLength(0.009);
h51->GetZaxis()->CenterTitle();

g_mag_eq->SetLineColor(kBlack);
g_mag_eq->Draw("C SAME");

TLatex *t510 = new TLatex();
   t510->SetTextSize(0.064);
   t510->DrawLatex(-191,90,"2019/06/14-07/14 14:00-16:00 LT 30-day median LLT map");
 //  t510->DrawLatex(-191,90,"2019/06/21-07/21 00:00-02:00LT 30-day median LLT map");

TLatex *t511 = new TLatex();
   t511->SetTextSize(0.05);
   t511->SetTextAngle(90);
   t511->DrawLatex(220.,-30.,"TEC (TECU)");

pad09->cd();
pad09->SetTopMargin(0.08);
pad09->SetBottomMargin(0.14);
pad09->SetLeftMargin(0.08);
pad09->SetRightMargin(0.09);
h52->Draw("colz ");
h52->GetXaxis()->SetNdivisions(210);
h52->GetYaxis()->SetNdivisions(410);
h52->GetXaxis()->SetTickLength(0.025);
h52->GetYaxis()->SetTickLength(0.015);
h52->GetXaxis()->SetLabelSize(0.05);
h52->GetYaxis()->SetLabelSize(0.05);
h52->GetXaxis()->SetLabelOffset(0.02);
h52->GetXaxis()->CenterTitle();
h52->GetYaxis()->CenterTitle();
h52->GetYaxis()->SetDecimals(kTRUE);
h52->GetXaxis()->SetRangeUser(95.,170.);
h52->GetYaxis()->SetRangeUser(-22.,18.);
//h32->GetXaxis()->SetRangeUser(-180., 180.);
//h32->GetYaxis()->SetRangeUser(-81.25, 86.25);

hm->Draw("CONT2 SAME");
hm->SetLineColor(kGray+3);
hm->SetLineWidth(0.1);

g00->SetMarkerStyle(29);
g00->SetMarkerSize(1.2);
g00->SetMarkerColor(kRed);
g00->Draw("P SAME");

//g02->SetLineColor(kRed);
//g02->Draw("SAME");

//g_mag_eq->SetLineColor(kBlack);
//g_mag_eq->Draw("C SAME");

TEllipse *el22 = new TEllipse(128.17,-0.52,12.,12.);
   el22->SetFillStyle(0);
   el22->SetLineStyle(2);
   el22->SetLineColor(kRed);
   el22->Draw();

gPad->Update();
TPaletteAxis *palette52 = (TPaletteAxis*)h52->GetListOfFunctions()->FindObject("palette");
if(palette52) {
   palette52->SetX1NDC(0.922); 
   palette52->SetX2NDC(0.935); 
   gPad->Modified(); 
}

h52->GetZaxis()->SetRangeUser(0.,3.);
//h52->GetZaxis()->SetRangeUser(-3.,-1.);
//h52->GetZaxis()->SetRangeUser(-3.,7.);
//h52->GetZaxis()->SetRangeUser(-4.,4.);
h52->GetZaxis()->SetNdivisions(204);
h52->GetZaxis()->SetLabelSize(0.047);
h52->GetZaxis()->SetLabelOffset(0.008);
h52->GetZaxis()->SetTickLength(0.009);
h52->GetZaxis()->CenterTitle();
h52->GetZaxis()->SetDecimals(kTRUE);

TLatex *t520 = new TLatex();
   t520->SetTextSize(0.064);
   t520->DrawLatex(90.,20.,"2019/07/13 14:00-16:00 LT difference LLT (|#DeltaTEC|>2.0) map");
//   t520->DrawLatex(90.5,20.,"2019/07/20 00:00-02:00LT difference LLT (|DTEC|>2.0) map");
//   t520->DrawLatex(-150,90.,"2019/07/13 14:00-16:00LT difference LLT map");
  // t520->DrawLatex(-150,90.,"2019/07/20 00:00-02:00LT difference LLT map");

TLatex *t521 = new TLatex();
   t521->SetTextSize(0.066);
   t521->DrawLatex(122.5,-30.5,"Longitude (#circ)");
  // t521->DrawLatex(-50,-115,"Longitude (#circ)");


TLatex *t522 = new TLatex();
   t522->SetTextSize(0.047);
   t522->SetTextAngle(90);
   t522->DrawLatex(180.,-5.5,"#DeltaTEC");
 // t522->DrawLatex(220.,-30.,"TEC(TECU)");

theApp.Run();

}//End main()

//function for calculating the quartile Q1, Q2, Q3 and IQR
float map_calc_quartile(int ngridx, int vind, int nq, float *qteccode)
{
//sort original array in an ascending order
float atemp;
for(int i=0; i<nq-1; i++){
   for(int j=i+1;j<nq;j++){
      if(qteccode[i] > qteccode[j]){
          atemp=qteccode[i];
          qteccode[i]=qteccode[j];
          qteccode[j]=atemp;
      }
}}
//for(int i=0; i<nq; i++){cout<<qteccode[i]<<";";}
//cout<<endl;

float q1=0.,q3=0.,q2=0.; //Q2:median

float qtemp1=1.0*(nq+1)/4.;
int qint1=int(qtemp1);
float qmod1=qtemp1-qint1;
if(qmod1==0.){q1=qteccode[qint1-1];}
else
{q1=qteccode[qint1-1]+qmod1*(qteccode[qint1]-qteccode[qint1-1]);}

float qtemp2=2.0*(nq+1)/4.;
int qint2=int(qtemp2);
float qmod2=qtemp2-qint2;
if(qmod2==0.){q2=qteccode[qint2-1];}
else
{q2=qteccode[qint2-1]+qmod2*(qteccode[qint2]-qteccode[qint2-1]);}

float qtemp3=3.0*(nq+1)/4.;
int qint3=int(qtemp3);
float qmod3=qtemp3-qint3;
if(qmod3==0.){q3=qteccode[qint3-1];}
else
{q3=qteccode[qint3-1]+qmod3*(qteccode[qint3]-qteccode[qint3-1]);}

float iqr=0.;
iqr=q3-q1;

if(ngridx==0){
vtec_30d_quartile_h02_on[vind][0]=q1;
vtec_30d_quartile_h02_on[vind][1]=q2;
vtec_30d_quartile_h02_on[vind][2]=q3;
vtec_30d_quartile_h02_on[vind][3]=iqr;
}

if(ngridx==1){
vtec_30d_quartile_h24_on[vind][0]=q1;
vtec_30d_quartile_h24_on[vind][1]=q2;
vtec_30d_quartile_h24_on[vind][2]=q3;
vtec_30d_quartile_h24_on[vind][3]=iqr;
}

if(ngridx==2){
vtec_30d_quartile_h46_on[vind][0]=q1;
vtec_30d_quartile_h46_on[vind][1]=q2;
vtec_30d_quartile_h46_on[vind][2]=q3;
vtec_30d_quartile_h46_on[vind][3]=iqr;
}

if(ngridx==111){return q1;}
if(ngridx==112){return q2;}
if(ngridx==113){return q3;}

//return iqr;
}

//function of bilinear interpolation
float map_bilinear_interpolation(float *datapoint4)
{
//calculate the center point value of the grid
float tec_center=0.;
for(int k=0;k<4;k++){tec_center=tec_center+datapoint4[k];}
tec_center=0.25*tec_center;

return tec_center;
}
