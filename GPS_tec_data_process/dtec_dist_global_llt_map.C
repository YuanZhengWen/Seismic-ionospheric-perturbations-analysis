//===============================================================
//Demeter Case Code "dtec_dist_global_llt_map.C" 
//TAO DAN 2016.09.02 updated

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

//g++ dtec_dist_global_llt_map.C -O2 -Wall -fPIC -pthread -m64 -I/home/wang/CERNroot/root/root/include/ -L/home/wang/CERNroot/root/root/lib/ -lCore -lCint -lTree -lRIO  -lHist -lGpad -lGraf -o dtec_dist_global_llt_map
#define PI 3.141592653589793238463 

float vtec_30d_quartile_h24_on[5183][4];
float vtec_30d_quartile_h46_on[5183][4];
float vtec_30d_quartile_h1214_on[5183][4];

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

int arr_h24_on[62],arr_h46_on[62],arr_h1214_on[62];
for(int i=0;i<31;i++){arr_h24_on[2*i]=468+12*i+3;arr_h24_on[2*i+1]=468+12*i+4;}
//for(int i=0;i<31;i++){arr_h24_on[2*i]=528+12*i+4;arr_h24_on[2*i+1]=528+12*i+5;}

for(int i=0;i<31;i++){arr_h46_on[2*i]=528+12*i+4;arr_h46_on[2*i+1]=528+12*i+5;}
//for(int i=0;i<31;i++){arr_h46_on[2*i]=552+12*i+0;arr_h46_on[2*i+1]=552+12*i+1;}
//for(int i=0;i<31;i++){arr_h46_on[2*i]=552+12*i+3;arr_h46_on[2*i+1]=552+12*i+4;}

for(int i=0;i<31;i++){arr_h1214_on[2*i]=552+12*i+3;arr_h1214_on[2*i+1]=552+12*i+4;}
//for(int i=0;i<31;i++){arr_h1214_on[2*i]=636+12*i+8;arr_h1214_on[2*i+1]=636+12*i+9;}

float vtec_max_h24_on[5183],vtec_max_h46_on[5183],vtec_max_h1214_on[5183];
float vdtec_h24_on[5183],vdtec_h46_on[5183],vdtec_h1214_on[5183];

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//loop for three TEC anomaly grid int 07/15
for(int ngrid=0;ngrid<3;ngrid++){
ifstream in2;
in2.open("/home/wang/Projects/EQ/Yahama/GPS_TEC2HR_IGS_Yahama_20190430_0724.txt");
//EQ:06-07-17 08:19:26UTC; +75:05-03; +60:05-18; +45:06-02; +30:06-17; -10:07-27  
//dtec extreme maximum time ponit: 03:00,05:00,13:00 07-15; 30d before 07-17: 06-17
//06/17(540)-07/17(912)--31d
int nrecord2=62;   //03:00([2-4]UTC) in per 31 days
float tecyear_on[nrecord2],tecsec_on[nrecord2],teccode_on[nrecord2][5183];

for(int i=0;i<nrecord2;i++){
    tecyear_on[i]=0.;tecsec_on[i]=0.;
    for(int j=0;j<5183;j++){teccode_on[i][j]=0.;}
}

char arr2[90000];
int m2=1;
int nrd2=0;
int nrd2_temp0=0;
float vtecyeartemp_on=0., vtecsectemp_on=0., vtectemp_on=0.;

 /*while(m2<119)
{
  in2.getline(arr2,90000);
  m2++;
}
 in2.getline(arr2,90000);*/ 

while(in2.good())
{
 int nrd2_temp1=0;

   if(ngrid==0){
     for(int k1=0;k1<62;k1++){
        if(nrd2==arr_h24_on[k1]){nrd2_temp1=nrd2;break;}
        else{nrd2_temp1=0;}
     }
   }
   if(ngrid==1){
     for(int k1=0;k1<62;k1++){
        if(nrd2==arr_h46_on[k1]){nrd2_temp1=nrd2;break;}
        else{nrd2_temp1=0;}
     }
   }
   if(ngrid==2){
     for(int k1=0;k1<62;k1++){
        if(nrd2==arr_h1214_on[k1]){nrd2_temp1=nrd2;break;}
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

//TEC at 23:00([22-24]UTC) 07-04; the global map: grid 5040=72*70
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
   vtec_max_on[l1]=teccode_c2h_on[28][l1];
  // if(ngrid==0){ vtec_max_on[l1]=teccode_c2h_on[28][l1];}
   //if(ngrid==1){ vtec_max_on[l1]=teccode_c2h_on[28][l1];}
 //  if(ngrid==2){ vtec_max_on[l1]=teccode_c2h_on[28][l1];}
   for(int l2=0;l2<31;l2++){
      tec_quartile_temp_on[l2]=0.;
      tec_quartile_temp_on[l2]=teccode_c2h_on[l2][l1];
   }

   map_calc_quartile(ngrid,l1,30,tec_quartile_temp_on);
   //cout<<vtec_max_on[l1]<<";"<<vtec_30d_quartile_on[l1][1]<<";"<<vtec_30d_quartile_on[l1][3]<<";"<<endl;
}
//cout<<vtec_max_on[0][0]<<endl;

if(ngrid==0){
for(int l2=0;l2<5183;l2++){
vtec_max_h24_on[l2]=vtec_max_on[l2];  
if(vtec_30d_quartile_h24_on[l2][3]!=0.){vdtec_h24_on[l2]=(vtec_max_h24_on[l2]-vtec_30d_quartile_h24_on[l2][1])/vtec_30d_quartile_h24_on[l2][3];}
}
}
if(ngrid==1){
for(int l2=0;l2<5183;l2++){
vtec_max_h46_on[l2]=vtec_max_on[l2];
if(vtec_30d_quartile_h46_on[l2][3]!=0.){vdtec_h46_on[l2]=(vtec_max_h46_on[l2]-vtec_30d_quartile_h46_on[l2][1])/vtec_30d_quartile_h46_on[l2][3];}
}
}
if(ngrid==2){
for(int l2=0;l2<5183;l2++){
vtec_max_h1214_on[l2]=vtec_max_on[l2];
if(vtec_30d_quartile_h1214_on[l2][3]!=0.){vdtec_h1214_on[l2]=(vtec_max_h1214_on[l2]-vtec_30d_quartile_h1214_on[l2][1])/vtec_30d_quartile_h1214_on[l2][3];}
}
}

}
//end loop+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

for(int l3=0;l3<5183;l3++){
if(vtec_max_h24_on[l3]==0. || vtec_max_h46_on[l3]==0. || vtec_max_h1214_on[l3]==0.){vtec_max_h24_on[l3]=1.0E-5;vtec_max_h46_on[l3]=1.0E-5;vtec_max_h1214_on[l3]=1.0E-5;}
for(int l4=0;l4<4;l4++){
if(vtec_30d_quartile_h24_on[l3][l4]==0. || vtec_30d_quartile_h46_on[l3][l4]==0. || vtec_30d_quartile_h1214_on[l3][l4]==0.){ 
 vtec_30d_quartile_h24_on[l3][l4]=1.0E-5;vtec_30d_quartile_h46_on[l3][l4]=1.0E-5;vtec_30d_quartile_h1214_on[l3][l4]=1.0E-5;
}
} 
   // if(vdtec_h46_on[l3]>=-2.||vdtec_h46_on[l3]<=-5.){vdtec_h46_on[l3]=-7;}
    //if(vdtec_h24_on[l3]<=2.){vdtec_h24_on[l3]=-7.;}
  //  if(vdtec_h46_on[l3]<=2.){vdtec_h46_on[l3]=-7.;}
//    if(vdtec_h1214_on[l3]>=-2.){vdtec_h1214_on[l3]=-7.;} //for changing the background color of the plot

  //  if(vdtec_h1214_on[l3]<=2.){vdtec_h1214_on[l3]=-7.;} 
}

TH2D *h00 = new TH2D("h00"," ",73,-182.5,182.5,71,-88.75,88.75);
TH2D *h01 = new TH2D("h01"," ",73,-182.5,182.5,71,-88.75,88.75);
TH2D *h02 = new TH2D("h02"," ",73,-182.5,182.5,71,-88.75,88.75);
TH2D *h10 = new TH2D("h10"," ",73,-182.5,182.5,71,-88.75,88.75);
TH2D *h11 = new TH2D("h11"," ",73,-182.5,182.5,71,-88.75,88.75);
TH2D *h12 = new TH2D("h12"," ",73,-182.5,182.5,71,-88.75,88.75);
TH2D *h20 = new TH2D("h20"," ",73,-182.5,182.5,71,-88.75,88.75);
TH2D *h21 = new TH2D("h21"," ",73,-182.5,182.5,71,-88.75,88.75);
TH2D *h22 = new TH2D("h22"," ",73,-182.5,182.5,71,-88.75,88.75);

float loneq_on[1],lateq_on[1],loneq_size_on[1],lateq_size_on[1];
loneq_on[0]=128.17;lateq_on[0]=-0.52;
loneq_size_on[0]=138.35;lateq_size_on[0]=-0.52;
float lon_zoom_on[5]={95.,170.,170.,95.,95.},lat_zoom_on[5]={-22.,-22.,18.,18.,-22.};


float mag_eq_lon[73];
float mag_eq_lat[73];

fstream in3;
in3.open("/home/wang/Projects/EQ/otherdata/Mag_Equator.txt");;
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
h00->Fill(teclon[k1],teclat[k2],vtec_max_h24_on[k1+73*k2]);
h10->Fill(teclon[k1],teclat[k2],vtec_max_h46_on[k1+73*k2]);
h20->Fill(teclon[k1],teclat[k2],vtec_max_h1214_on[k1+73*k2]);

h01->Fill(teclon[k1],teclat[k2],vtec_30d_quartile_h24_on[k1+73*k2][1]);
h11->Fill(teclon[k1],teclat[k2],vtec_30d_quartile_h46_on[k1+73*k2][1]);
h21->Fill(teclon[k1],teclat[k2],vtec_30d_quartile_h1214_on[k1+73*k2][1]);

if(vdtec_h24_on[k1+73*k2]>2. || vdtec_h24_on[k1+73*k2]<-2.){
h02->Fill(teclon[k1],teclat[k2],vdtec_h24_on[k1+73*k2]);
//if(teclon[k1]>95. && teclon[k1]<170. && teclat[k2]>-22. && teclat[k2]<18.){cout<<"1st "<<(vdtec_h24_on[k1+73*k2] - 2) /2 * 100<<"%"<<endl;}
}
if(vdtec_h46_on[k1+73*k2]>2. || vdtec_h46_on[k1+73*k2]<-2.){
h12->Fill(teclon[k1],teclat[k2],vdtec_h46_on[k1+73*k2]);
//if(teclon[k1]>95. && teclon[k1]<170. && teclat[k2]>-22. && teclat[k2]<18.){cout<<"2nd "<<(vdtec_h46_on[k1+73*k2] - 2) /2 * 100<<"%"<<endl;}
}
if(vdtec_h1214_on[k1+73*k2]>2. || vdtec_h1214_on[k1+73*k2]<-2){
h22->Fill(teclon[k1],teclat[k2],vdtec_h1214_on[k1+73*k2]);
//if(teclon[k1]>95. && teclon[k1]<170. && teclat[k2]>-22. && teclat[k2]<18.){cout<<"3rd "<<(vdtec_h1214_on[k1+73*k2] - 2) /2 * 100<<"%"<<endl;}
//if(teclon[k1]>95. && teclon[k1]<170. && teclat[k2]>-22. && teclat[k2]<18.){cout<<"3rd "<<(vdtec_h1214_on[k1+73*k2] +  2) /2 * 100<<"%"<<endl;}
}
}}   

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

pad01->cd();
pad01->SetTopMargin(0.08);
pad01->SetBottomMargin(0.11);
pad01->SetLeftMargin(0.10);
pad01->SetRightMargin(0.09);
h00->Draw("colz ");
h00->GetXaxis()->SetNdivisions(210);
h00->GetYaxis()->SetNdivisions(410);
h00->GetXaxis()->SetTickLength(0.025);
h00->GetYaxis()->SetTickLength(0.015);
h00->GetXaxis()->SetLabelSize(0.05);
h00->GetYaxis()->SetLabelSize(0.05);
h00->GetXaxis()->SetLabelOffset(0.02);
h00->GetXaxis()->CenterTitle();
h00->GetYaxis()->CenterTitle();
h00->GetYaxis()->SetDecimals(kTRUE);
h00->GetXaxis()->SetRangeUser(-180., 180.);
h00->GetYaxis()->SetRangeUser(-81.25, 86.25);

hm->Draw("CONT2 SAME");
hm->SetLineColor(kGray+3);
hm->SetLineWidth(0.1);

g02->SetLineColor(kRed);
g02->Draw("SAME");

g00->SetMarkerStyle(29);
g00->SetMarkerSize(1.2);
g00->SetMarkerColor(kRed);
g00->Draw("P SAME");

g_mag_eq->SetLineColor(kBlack);
g_mag_eq->Draw("C SAME");

gPad->Update();
TPaletteAxis *palette00 = (TPaletteAxis*)h00->GetListOfFunctions()->FindObject("palette");
if(palette00) {
   palette00->SetX1NDC(0.922); 
   palette00->SetX2NDC(0.935); 
   gPad->Modified(); 
}

h00->GetZaxis()->SetRangeUser(0,30);
h00->GetZaxis()->SetNdivisions(206);
h00->GetZaxis()->SetLabelSize(0.049);
h00->GetZaxis()->SetLabelOffset(0.008);
h00->GetZaxis()->SetTickLength(0.009);
h00->GetZaxis()->CenterTitle();

TLatex *t000 = new TLatex();
   t000->SetTextSize(0.064);
   t000->DrawLatex(-172,90,"2019/07/06 06:00-08:00 UT (14:00-16:00 LT) LLT map");

TLatex *t001 = new TLatex();
   t001->SetTextSize(0.07);
   t001->SetTextAngle(90);
   t001->DrawLatex(-213.,-35.,"Latitude (#circ)");

TLatex *t002 = new TLatex();
   t002->SetTextSize(0.05);
   t002->SetTextAngle(90);
   t002->DrawLatex(220.,-30.,"TEC (TECU)");

TLatex *t003 = new TLatex();
   t003->SetTextSize(0.07);
   t003->DrawLatex(-225.,90.,"(a)");

pad02->cd();
pad02->SetTopMargin(0.08);
pad02->SetBottomMargin(0.11);
pad02->SetLeftMargin(0.10);
pad02->SetRightMargin(0.09);
h01->Draw("colz ");
h01->GetXaxis()->SetNdivisions(210);
h01->GetYaxis()->SetNdivisions(410);
h01->GetXaxis()->SetTickLength(0.025);
h01->GetYaxis()->SetTickLength(0.015);
h01->GetXaxis()->SetLabelSize(0.05);
h01->GetYaxis()->SetLabelSize(0.05);
h01->GetXaxis()->SetLabelOffset(0.02);
h01->GetXaxis()->CenterTitle();
h01->GetYaxis()->CenterTitle();
h01->GetYaxis()->SetDecimals(kTRUE);
h01->GetXaxis()->SetRangeUser(-180., 180.);
h01->GetYaxis()->SetRangeUser(-81.25, 86.25);

hm->Draw("CONT2 SAME");
hm->SetLineColor(kGray+3);
hm->SetLineWidth(0.1);

g02->SetLineColor(kRed);
g02->Draw("SAME");

g00->SetMarkerStyle(29);
g00->SetMarkerSize(1.2);
g00->SetMarkerColor(kRed);
g00->Draw("P SAME");

g_mag_eq->SetLineColor(kBlack);
g_mag_eq->Draw("C SAME");

gPad->Update();
TPaletteAxis *palette01 = (TPaletteAxis*)h01->GetListOfFunctions()->FindObject("palette");
if(palette01) {
   palette01->SetX1NDC(0.922); 
   palette01->SetX2NDC(0.935); 
   gPad->Modified(); 
}

h01->GetZaxis()->SetRangeUser(0,30);
h01->GetZaxis()->SetNdivisions(206);
h01->GetZaxis()->SetLabelSize(0.049);
h01->GetZaxis()->SetLabelOffset(0.008);
h01->GetZaxis()->SetTickLength(0.009);
h01->GetZaxis()->CenterTitle();

TLatex *t010 = new TLatex();
   t010->SetTextSize(0.064);
   t010->DrawLatex(-190,90,"2019/06/07-07/07 06:00-08:00 UT 30-day median LLT map");

TLatex *t011 = new TLatex();
   t011->SetTextSize(0.07);
   t011->SetTextAngle(90);
   t011->DrawLatex(-213.,-35.,"Latitude (#circ)");

TLatex *t012 = new TLatex();
   t012->SetTextSize(0.05);
   t012->SetTextAngle(90);
   t012->DrawLatex(220.,-30.,"TEC (TECU)");

TLatex *t013 = new TLatex();
   t013->SetTextSize(0.07);
   t013->DrawLatex(-225.,90.,"(b)");

pad03->cd();
pad03->SetTopMargin(0.08);
pad03->SetBottomMargin(0.14);
pad03->SetLeftMargin(0.10);
pad03->SetRightMargin(0.09);
h02->Draw("colz ");
h02->GetXaxis()->SetNdivisions(210);
h02->GetYaxis()->SetNdivisions(410);
h02->GetXaxis()->SetTickLength(0.025);
h02->GetYaxis()->SetTickLength(0.015);
h02->GetXaxis()->SetLabelSize(0.05);
h02->GetYaxis()->SetLabelSize(0.05);
h02->GetXaxis()->SetLabelOffset(0.02);
h02->GetXaxis()->CenterTitle();
h02->GetYaxis()->CenterTitle();
h02->GetYaxis()->SetDecimals(kTRUE);
h02->GetXaxis()->SetRangeUser(95., 170.);
h02->GetYaxis()->SetRangeUser(-22., 18.);
//h02->GetXaxis()->SetRangeUser(-180., 180.);
//h02->GetYaxis()->SetRangeUser(-81.25, 86.25);


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
TPaletteAxis *palette02 = (TPaletteAxis*)h02->GetListOfFunctions()->FindObject("palette");
if(palette02) {
   palette02->SetX1NDC(0.922); 
   palette02->SetX2NDC(0.935); 
   gPad->Modified(); 
}

h02->GetZaxis()->SetRangeUser(0.,3.);
//h02->GetZaxis()->SetRangeUser(-3.,7.);
//h02->GetZaxis()->SetRangeUser(-4.,8.);
h02->GetZaxis()->SetNdivisions(204);
h02->GetZaxis()->SetLabelSize(0.049);
h02->GetZaxis()->SetLabelOffset(0.008);
h02->GetZaxis()->SetTickLength(0.009);
h02->GetZaxis()->CenterTitle();

TLatex *t020 = new TLatex();
   t020->SetTextSize(0.064);
   t020->DrawLatex(90.5,20.,"2019/07/06 06:00-08:00 UT difference LLT (|#DeltaTEC|>2.0) map");
   //t020->DrawLatex(-150,90.,"2019/07/06 06:00-08:00UT difference LLT map");
   
TLatex *t021 = new TLatex();
   t021->SetTextSize(0.07);
   t021->SetTextAngle(90);
   t021->DrawLatex(86.5,-12,"Latitude (#circ)");
  // t021->DrawLatex(-213,-35,"Latitude (#circ)");

TLatex *t022 = new TLatex();
   t022->SetTextSize(0.066);
   t022->DrawLatex(122.5,-30.5,"Longitude (#circ)");
   //t022->DrawLatex(-50,-105,"Longitude (#circ)");

TLatex *t023 = new TLatex();
   t023->SetTextSize(0.05);
   t023->SetTextAngle(90);
   t023->DrawLatex(180.,-5.5,"#DeltaTEC");
   //t023->DrawLatex(220.,-30,"TEC(TECU)");

TLatex *t024 = new TLatex();
   t024->SetTextSize(0.075);
   t024->DrawLatex(83.,20.,"(d)");
  // t024->DrawLatex(-225,90.,"(c)");

//++++++++++++++++
pad04->cd();
pad04->SetTopMargin(0.08);
pad04->SetBottomMargin(0.11);
pad04->SetLeftMargin(0.08);
pad04->SetRightMargin(0.09);
h10->Draw("colz ");
h10->GetXaxis()->SetNdivisions(210);
h10->GetYaxis()->SetNdivisions(410);
h10->GetXaxis()->SetTickLength(0.025);
h10->GetYaxis()->SetTickLength(0.015);
h10->GetXaxis()->SetLabelSize(0.05);
h10->GetYaxis()->SetLabelSize(0.05);
h10->GetXaxis()->SetLabelOffset(0.02);
h10->GetXaxis()->CenterTitle();
h10->GetYaxis()->CenterTitle();
h10->GetYaxis()->SetDecimals(kTRUE);
h10->GetXaxis()->SetRangeUser(-180., 180.);
h10->GetYaxis()->SetRangeUser(-81.25, 86.25);

hm->Draw("CONT2 SAME");
hm->SetLineColor(kGray+3);
hm->SetLineWidth(0.1);

g02->SetLineColor(kRed);
g02->Draw("SAME");

g00->SetMarkerStyle(29);
g00->SetMarkerSize(1.2);
g00->SetMarkerColor(kRed);
g00->Draw("P SAME");

g_mag_eq->SetLineColor(kBlack);
g_mag_eq->Draw("C SAME");

gPad->Update();
TPaletteAxis *palette10 = (TPaletteAxis*)h10->GetListOfFunctions()->FindObject("palette");
if(palette10) {
   palette10->SetX1NDC(0.922); 
   palette10->SetX2NDC(0.935); 
   gPad->Modified(); 
}

h10->GetZaxis()->SetRangeUser(0,30);
h10->GetZaxis()->SetNdivisions(206);
h10->GetZaxis()->SetLabelSize(0.049);
h10->GetZaxis()->SetLabelOffset(0.008);
h10->GetZaxis()->SetTickLength(0.009);
h10->GetZaxis()->CenterTitle();

TLatex *t100 = new TLatex();
   t100->SetTextSize(0.064);
   t100->DrawLatex(-173,90,"2019/07/11 08:00-10:00 UT (16:00-18:00 LT) LLT map");

TLatex *t101 = new TLatex();
   t101->SetTextSize(0.05);
   t101->SetTextAngle(90);
   t101->DrawLatex(220.,-30.,"TEC (TECU)");

pad05->cd();
pad05->SetTopMargin(0.08);
pad05->SetBottomMargin(0.11);
pad05->SetLeftMargin(0.08);
pad05->SetRightMargin(0.09);
h11->Draw("colz ");
h11->GetXaxis()->SetNdivisions(210);
h11->GetYaxis()->SetNdivisions(410);
h11->GetXaxis()->SetTickLength(0.025);
h11->GetYaxis()->SetTickLength(0.015);
h11->GetXaxis()->SetLabelSize(0.05);
h11->GetYaxis()->SetLabelSize(0.05);
h11->GetXaxis()->SetLabelOffset(0.02);
h11->GetXaxis()->CenterTitle();
h11->GetYaxis()->CenterTitle();
h11->GetYaxis()->SetDecimals(kTRUE);
h11->GetXaxis()->SetRangeUser(-180., 180.);
h11->GetYaxis()->SetRangeUser(-81.25, 86.25);

hm->Draw("CONT2 SAME");
hm->SetLineColor(kGray+3);
hm->SetLineWidth(0.1);

g02->SetLineColor(kRed);
g02->Draw("SAME");

g00->SetMarkerStyle(29);
g00->SetMarkerSize(1.2);
g00->SetMarkerColor(kRed);
g00->Draw("P SAME");

g_mag_eq->SetLineColor(kBlack);
g_mag_eq->Draw("C SAME");

gPad->Update();
TPaletteAxis *palette11 = (TPaletteAxis*)h11->GetListOfFunctions()->FindObject("palette");
if(palette11) {
   palette11->SetX1NDC(0.922); 
   palette11->SetX2NDC(0.935); 
   gPad->Modified(); 
}

h11->GetZaxis()->SetRangeUser(0,30);
h11->GetZaxis()->SetNdivisions(206);
h11->GetZaxis()->SetLabelSize(0.049);
h11->GetZaxis()->SetLabelOffset(0.008);
h11->GetZaxis()->SetTickLength(0.009);
h11->GetZaxis()->CenterTitle();

TLatex *t110 = new TLatex();
   t110->SetTextSize(0.064);
   t110->DrawLatex(-191,90,"2019/06/12-07/12 08:00-10:00 UT 30-day median LLT map");

TLatex *t111 = new TLatex();
   t111->SetTextSize(0.05);
   t111->SetTextAngle(90);
   t111->DrawLatex(220.,-30.,"TEC (TECU)");

pad06->cd();
pad06->SetTopMargin(0.08);
pad06->SetBottomMargin(0.14);
pad06->SetLeftMargin(0.08);
pad06->SetRightMargin(0.09);
h12->Draw("colz ");
h12->GetXaxis()->SetNdivisions(210);
h12->GetYaxis()->SetNdivisions(410);
h12->GetXaxis()->SetTickLength(0.025);
h12->GetYaxis()->SetTickLength(0.015);
h12->GetXaxis()->SetLabelSize(0.05);
h12->GetYaxis()->SetLabelSize(0.05);
h12->GetXaxis()->SetLabelOffset(0.02);
h12->GetXaxis()->CenterTitle();
h12->GetYaxis()->CenterTitle();
h12->GetYaxis()->SetDecimals(kTRUE);
h12->GetXaxis()->SetRangeUser(95., 170.);
h12->GetYaxis()->SetRangeUser(-22., 18.);
//h12->GetXaxis()->SetRangeUser(-180., 180.);
//h12->GetYaxis()->SetRangeUser(-81.25, 86.25);

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
TPaletteAxis *palette12 = (TPaletteAxis*)h12->GetListOfFunctions()->FindObject("palette");
if(palette12) {
   palette12->SetX1NDC(0.922); 
   palette12->SetX2NDC(0.935); 
   gPad->Modified(); 
}

h12->GetZaxis()->SetRangeUser(0.,3.);
//h12->GetZaxis()->SetRangeUser(-3.,7.);
//h12->GetZaxis()->SetRangeUser(-4.,8.);
h12->GetZaxis()->SetNdivisions(204);
h12->GetZaxis()->SetLabelSize(0.049);
h12->GetZaxis()->SetLabelOffset(0.008);
h12->GetZaxis()->SetTickLength(0.009);
h12->GetZaxis()->CenterTitle();

TLatex *t120 = new TLatex();
   t120->SetTextSize(0.064);
   t120->DrawLatex(90.5,20.,"2019/07/11 08:00-10:00 UT difference LLT (|#DeltaTEC|>2.0) map");
  // t120->DrawLatex(-150,90.,"2019/07/11 08:00-10:00UT difference LLT map");

TLatex *t122 = new TLatex();
   t122->SetTextSize(0.066);
   t122->DrawLatex(122.5,-30.5,"Longitude (#circ)");
  //  t122->DrawLatex(-50,-105,"Longitude (#circ)");

TLatex *t123 = new TLatex();
   t123->SetTextSize(0.05);
   t123->SetTextAngle(90);
   t123->DrawLatex(180.,-5.5,"#DeltaTEC");
  // t123->DrawLatex(220,-30,"TEC(TECU)");

//++++++++++++++++
pad07->cd();
pad07->SetTopMargin(0.08);
pad07->SetBottomMargin(0.11);
pad07->SetLeftMargin(0.08);
pad07->SetRightMargin(0.09);
h20->Draw("colz ");
h20->GetXaxis()->SetNdivisions(210);
h20->GetYaxis()->SetNdivisions(410);
h20->GetXaxis()->SetTickLength(0.025);
h20->GetYaxis()->SetTickLength(0.015);
h20->GetXaxis()->SetLabelSize(0.05);
h20->GetYaxis()->SetLabelSize(0.05);
h20->GetXaxis()->SetLabelOffset(0.02);
h20->GetXaxis()->CenterTitle();
h20->GetYaxis()->CenterTitle();
h20->GetYaxis()->SetDecimals(kTRUE);
h20->GetXaxis()->SetRangeUser(-180., 180.);
h20->GetYaxis()->SetRangeUser(-81.25, 86.25);

hm->Draw("CONT2 SAME");
hm->SetLineColor(kGray+3);
hm->SetLineWidth(0.1);

g02->SetLineColor(kRed);
g02->Draw("SAME");

g00->SetMarkerStyle(29);
g00->SetMarkerSize(1.2);
g00->SetMarkerColor(kRed);
g00->Draw("P SAME");

g_mag_eq->SetLineColor(kBlack);
g_mag_eq->Draw("C SAME");

gPad->Update();
TPaletteAxis *palette20 = (TPaletteAxis*)h20->GetListOfFunctions()->FindObject("palette");
if(palette20) {
   palette20->SetX1NDC(0.922); 
   palette20->SetX2NDC(0.935); 
   gPad->Modified(); 
}

h20->GetZaxis()->SetRangeUser(0,30);
h20->GetZaxis()->SetNdivisions(206);
h20->GetZaxis()->SetLabelSize(0.049);
h20->GetZaxis()->SetLabelOffset(0.008);
h20->GetZaxis()->SetTickLength(0.009);
h20->GetZaxis()->CenterTitle();

TLatex *t200 = new TLatex();
   t200->SetTextSize(0.064);
   t200->DrawLatex(-173,90,"2019/07/13 06:00-08:00 UT (14:00-16:00 LT) LLT map");
 //  t200->DrawLatex(-173,90,"2019/07/20 16:00-18:00UT (00:00-02:00LT) LLT map");
   
TLatex *t201 = new TLatex();
   t201->SetTextSize(0.05);
   t201->SetTextAngle(90);
   t201->DrawLatex(220.,-30.,"TEC (TECU)");

pad08->cd();
pad08->SetTopMargin(0.08);
pad08->SetBottomMargin(0.11);
pad08->SetLeftMargin(0.08);
pad08->SetRightMargin(0.09);
h21->Draw("colz ");
h21->GetXaxis()->SetNdivisions(210);
h21->GetYaxis()->SetNdivisions(410);
h21->GetXaxis()->SetTickLength(0.025);
h21->GetYaxis()->SetTickLength(0.015);
h21->GetXaxis()->SetLabelSize(0.05);
h21->GetYaxis()->SetLabelSize(0.05);
h21->GetXaxis()->SetLabelOffset(0.02);
h21->GetXaxis()->CenterTitle();
h21->GetYaxis()->CenterTitle();
h21->GetYaxis()->SetDecimals(kTRUE);
h21->GetXaxis()->SetRangeUser(-180., 180.);
h21->GetYaxis()->SetRangeUser(-81.25, 86.25);

hm->Draw("CONT2 SAME");
hm->SetLineColor(kGray+3);
hm->SetLineWidth(0.1);

g02->SetLineColor(kRed);
g02->Draw("SAME");

g00->SetMarkerStyle(29);
g00->SetMarkerSize(1.2);
g00->SetMarkerColor(kRed);
g00->Draw("P SAME");

g_mag_eq->SetLineColor(kBlack);
g_mag_eq->Draw("C SAME");

gPad->Update();
TPaletteAxis *palette21 = (TPaletteAxis*)h21->GetListOfFunctions()->FindObject("palette");
if(palette21) {
   palette21->SetX1NDC(0.922); 
   palette21->SetX2NDC(0.935); 
   gPad->Modified(); 
}

h21->GetZaxis()->SetRangeUser(0,30);
h21->GetZaxis()->SetNdivisions(206);
h21->GetZaxis()->SetLabelSize(0.049);
h21->GetZaxis()->SetLabelOffset(0.008);
h21->GetZaxis()->SetTickLength(0.009);
h21->GetZaxis()->CenterTitle();

TLatex *t210 = new TLatex();
   t210->SetTextSize(0.064);
   t210->DrawLatex(-191,90,"2019/06/14-07/14 06:00-08:00 UT 30-day median LLT map");
 //  t210->DrawLatex(-191,90,"2019/06/21-07/21 16:00-18:00UT 30-day median LLT map");

TLatex *t211 = new TLatex();
   t211->SetTextSize(0.05);
   t211->SetTextAngle(90);
   t211->DrawLatex(220.,-30.,"TEC (TECU)");

pad09->cd();
pad09->SetTopMargin(0.08);
pad09->SetBottomMargin(0.14);
pad09->SetLeftMargin(0.08);
pad09->SetRightMargin(0.09);
h22->Draw("colz ");
h22->GetXaxis()->SetNdivisions(210);
h22->GetYaxis()->SetNdivisions(410);
h22->GetXaxis()->SetTickLength(0.025);
h22->GetYaxis()->SetTickLength(0.015);
h22->GetXaxis()->SetLabelSize(0.05);
h22->GetYaxis()->SetLabelSize(0.05);
h22->GetXaxis()->SetLabelOffset(0.02);
h22->GetXaxis()->CenterTitle();
h22->GetYaxis()->CenterTitle();
h22->GetYaxis()->SetDecimals(kTRUE);
h22->GetXaxis()->SetRangeUser(95., 170.);
h22->GetYaxis()->SetRangeUser(-22., 18.);
//h22->GetXaxis()->SetRangeUser(-180., 180.);
//h22->GetYaxis()->SetRangeUser(-81.25, 86.25);


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
TPaletteAxis *palette22 = (TPaletteAxis*)h22->GetListOfFunctions()->FindObject("palette");
if(palette22) {
   palette22->SetX1NDC(0.922); 
   palette22->SetX2NDC(0.935); 
   gPad->Modified(); 
}

h22->GetZaxis()->SetRangeUser(0.,3.);
//h22->GetZaxis()->SetRangeUser(-3,-1);
//h22->GetZaxis()->SetRangeUser(-3,7);
//h22->GetZaxis()->SetRangeUser(-4,8);
h22->GetZaxis()->SetNdivisions(204);
h22->GetZaxis()->SetLabelSize(0.049);
h22->GetZaxis()->SetLabelOffset(0.008);
h22->GetZaxis()->SetTickLength(0.009);
h22->GetZaxis()->CenterTitle();
h22->GetZaxis()->SetDecimals(kTRUE);

TLatex *t220 = new TLatex();
   t220->SetTextSize(0.064);
   t220->DrawLatex(90.,20.,"2019/07/13 06:00-08:00 UT difference LLT (|#DeltaTEC|>2.0) map");
 //  t220->DrawLatex(90.5,20.,"2019/07/20 16:00-18:00UT difference LLT (|#delta_{TEC}|>2.0) map");
//   t220->DrawLatex(-150,90.,"2019/07/13 06:00-08:00UT difference LLT map");
//    t220->DrawLatex(-150,90.,"2019/07/20 16:00-18:00UT difference LLT map");



TLatex *t221 = new TLatex();
   t221->SetTextSize(0.066);
   t221->DrawLatex(122.5,-30.5,"Longitude (#circ)");
  // t221->DrawLatex(-50,-105,"Longitude (#circ)");
   

TLatex *t222 = new TLatex();
   t222->SetTextSize(0.047);
   t222->SetTextAngle(90);
   t222->DrawLatex(180.,-5.5,"#DeltaTEC");
   //t222->DrawLatex(220.,-20,"TEC(TECU)");

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
vtec_30d_quartile_h24_on[vind][0]=q1;
vtec_30d_quartile_h24_on[vind][1]=q2;
vtec_30d_quartile_h24_on[vind][2]=q3;
vtec_30d_quartile_h24_on[vind][3]=iqr;
}

if(ngridx==1){
vtec_30d_quartile_h46_on[vind][0]=q1;
vtec_30d_quartile_h46_on[vind][1]=q2;
vtec_30d_quartile_h46_on[vind][2]=q3;
vtec_30d_quartile_h46_on[vind][3]=iqr;
}

if(ngridx==2){
vtec_30d_quartile_h1214_on[vind][0]=q1;
vtec_30d_quartile_h1214_on[vind][1]=q2;
vtec_30d_quartile_h1214_on[vind][2]=q3;
vtec_30d_quartile_h1214_on[vind][3]=iqr;
}

return iqr;
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
