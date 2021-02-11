//===============================================================
//Demeter Case Code "dtec_solar_magnetic_anomaly.C" 
//TAO DAN 2016.12.29 updated
//Wang 2020.08.11 updated 

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

using namespace std;

//g++ dtec_solar_magnetic_anomaly.C -O2 -Wall -fPIC -pthread -m64 -I/home/wang/CERNroot/root/root/include/ -L/home/wang/CERNroot/root/root/lib/ -lGraf -lGraf3d -lCore -lCint -lTree -lRIO  -lHist -lGpad -o dtec_solar_magnetic_anomaly

#define PI 3.141592653589793238463

static float vquartile_on[56][12][4]; 
static float vquartile_off[56][12][4]; 

int main(int argc, char** argv) 
{
float calc_quartile(int sflag, int nh, int nd, int nq, float *qteccode); //sflag=1 -NWC_ON;sflag=0 -NWC_OFF;
float bilinear_interpolation(int lflag, float *datapoint4); //lflag=1 -NWC_ON;lflag=0 -NWC_OFF;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//read index F10.7,Dst,Ap,AE, Kp
TH2D *h00 = new TH2D("h00"," ",56,0,56,12,0,24);
TH2D *h01 = new TH2D("h01"," ",56,0,56,12,0,24);
TH2D *h04 = new TH2D("h04"," ",56,0,56,12,0,24);

//NWC ON; EQ7.7-(-9.284,107.419)-20.0km-south of Java-06/07/17-08:19:26UTC
//
//EQ7.1-(1.621,126.416)-138km-Bitung Indonesia -2019/11/14-16:17:4.UTC
ifstream in0;
in0.open("/home/wang/Projects/EQ/Yahama/OMNI2_H0_MRG1HR_Yahama_20190430_0724.txt");
//EQ:06-07-17 08:19:26UTC; +75:05-03; +60:05-18; +45:06-02; +30:06-17; -10:07-27 
//
//EQ:2020/11/14 16:17:40UTC +75:08-30; +60:09-15;+45:09-30;+30:10-15;-10:11-24
int nrecord0=2064;//(85 + 1)*24 = 2064
float indyear77[nrecord0],indsec77[nrecord0],f10777[nrecord0],Kp77[nrecord0],Dst77[nrecord0],AE77[nrecord0],Ap77[nrecord0];

int nrd0=0; 

 while(in0.good())
{
 in0>>indyear77[nrd0]>>indsec77[nrd0]>>f10777[nrd0]>>Kp77[nrd0]>>Dst77[nrd0]>>AE77[nrd0]>>Ap77[nrd0];

 //2006/05/03-123th
 // indsec77[nrd0]=indsec77[nrd0]-10283400;
 Kp77[nrd0]=Kp77[nrd0]/10.;
 nrd0++;
} 
 in0.close();
//cout<<indsec77[2063]<<";"<<Dst77[2063]<<";"<<Ap77[2063]<<";"<<f10777[2063]<<";"<<Kp77[2063]<<";"<<AE77[2063]<<endl;

float dst_on[1032],kp_on[1032],f107_on[1032];
for(int i=0;i<1032;i++){
   dst_on[i]=0.;kp_on[i]=0.;f107_on[i]=0.;
   dst_on[i]=0.5*(Dst77[2*i]+Dst77[2*i+1]);
   kp_on[i]=0.5*(Kp77[2*i]+Kp77[2*i+1]);
   f107_on[i]=0.5*(f10777[2*i]+f10777[2*i+1]);
}

float dst_2h_on[56][12],kp_2h_on[56][12],f107_2h_on[56][12]; //56-day Dst and Kp value and 12 time periods per day
for(int k1=0;k1<56;k1++){
for(int k2=0;k2<12;k2++){
dst_2h_on[k1][k2]=0.;kp_2h_on[k1][k2]=0.;f107_2h_on[k1][k2]=0.;

dst_2h_on[k1][k2]=dst_on[360+12*k1+k2];
kp_2h_on[k1][k2]=kp_on[360+12*k1+k2];
f107_2h_on[k1][k2]=f107_on[360+12*k1+k2];

//cout<<f107_2h_on[k1][k2]<<";";
//if(k1==42 || k1==43){cout<<dst_2h_on[k1][k2]<<";";}
//if(k1==42 || k1==43){cout<<kp_2h_on[k1][k2]<<";";}

h00->SetBinContent(k1+1,12-k2,dst_2h_on[k1][k2]);
if(kp_2h_on[k1][k2]!=0.){h01->SetBinContent(k1+1,12-k2,kp_2h_on[k1][k2]);}
else
  {h01->SetBinContent(k1+1,12-k2,1.0E-7);}  //fill color for those grid which value=0.
h04->SetBinContent(k1+1,12-k2,f107_2h_on[k1][k2]);
}}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//read GPS-TEC-2HR data 30d-before and 10d-after
//float loneqoff=107.419;
//float lateqoff=-5.859;
//float loneqon=107.419;
//float lateqon=-9.284;
//float deg2rad=PI/180.;
//
//lat = 1.621
//lon = 126.416
float teclat[71],teclon[73];
for(int i=0;i<71;i++){teclat[i]=87.5-2.5*i;}
for(int i=0;i<73;i++){teclon[i]=-180+5.0*i;}

//GPS-TEC: NWC ON; EQ7.7-(107.419,-9.284)-20.0km-south of Java-06/07/17-08:19:26UTC
//GPS-TEC EQ7.1-(126.416,1.621)-138km-Bitung Indonesia-2020/11/14-16:17:40
ifstream in2;
in2.open("/home/wang/Projects/EQ/Yahama/GPS_TEC2HR_IGS_Yahama_20190430_0724.txt");
//EQ:06-07-17 08:19:26UTC; +75:05-03; +60:05-18; +45:06-02; +30:06-17; -10:07-27  
//EQ:2020/11/14 16:17:40UTC +75:08-31; +60:09-15;+45:10-15;+30:10-30;-10:11-24
int nrecord2=1033;
float tecyear_on[nrecord2],tecsec_on[nrecord2],teccode_on[nrecord2][4];

int nrd2=0;
float vtectemp_on=0.;

while(in2.good())
{
 in2>>tecyear_on[nrd2]>>tecsec_on[nrd2];
 for(int k2=0;k2<5183;k2++){
    if(k2!=2619 && k2!=2692 && k2!=2620 && k2!=2693){in2>>vtectemp_on;} //or 4085 4086 4156 4157
    else
    {
    //four data points nearest the epicenter (107.419,-9.284) (126,1)
        if(k2==2619){in2>>teccode_on[nrd2][0];}  //
        if(k2==2692){in2>>teccode_on[nrd2][1];}  //
        if(k2==2620){in2>>teccode_on[nrd2][2];}  //
        if(k2==2693){in2>>teccode_on[nrd2][3];}  //
    }
 }
 
 nrd2++;
} 
in2.close();
//cout<<tecyear_on[1032]<<";"<<tecsec_on[1032]<<";"<<teccode_on[1032][3]<<";"<<teccode_on[1032][2]<<";"<<endl;

float teccode_c2h_on[nrecord2-1][4],teccode_grid_on[nrecord2-1];
for(int i=0;i<nrecord2-1;i++){
   float arr_bilinear_inter_temp_on[4];
   for(int j=0;j<4;j++){
      teccode_c2h_on[i][j]=0.5*(teccode_on[i][j]+teccode_on[i+1][j]); //mid-piont value of 2h time period
      arr_bilinear_inter_temp_on[j]=0.;arr_bilinear_inter_temp_on[j]=teccode_c2h_on[i][j];
   }
 
   teccode_grid_on[i]=bilinear_interpolation(1,arr_bilinear_inter_temp_on);
}

//compute the M,UQ, LQ, they are the 30-day running median, upper quartile and lower quartile
//From 45 days before and 10 days after the EQ:06-07-17
//EQ:06-07-17 08:19:26UTC; +75:05-03; +60:05-18; +45:06-02; +30:06-17; -10:07-27
//
//From 45 days before and 10 days after the EQ:2019/11/14
//EQ:2020/11/14 16:17:40UTC +75:08-30; +60:09-15;+45:10-15;+30:10-30;-10:11-24
float tec_daily_on[86][12]; //12 time periods in per day
for(int k1=0;k1<86;k1++){
for(int k2=0;k2<12;k2++){
tec_daily_on[k1][k2]=0.;

tec_daily_on[k1][k2]=teccode_grid_on[12*k1+k2];
}}

float tec_ob_on[56][12]; //56-day observed TEC value and 12 time periods per day
for(int k1=0;k1<56;k1++){
for(int k2=0;k2<12;k2++){
tec_ob_on[k1][k2]=0.;

tec_ob_on[k1][k2]=teccode_grid_on[360+12*k1+k2];
}}

int nq1_0=0;
float tec_quartile_temp_on[30]; 
for(int l1=0;l1<12;l1++){
for(int l2=0;l2<56;l2++){
   for(int l3=0;l3<30;l3++){
      tec_quartile_temp_on[l3]=0.;
      tec_quartile_temp_on[l3]=tec_daily_on[l2+l3][l1];
   }

   calc_quartile(1,l1,l2,30,tec_quartile_temp_on);  //call function for calculating the 30-day running quartile values
   //cout<<tec_ob_on[l2][l1]<<";"<<vquartile_on[l2][l1][1]<<";"<<vquartile_on[l2][l1][3]<<";"<<endl;
   if(vquartile_on[l2][l1][0]==0.){nq1_0++;}
}}
//cout<<nq1_0<<endl;

float deq_on[1],heq_on[1];
deq_on[0]=45.5;
heq_on[0]=24.-(17.0+10./60.+0./3600.);

TH2D *h02 = new TH2D("h02"," ",56,0,56,12,0,24);
TGraph *g02 = new TGraph(1,deq_on,heq_on);
TH2D *h03 = new TH2D("h03"," ",56,0,56,12,0,24);

//int ng2=0;
float dtec_on[56][12],abs_dtec_on[56][12],teclb_on[56][12],tecub_on[56][12];
for(int k1=0;k1<56;k1++){
for(int k2=0;k2<12;k2++){
    dtec_on[k1][k2]=0.;abs_dtec_on[k1][k2]=0.;
    
    teclb_on[k1][k2]=vquartile_on[k1][k2][1]-2.0*(vquartile_on[k1][k2][1]-vquartile_on[k1][k2][0]);
    tecub_on[k1][k2]=vquartile_on[k1][k2][1]+2.0*(vquartile_on[k1][k2][2]-vquartile_on[k1][k2][1]);
    //if(tec_ob_on[k1][k2]>tecub_on[k1][k2] || tec_ob_on[k1][k2]<teclb_on[k1][k2]){ng2++;cout<<k1<<";"<<k2<<";"<<endl;}

    dtec_on[k1][k2]=(tec_ob_on[k1][k2]-vquartile_on[k1][k2][1])/vquartile_on[k1][k2][3];
    h02->SetBinContent(k1+1,12-k2,dtec_on[k1][k2]);
//	 if(dtec_on[k1][k2] < 0){dtec_on[k1][k2] = 0;}
   abs_dtec_on[k1][k2]=abs(dtec_on[k1][k2]);  
    //if(abs_dtec_on[k1][k2]>2.0){ng2++;cout<<k1<<";"<<k2<<";"<<dtec_on[k1][k2]<<endl;}
    if(abs_dtec_on[k1][k2]>3.5){abs_dtec_on[k1][k2]=0;}
    if(abs_dtec_on[k1][k2]>2. && dst_2h_on[k1][k2]>-30. && kp_2h_on[k1][k2]<3.0 && f107_2h_on[k1][k2]<100.0){
       h03->SetBinContent(k1+1,12-k2,dtec_on[k1][k2]);
   //   cout<<k1<<";"<<k2<<";"<<dtec_on[k1][k2]<<";"<<endl;
    }
    else
    {h03->SetBinContent(k1+1,12-k2,1.0E-7);}  //fill color for those grid which value=0.
}}
//cout<<ng2<<endl;
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++





//plotting================================================================
gStyle->SetOptStat(000);
gStyle->SetPalette(55);
//gStyle->SetTitleFontSize(0.0);
gStyle->SetErrorX(0);

TApplication theApp("App", &argc, argv);

TCanvas *c0 = new TCanvas("c0","TEC_ANOMALY",800,500);
c0->Connect("Closed()", "TApplication", &theApp, "Terminate()");

TPad *pad01 = new TPad("pad01","pad01",0.01,0.50,0.49,0.99);
pad01->Draw();
TPad *pad02 = new TPad("pad02","pad02",0.5,0.50,0.99,0.99);
pad02->Draw();
TPad *pad03 = new TPad("pad03","pad03",0.01,0.01,0.49,0.5);
pad03->Draw();
TPad *pad04 = new TPad("pad04","pad04",0.5,0.01,0.99,0.5);
pad04->Draw();

pad01->cd();
pad01->SetTopMargin(0.15);
pad01->SetBottomMargin(0.1);
pad01->SetLeftMargin(0.12);
h00->Draw("colz ");
h00->GetZaxis()->SetRangeUser(-40,40);
h00->GetXaxis()->SetTickLength(0.);
h00->GetYaxis()->SetTickLength(0.);
h00->GetXaxis()->SetLabelSize(0.);
h00->GetYaxis()->SetLabelSize(0.);
h00->GetXaxis()->SetTitleOffset(1.0);
h00->GetXaxis()->SetTitleSize(0.055);
h00->GetXaxis()->CenterTitle();
h00->GetYaxis()->CenterTitle();
h00->GetYaxis()->SetDecimals(kTRUE);

TLatex *t000 = new TLatex();
   t000->SetTextSize(0.06);
   t000->DrawLatex(24,25,"Dst index");

TLatex *t001 = new TLatex();
   t001->SetTextSize(0.06);
   t001->SetTextAngle(90);
   t001->DrawLatex(-4.5,8,"Time (UT)");

TLatex *t002 = new TLatex();
   t002->SetTextSize(0.06);
   t002->DrawLatex(-6.5,25,"(a)");

TLatex *t003 = new TLatex();
   t003->SetTextSize(0.048);
   t003->DrawLatex(56.5,24.5,"(nT)");

TGaxis *axis000 = new TGaxis(-0.5,0.,55.5,0.,-46,10,1010,"+RS");
axis000->SetLabelOffset(0.007);
axis000->SetTickSize(0.035);
axis000->Draw();

TGaxis *axis001 = new TGaxis(0.0000001,24.,0.0,0.,0,24,12,"+RS");
axis001->SetLabelOffset(-0.0035);
axis001->SetTickSize(0.02);
axis001->Draw();

g02->SetMarkerStyle(29);
g02->SetMarkerSize(1.0);
g02->SetMarkerColor(kRed);
g02->Draw("P SAME");

TLatex *t004 = new TLatex();
   t004->SetTextSize(0.049);
   t004->SetTextColor(kRed);
   t004->DrawLatex(45.5,4.8,"EQ time");

gPad->Update();
TPaletteAxis *palette00 = (TPaletteAxis*)h00->GetListOfFunctions()->FindObject("palette");
if(palette00) {
   palette00->SetX1NDC(0.92); 
   palette00->SetX2NDC(0.935); 
   gPad->Modified(); 
}
	
h00->GetZaxis()->SetLabelSize(0.040);
h00->GetZaxis()->SetTickLength(0.01);
h00->GetZaxis()->CenterTitle();

pad02->cd();
pad02->SetTopMargin(0.15);
pad02->SetBottomMargin(0.1);
pad02->SetLeftMargin(0.12);
h01->Draw("colz ");
h01->GetXaxis()->SetTickLength(0.);
h01->GetYaxis()->SetTickLength(0.);
h01->GetXaxis()->SetLabelSize(0.);
h01->GetYaxis()->SetLabelSize(0.);
h01->GetXaxis()->SetTitleOffset(1.0);
h01->GetXaxis()->SetTitleSize(0.055);
h01->GetXaxis()->CenterTitle();
h01->GetYaxis()->CenterTitle();
h01->GetYaxis()->SetDecimals(kTRUE);

TLatex *t010 = new TLatex();
   t010->SetTextSize(0.06);
   t010->DrawLatex(25,25,"Kp index");

TLatex *t011 = new TLatex();
   t011->SetTextSize(0.06);
   t011->SetTextAngle(90);
   t011->DrawLatex(-4.5,8,"Time (UT)");

TLatex *t012 = new TLatex();
   t012->SetTextSize(0.06);
   t012->DrawLatex(-8.3,25,"(b)");

TGaxis *axis010 = new TGaxis(-0.5,0.,55.5,0.,-46,10,1010,"+RS");
axis010->SetLabelOffset(0.007);
axis010->SetTickSize(0.035);
axis010->Draw();

TGaxis *axis011 = new TGaxis(0.0000001,24.,0.0,0.,0,24,12,"+RS");
axis011->SetLabelOffset(-0.0035);
axis011->SetTickSize(0.02);
axis011->Draw();

g02->SetMarkerStyle(29);
g02->SetMarkerSize(1.0);
g02->SetMarkerColor(kRed);
g02->Draw("P SAME");

TLatex *t013 = new TLatex();
   t013->SetTextSize(0.049);
   t013->SetTextColor(kRed);
   t013->DrawLatex(45.5,4.8,"EQ time");

gPad->Update();
TPaletteAxis *palette01 = (TPaletteAxis*)h01->GetListOfFunctions()->FindObject("palette");
if(palette01) {
   palette01->SetX1NDC(0.92); 
   palette01->SetX2NDC(0.935); 
   gPad->Modified(); 
}
	
h01->GetZaxis()->SetLabelSize(0.040);
h01->GetZaxis()->SetTickLength(0.01);
h01->GetZaxis()->CenterTitle();

pad03->cd();
pad03->SetTopMargin(0.08);
pad03->SetBottomMargin(0.16);
pad03->SetLeftMargin(0.12);
h04->Draw("colz ");
h04->GetXaxis()->SetTickLength(0.);
h04->GetYaxis()->SetTickLength(0.);
h04->GetXaxis()->SetLabelSize(0.);
h04->GetYaxis()->SetLabelSize(0.);
h04->GetXaxis()->CenterTitle();
h04->GetYaxis()->CenterTitle();
h04->GetYaxis()->SetDecimals(kTRUE);

TLatex *t040 = new TLatex();
   t040->SetTextSize(0.06);
   t040->DrawLatex(23,25,"F10.7 index");

TLatex *t041 = new TLatex();
   t041->SetTextSize(0.06);
   t041->SetTextAngle(90);
   t041->DrawLatex(-4.5,8,"Time (UT)");

TLatex *t042 = new TLatex();
   t042->SetTextSize(0.057);
   t042->DrawLatex(2.5,-3.5,"Day relative to the earthquake (M7.2-Laiwui) day");

TLatex *t043 = new TLatex();
   t043->SetTextSize(0.06);
   t043->DrawLatex(-6.5,25,"(c)");

TLatex *t044 = new TLatex();
   t044->SetTextSize(0.048);
   t044->DrawLatex(56.1,25.0,"(sfu)");

TGaxis *axis040 = new TGaxis(-0.5,0.,55.5,0.,-46,10,1010,"+RS");
axis040->SetLabelOffset(0.007);
axis040->SetTickSize(0.035);
axis040->Draw();

TGaxis *axis041 = new TGaxis(0.0000001,24.,0.0,0.,0,24,12,"+RS");
axis041->SetLabelOffset(-0.0035);
axis041->SetTickSize(0.02);
axis041->Draw();

g02->SetMarkerStyle(29);
g02->SetMarkerSize(1.0);
g02->SetMarkerColor(kRed);
g02->Draw("P SAME");

TLatex *t025 = new TLatex();
   t025->SetTextSize(0.049);
   t025->SetTextColor(kRed);
   t025->DrawLatex(45.5,4.8,"EQ time");

gPad->Update();
TPaletteAxis *palette04 = (TPaletteAxis*)h04->GetListOfFunctions()->FindObject("palette");
if(palette04) {
   palette04->SetX1NDC(0.92); 
   palette04->SetX2NDC(0.935); 
   gPad->Modified(); 
}


h04->GetZaxis()->SetRangeUser(40,100);
h04->GetZaxis()->SetDecimals(kTRUE);	
h04->GetZaxis()->SetLabelSize(0.040);
h04->GetZaxis()->SetTickLength(0.01);
h04->GetZaxis()->CenterTitle();

pad04->cd();
pad04->SetTopMargin(0.08);
pad04->SetBottomMargin(0.16);
pad04->SetLeftMargin(0.12);
h03->Draw("colz ");
h03->GetXaxis()->SetTickLength(0.);
h03->GetYaxis()->SetTickLength(0.);
h03->GetXaxis()->SetLabelSize(0.);
h03->GetYaxis()->SetLabelSize(0.);
h03->GetXaxis()->CenterTitle();
h03->GetYaxis()->CenterTitle();
h03->GetYaxis()->SetDecimals(kTRUE);

TLatex *t030 = new TLatex();
   t030->SetTextSize(0.06);
   //t030->DrawLatex(16.,25,"TEC anomaly(|DTEC>2.0|)");
   t030->DrawLatex(16.,25,"TEC anomaly (|#DeltaTEC|>2.0)");

TLatex *t0300 = new TLatex();
   t0300->SetTextSize(0.045);
  // t0300->DrawLatex(26,25.08,"(Dst>-30nT,Kp<2.3 and |DTEC|>2.0)");

TLatex *t031 = new TLatex();
   t031->SetTextSize(0.06);
   t031->SetTextAngle(90);
   t031->DrawLatex(-5.7,7,"Time (#color[2]{LT}-UT)");

TLatex *t032 = new TLatex();
   t032->SetTextSize(0.057);
   t032->DrawLatex(3.5,-3.5,"Day relative to the earthquake (M7.2-Laiwui) day");

TLatex *t033 = new TLatex();
   t033->SetTextSize(0.06);
   t033->DrawLatex(-8.3,25,"(d)");

TLatex *t034 = new TLatex();
   t034->SetTextSize(0.044);
   t034->DrawLatex(55.,25.0,"(#DeltaTEC)");

TGaxis *axis030 = new TGaxis(-0.5,0.,55.5,0.,-46,10,1010,"+RS");
axis030->SetLabelOffset(0.007);
axis030->SetTickSize(0.035);
axis030->Draw();

TGaxis *axis031 = new TGaxis(0.0000001,24.,0.0,0.,0,24,12,"+RS");
axis031->SetLabelOffset(-0.0030);
axis031->SetTickSize(0.02);
axis031->Draw();

TGaxis *axis032 = new TGaxis(-3.0000001,23.999,-3.0,7,8,25,209,"+RS");
//TGaxis *axis032 = new TGaxis(-3.0000001,23.999,-3.0,8,9,25,208,"+RS");
axis032->SetLabelOffset(-0.0035);
axis032->SetTickSize(0.02);
axis032->SetLineColor(kRed);
axis032->SetLabelColor(kRed);
axis032->Draw();

TGaxis *axis033 = new TGaxis(-3.0000001,7.,-3.0,0.,1,8,204,"+RS");
//TGaxis *axis033 = new TGaxis(-3.0000001,8.,-3.0,0.,1,9,204,"+RS");
axis033->SetLabelOffset(-0.0035);
axis033->SetTickSize(0.06);
axis033->SetLineColor(kRed);
axis033->SetLabelColor(kRed);
axis033->Draw();

g02->SetMarkerStyle(29);
g02->SetMarkerSize(1.0);
g02->SetMarkerColor(kRed);
g02->Draw("P SAME");

TLatex *t035 = new TLatex();
   t035->SetTextSize(0.05);
   t035->SetTextColor(kRed);
   t035->DrawLatex(-5.,24.8,"#color[2]{LT}#color[4]{=}#color[1]{UT}#color[4]{+8h}");
  // t035->DrawLatex(-5.,24.8,"#color[2]{LT}#color[4]{=}#color[1]{UT}#color[4]{+9h}");

TLatex *t036 = new TLatex();
   t036->SetTextSize(0.049);
   t036->SetTextColor(kRed);
   t036->DrawLatex(45.5,4.8,"EQ time");

gPad->Update();
TPaletteAxis *palette03 = (TPaletteAxis*)h03->GetListOfFunctions()->FindObject("palette");
if(palette03){
   palette03->SetX1NDC(0.92); 
   palette03->SetX2NDC(0.935); 
   gPad->Modified(); 
}

h03->GetZaxis()->SetRangeUser(-3.,3.);
h03->GetZaxis()->SetDecimals(kTRUE);
h03->GetZaxis()->SetLabelSize(0.040);
h03->GetZaxis()->SetTickLength(0.01);
h03->GetZaxis()->CenterTitle();

theApp.Run();

}//End main()

//================================================================================
//================================================================================






//function for calculating the quartile Q1, Q2, Q3 and IQR
float calc_quartile(int sflag, int nh, int nd, int nq, float *qteccode)
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

float iqr=q3-q1;

if(sflag==1){
vquartile_on[nd][nh][0]=q1;
vquartile_on[nd][nh][1]=q2;
vquartile_on[nd][nh][2]=q3;
vquartile_on[nd][nh][3]=iqr;
}
else
{
vquartile_off[nd][nh][0]=q1;
vquartile_off[nd][nh][1]=q2;
vquartile_off[nd][nh][2]=q3;
vquartile_off[nd][nh][3]=iqr;
}
//cout<<nd<<";"<<nh<<"|"<<vquartile_on[nd][nh][0]<<";"<<vquartile_on[nd][nh][1]<<";"<<vquartile_on[nd][nh][2]<<";"<<vquartile_on[nd][nh][3]<<";"<<endl;
//if(q1==0.){
//   for(int i=0;i<30;i++){cout<<qteccode[i]<<";";}
//   cout<<endl;
//}

return iqr;
}

//function of bilinear interpolation
float bilinear_interpolation(int lflag, float *datapoint4)
{
//four data points nearest the epicenter (107.419,-9.284) and (107.419,-5.859)
//EQ_ON:Q12(105,-7.5)-datapoint4[0]; Q11(105,-10.)-datapoint4[1]; Q22(110,-7.5)-datapoint4[2]; Q21(110,-10.)-datapoint4[3];
//EQ_OFF:Q12(105,-5.0)-datapoint4[0]; Q11(105,-7.5)-datapoint4[1]; Q22(110,-5.0)-datapoint4[2]; Q21(110,-7.5)-datapoint4[3];
float x1,x2,x,y1,y2,y;
if(lflag==1){x1=125.;x2=130.;y1=0.;y2=-2.5;x=128.17;y=-0.52;}
else{x1=125.;x2=130.;y1=0.;y2=2.5;x=126.416;y=1.621;}

//R1(x,y1),R2(x,y2)-linear interpolation on x axis
float tecr1=datapoint4[1]*(x2-x)/(x2-x1)+datapoint4[3]*(x-x1)/(x2-x1);
float tecr2=datapoint4[0]*(x2-x)/(x2-x1)+datapoint4[2]*(x-x1)/(x2-x1);

//Peq(x,y)-linear interpolation on y axis
float tecpeq=tecr1*(y2-y)/(y2-y1)+tecr2*(y-y1)/(y2-y1);

return tecpeq;
}
