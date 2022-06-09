#include <iostream>
#include <algorithm>
#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TH3D.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TGraph.h>
#include <TVector3.h>
#include <TMultiGraph.h>
#include <TEllipse.h>
#include <TRandom.h>
#include <TLine.h>

#include "Garfield/AvalancheMC.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/SolidBox.hh"
#include "Garfield/Component.hh"
#include "Garfield/GeometrySimple.hh"
#include "Garfield/ComponentConstant.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/Plotting.hh"
#include "Garfield/SolidTube.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/ViewCell.hh"
#include "Garfield/ViewGeometry.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewSignal.hh"
#include "Garfield/Random.hh"
#include "TProfile.h"

// unit convention: 
// length: cm
// time: ns
// Voltage: V
// Charge: fC

// COMET
const double d_sense = 0.0025;  //diameter of wire
const double d_field = 0.0126;  //diameter of wire
const double PI = 3.141592653589793;

using namespace std;
using namespace Garfield;

class UserAnalysis{
  private:
    MediumMagboltz* fGas;
    GeometrySimple* fGeo;
    ComponentAnalyticField* fCompE_Field;
    Sensor* fSensor_e;
    double fSize_box;


    double ix;
    double emag;
    double vmag;
    double dfl;
    double dft;
    double alpha;
    double attach;


    TLine *lTrack;
    TCanvas* fMyCanvas;
    TH3D* hElectricFieldxy;
    TH2D* hElTrack;
    TH2D* hXT;
    TH2D* hTX;
    TH1D* hNElectrons;
    TH1D* hTrackLength;
    TH1D* hEdep;
    TH1D* hDCAandCIonPos;
    TGraph* gMinRIonisedPt;
    TH2D* DriftLine;
    TH1D* hGain;
    //output
    TFile* fFile_out;
    TTree* fTree_out;
    // Initial position of ionized electron
    double fInPosX;
    double fInPosY;
    double fInPosZ;
    double fInPosT;
    double fDriftTime;
    double fDriftTime_uniSmear;
    int fTDC;
    double fDCA;
    double fIonDist2Wire;
    double fED;
    std::vector<double> fDel_transE;
    std::vector<double> fDel_KE;
    std::vector<double> fDel_y;
    std::vector<double> fDel_x;
    double fTrackLength;
    double fDiffL;
    double fDiffT;
    double fDv;
    double fEMag;
    double fSignal;

    int fNElectron;
    int fNElectron_in_col;
    int fHV;
    int fNEvent;
    int fFileId;
    TString fOutDir;
    TString fSuffix;
    TString fPrimaryParticle;
    double fPrimaryMomentum;
    double fWidth;
    double fHeight;
  public:
    UserAnalysis(int HV,int nevents,int fileId,TString outDir, TString suffix, TString primary, double momentum):
      fGas(NULL),fGeo(NULL),fCompE_Field(NULL),fSensor_e(NULL),fMyCanvas(NULL),lTrack(NULL),
      hElectricFieldxy(NULL),fPrimaryParticle(primary), fPrimaryMomentum(momentum)
  {
    fSuffix = suffix;
    fFileId = fileId;
    fHV = HV;
    fNEvent = nevents;
    fOutDir = outDir;
    printf("%s \n",fOutDir.Data());
    fSize_box = 147.42;
    fMyCanvas = new TCanvas();
    /* hElTrack = new TH2D("e_track","e_track",1e4,-2.0,2.0,1e4,-2.0,2.0);
    hXT = new TH2D("XT","XT",100 ,-10.0, 10.0, 256, 0 , 512 );
    hXT->GetXaxis()->SetTitle("Drift distance [mm]");
    hXT->GetYaxis()->SetTitle("Drift time [ns]");
    hTX = new TH2D("TX","TX", 256, 0 , 512,100 ,-10.0, 10.0 );
    hTX->GetXaxis()->SetTitle("Drift time [ns]");
    hTX->GetYaxis()->SetTitle("Drift distance [mm]");

    hGain = new TH1D("hGain","hGain",200,0,50000);
    
    DriftLine = new TH2D("driftline","drift time - MC",1e4,-2.0 ,2.0,1e4,-2.0,2.0 );
    hNElectrons = new TH1D("hNElectrons","hNElectrons",200,0,200);
    hTrackLength = new TH1D("hTrackLength","hTrackLength",100,0,2.0);
    hEdep = new TH1D("hEdep","hEdep",100,0,10);
    hEdep -> GetXaxis()->SetTitle("energy loss [keV]/[cm]");
    gMinRIonisedPt = new TGraph(1);
    hDCAandCIonPos = new TH1D("hDCAandCIonPos","hDCAandCIonPos",100,-1,1);
    */
    std::cout << "check " << fSuffix << std::endl;
    fFile_out = new TFile(Form("%s/output_%d_%s.root",fOutDir.Data(),fileId,fSuffix.Data()),"RECREATE");
    /*fTree_out = new TTree("t","t");
    fTree_out->Branch("inPosX",&fInPosX);
    fTree_out->Branch("inPosY",&fInPosY);
    fTree_out->Branch("inPosZ",&fInPosZ);
    fTree_out->Branch("inPosT",&fInPosT);
    fTree_out->Branch("driftTime",&fDriftTime);
    fTree_out->Branch("driftTime_uniSmear",&fDriftTime_uniSmear);
    fTree_out->Branch("TDC",&fTDC);
    fTree_out->Branch("dca",&fDCA);
    fTree_out->Branch("ionDist2Wire",&fIonDist2Wire);
    fTree_out->Branch("ED",&fED);
    fTree_out->Branch("trackLength",&fTrackLength);
    // fTree_out->Branch("diff_L",&fDiffL);
    // fTree_out->Branch("diff_T",&fDiffT);
    // fTree_out->Branch("dv",&fDv);
    // fTree_out->Branch("emag",&fEMag);
    // fTree_out->Branch("delta_transE",&fDel_transE);
    // fTree_out->Branch("delta_ke",&fDel_KE);
    fTree_out->Branch("nElectron",&fNElectron);
    fTree_out->Branch("nElectron_in_col",&fNElectron_in_col); */
  }
    void SetupGas(TString gasFile, TString ionMobilityFile){
      fGas = new MediumMagboltz();
      const double pressure = 750.;
      const double temperature = 293.15;
      fGas->LoadGasFile(gasFile.Data());
      fGas->LoadIonMobility(ionMobilityFile.Data());
      fGas->SetMaxElectronEnergy(200.);
      fGas->Initialise(true);
      fGas->PrintGas();
    }

    void SetupGeo(TString geomFile){
      fGeo = new GeometrySimple();

      //Wire dimension  [cm]
      SolidBox* box = new SolidBox(0, 0, 0, 10, 10, fSize_box/2);
      fGeo->AddSolid(box, fGas);

      cout<<"\033[32m";
      cout<<fGeo->GetSolid(0,0,-73)<<"\n";
      cout<<fGeo->IsInside(0,0,-73)<<"\n";
      cout<<"\033[0m";

      fCompE_Field = new ComponentAnalyticField();
      fCompE_Field->SetGeometry(fGeo);

      const int senseHV = fHV;
      //const int senseHV = 500;
      //const int senseHV = 1800;
      const int fieldHV = 0;

      // load the geometry file
      TFile * ifile = new TFile(geomFile);
      TTree * itree = (TTree*) ifile->Get("t");
      double x,y;
      int type;
      itree->SetBranchAddress("x",&x);
      itree->SetBranchAddress("y",&y);
      itree->SetBranchAddress("type",&type);
      for (int i = 0; i<itree->GetEntries(); i++){
	itree->GetEntry(i);
	if (type==0){ // field wire
	  fCompE_Field->AddWire( x,  y, d_field, fieldHV, Form("f%d",i), 147.42);
	}
	else if (type==1){ // sense wire
	  fCompE_Field->AddWire( x,  y, d_sense, senseHV, Form("s%d",i), 147.42);
          fCompE_Field->AddReadout(Form("s%d",i));
	  
	}
      }


      fCompE_Field->SetMagneticField(0,0,0);

      //view geometry
      ViewGeometry *geoView;
      geoView->SetGeometry(fGeo);
      geoView->Plot3d();
    }

    void SetupSensor(){
      const Double_t tMin =  0.;     // [ns]
      const Double_t tMax =  1000.;     // [ns]
      const Double_t tStep = 0.02;   // 0.02 [ns]
      const Int_t nTimeBins = (Int_t) ((tMax - tMin)/tStep);

      fSensor_e = new Sensor();
      fSensor_e->AddComponent(fCompE_Field);
      fSensor_e->AddElectrode(fCompE_Field, "s24");
      fSensor_e->SetTimeWindow(0., tStep, nTimeBins);
      fSensor_e->SetArea(-0.8,-0.8,-73.71,0.8,0.8,73.71);
    }

    void MakeTrack(){

      TrackHeed* track = new TrackHeed();
      track->SetSensor(fSensor_e);
      track->SetParticle(fPrimaryParticle.Data());
      cout<<"Primary Momentum "<<fPrimaryMomentum<<" eV"<<endl;
      track->SetMomentum(fPrimaryMomentum);

      const int nEvents = fNEvent;

      track->EnableDebugging();

      AvalancheMicroscopic* drift = new AvalancheMicroscopic();
      drift->SetTimeWindow(0.,1000.);
      //drift->SetTimeSteps(0.02);
      drift->SetSensor(fSensor_e);
      drift->EnableSignalCalculation();

      //make tracks
      TH1D *hTimeAllClusters = new TH1D("hTimeAllClusters","hTimeAllClusters",1000,0,1000);
      TH2D *hXTAllClusters = new TH2D("hXTAllClusters","hXTAllClusters",100,-10,10,1000,0,1000);
      gRandom->SetSeed(fFileId);
      int from = fFileId * nEvents;
      int to = (fFileId+1) * nEvents;
      cout<<"From "<<from<<" to "<<to<<endl;
      //fFile_out->cd();
      for(int i = from; i< to; i++){

	if(i==1){ track->DisableDebugging(); }
	cout<<i<<"/"<<nEvents<<"\n";

	double alpha = PI/2.;   // Incident Angle FIXME

	double x0=gRandom->Uniform(-0.8,0.8);
	double y0=gRandom->Uniform(-0.8,0.8);
	double z0=gRandom->Uniform(-73.71,73.71);
	double t0=0.;
	double dx0=gRandom->Uniform(-1.,1.);
	double dy0=gRandom->Uniform(-1.,1.);
	double dz0=gRandom->Uniform(-1.,1.);

	track->NewTrack(x0,y0,z0,t0,dx0,dy0,dz0);

	double xc=0., yc=0., zc=0., tc=0.;
	int nc=0;
	double ec=0.;
	double extra=0.;
	double esum=0.;
	double nsum=0.;
	double nsum_in_col=0.;
	TVector3 minR_point;
	double minEnergy=1e5;
	double minTime=1e5;
	double minR=1e10;
	track->EnableDeltaElectronTransport();
	fDel_transE.clear();
	fDel_KE.clear();
	fDel_y.clear();
	fDel_x.clear();

	std::cout<<"Get cluster of ionized electrons"<<std::endl;
	while(track->GetCluster(xc,yc,zc,tc,nc,ec,extra)){
	  double r=sqrt(xc*xc+yc*yc+zc*zc);
	  double x,y,z,t;
	  int nofdlp;
	  if(fabs(yc)<=0.8 && fabs(xc)<=0.8){
	    int i=0;
	    double del_dx,del_dy,del_dz;
	    double te,ke;
	    double del_x,del_y,del_z;
	    track->GetElectron(i,del_x,del_y,del_z,te,ke,del_dx,del_dy,del_dz);

	    fDel_KE.push_back(ke);
	    fDel_transE.push_back(ec);
	    fDel_x.push_back(del_x);
	    fDel_y.push_back(del_y);
	    esum+=ec;
	    nsum++;
	    nsum_in_col+=nc;
	    //drift->DriftElectron(xc,yc,zc,tc);
	    drift->AvalancheElectron(del_x,del_y,del_z,te,ke,del_dx,del_dy,del_dz);
	    i++;
	  }

	}

	// get signals
	for(int j=0;j<100000;j++){

        fSignal += fSensor_e->GetSignal("s24",j);
	}

        std::cout<<"\033[32m";
	std::cout<<"number of ionized electrons"<<nsum_in_col<<std::endl;
	std::cout<<"Signal from the s24 wire"<<fSignal<<std::endl;
        std::cout<<"\033[0m";
	// view signals

       	/* gMinRIonisedPt->SetPoint(0,minR_point.X(),minR_point.Y());

	drift->DriftElectron(minR_point.X(), minR_point.Y(), minR_point.Z(),minTime);
	int nofdlp=drift->GetNumberOfDriftLinePoints();
	double dt = minTime;
	double x,y,z,t;
	drift->GetDriftLinePoint(nofdlp-1, x, y, z, t);
	double dt_final = t - dt;
	double xe2, ye2, ze2, te2;
	int status;
	drift->GetElectronEndpoint(0, x, y, z, t, xe2, ye2, ze2, te2, status);

	hXT->Fill(minR_point.X()*10,dt_final);
	hTX->Fill(dt_final,minR_point.X()*10);

	//printf("Cluster %f %f %f : time %f \n",xc,yc,zc,tc);
	double tl = 1.6;
	hNElectrons->Fill(nsum/tl);
	hEdep->Fill(esum * 1.e-3/tl);
	fInPosX = x;
	fInPosY = y;
	fInPosZ = z;
	fInPosT = t;
	fDriftTime = dt_final;
	fDriftTime_uniSmear = gRandom->Uniform(-4,4)+dt_final;
	//fTDC = std::floor(fDriftTime_uniSmear/0.96);
	fTDC = std::floor(fDriftTime_uniSmear/1.04);
	fDCA = sqrt(x0*x0+y0*y0+z0*z0);
	fIonDist2Wire = minR;
	fED = esum;
	fNElectron = nsum;
	fNElectron_in_col = nsum_in_col;
	fTrackLength = tl;
	fTree_out->Fill(); */

	fSensor_e->ClearSignal();
	fSignal = 0.;
      }

      /* hTimeAllClusters->Write();
      fTree_out->Write();
      fFile_out->Close();*/
	TCanvas *c_view = new TCanvas("view","", 1400, 600);
	c_view->cd();
	//ViewDrift *trackView;
	ViewDrift *driftView;
	//track->EnablePlotting(trackView);
	drift->EnablePlotting(driftView);
	driftView->SetArea(-0.8,-0.8,0.8,0.8);
	driftView->Plot(true);
    }


    TVector3 GetClosestIonizedElectron(TrackHeed* track){
      TVector3 point;
      double minR = 1e10;
      // Cluster coordinates
      double xc = 0., yc = 0., zc = 0., tc = 0.;
      // Number of electrons produced in a collision
      int nc = 0;
      // Energy loss in a collision
      double ec = 0.;
      // Dummy variable (not used at present)
      double extra = 0.;
      while (track->GetCluster(xc, yc, zc, tc, nc, ec, extra)) {
	double r = sqrt(xc*xc+yc*yc+zc*zc);
	if(r<minR){
	  point.SetX(xc) ;
	  point.SetY(yc) ;
	  point.SetZ(zc) ;
	}
      }
      return point;
    }

    void DrawFieldxy(){
      int binx=100; double minx=-3.0; double maxx=3.0; //FIXME
      int biny=100; double miny=-3.0; double maxy=3.0; //FIXME
      int binz=100; double minz=-74.0; double maxz=74.0; //FIXME
      hElectricFieldxy = new TH3D("hElectricFieldxy","hElectricFieldxy",
	  binx,minx,maxx,
	  biny,miny,maxy,
	  binz,minz,maxz
	  );
      double ex,ey,ez,v;
      int status;
      Medium * m = nullptr;
      for(int i=0;i<binx;i++){
	ix = minx+(maxx-minx)*i/(binx-1);
	for(int j=0;j<biny;j++){
	  double iy = miny+(maxy-miny)*j/(biny-1);
	  for(int k=0;k<binz;k++){
	    double iz = minz+(maxz-minz)*k/(binz-1);
	    fCompE_Field->ElectricField(ix,iy,iz,
		ex,ey,ez,v,
		m,
		status
		);
	    //printf("(%f %f %f) E (%f %f %f) V %f status %d \n",
	    //       ix,iy,iz,
	    //       ex,ey,ez,v,
	    //       status);
	    double emag = sqrt(ex*ex+ey*ey+ez*ez);
	    hElectricFieldxy->SetBinContent(i+1,j+1,k+1,v);
	  }
	}
      }
      fFile_out->cd();
      hElectricFieldxy->Write();
    }

#if 0
    void PrintWire(){
      const int nWires = fCompE_Field->GetNumberOfWires();

      double ex,ey,ez,v;
      int status;
      Medium * m;

      int binx=1000; double minx=-0.83745; double maxx=0.83745;
      for(int i=0;i<binx;i++){
	ix = minx+(maxx-minx)*i/(binx-1);
	fCompE_Field->ElectricField(ix,0,0,
	    ex,ey,ez,v,
	    m,
	    status
	    );

	emag = sqrt(ex*ex+ey*ey+ez*ez);
	if(m){
	  double vx,vy,vz;
	  m->ElectronVelocity(ex,ey,ez, 0,0,0, vx,vy,vz);
	  m->ElectronDiffusion(ex,ey,ez, 0,0,0, dfl, dft );
	  m->ElectronTownsend(ex,ey,ez, 0,0,0, alpha);
	  m->ElectronAttachment(ex,ey,ez, 0,0,0, attach);
	  vmag = sqrt(vx*vx+vy*vy+vz*vz);
	  fTree_out->Fill();
	}
      }

      fFile_out->cd();
      fTree_out->Write();
    }
#endif

    void Finalize(){
      fFile_out->Close();
    }
};

int main(int argc, char * argv[]) {

  TString gasFile = argv[1];
  TString ionMobFile = argv[2];
  TString geomFile = argv[3];
  int hv = atoi(argv[4]);
  int nevents = atoi(argv[5]);
  int fileId = atoi(argv[6]);
  TString outDir=argv[7];
  TString suffix=argv[8];
  TString primary_particle = argv[9];
  double primary_momentum = atof(argv[10]);
  if(argc<9){
    printf("/GainSaturation <gasfile><ionMobfile><geomFile><highvoltage><nevents><fileId><outDir><suffix><Prim_particle><Prim_momentum[eV]>");
    return 0;
  }
  for(int i=0;i<argc;i++){
    std::cout << i << " " << argv[i] << std::endl;
  }
  UserAnalysis *userCode = new UserAnalysis(
      hv,nevents,fileId,outDir,suffix,
      primary_particle,primary_momentum
      );
  userCode->SetupGas(gasFile,ionMobFile);
  userCode->SetupGeo(geomFile);
  userCode->SetupSensor();
  userCode->DrawFieldxy();
  //userCode->MakeTrack();
  userCode->Finalize();

  return 1;
}
