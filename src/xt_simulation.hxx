#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TApplication.h"

#include "Garfield/MediumMagboltz.hh"
#include "Garfield/GeometryRoot.hh"
#include "Garfield/GeometrySimple.hh"
#include "Garfield/SolidBox.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/TrackPAI.hh"
#include "Garfield/DriftLineRKF.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/ViewGeometry.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewCell.hh"
#include "Garfield/ViewDrift.hh"

using namespace Garfield;

const double d_sense = 0.0025;  //diameter of wire
const double d_field = 0.0126;  //diameter of wire

class simXT
{
public:
        simXT(int, int, int, std::string, double);
	~simXT();

	void SetupGas(std::string, std::string);
	void SetupGeo(std::string);
	void SetupSensor();
	void MakeTrack(int nTrack = 1);
	void DrawField();

private:
        MediumMagboltz* fGas;
	GeometrySimple* fGeo;
	ComponentAnalyticField* fCompE_field;
	Sensor* fSensor_e;
	double fSize_box;

	TCanvas* fMyCanvas;
	std::string fPrimaryParticle;
	int fPrimaryMomentum;
	int fFileId;
	int fHV;
	int fNEvent;
	TFile* fFile_out;
	TTree* fTree_out;

	double wire_x = -1e5;
	double wire_y = -1e5;
	int wire_type = -1e5;
		
	std::vector<double> f_x0;
	std::vector<double> f_y0;
	std::vector<double> f_z0;
	std::vector<double> f_t0;	
	std::vector<double> f_x1;
	std::vector<double> f_y1;
	std::vector<double> f_z1;
	std::vector<double> f_t1;	
	std::vector<double> f_ec_dis;
	std::vector<double> f_e_cls;
	std::vector<int> f_n_ele;
	int f_n_cls;
	int f_Nt;                             //Number of electrons generated in a track
	int f_index_e;
	double f_edep;
        double f_DriftDistance;
        double f_DriftTime;
	double f_t_wire;
	double f_DOCA;
	int f_Status;

};
