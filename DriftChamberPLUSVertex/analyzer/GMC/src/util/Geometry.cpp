#include "util/Geometry.h"
#include "TVector3.h"
#include <iostream>
#include "generated/GMCAnalyzer.h"
#include "generated/GMCGlobalSteering.h"

//change if you want to rotate the detector by ROTANGLE
#define ROTATE 0
#define ROTANGLE 0.21
//change if you want that the central tube has a different staggering
#define stagCentral 0
#define dstagCentral 2.5

using namespace std;

ClassImp(Geometry)

Geometry* Geometry::Instance = 0;

Geometry::Geometry()
{

  TVector3 c;
  Double_t R = gAnalyzer->GetGSP()->Gettube_radius();
  
  wire_posX[0] = 9999999;
  wire_posY[0] = 9999999;
  wire_posZ[0] = 9999999;
  
  for(int i=1; i<=16; i++){
    //even tubes
    if(i%2==0){
      if(i < 9){
	c.SetXYZ(0.,2.*R,5.*R-2*(i-2)*R); //up-tube
      }
      else {
	c.SetXYZ(0.,0.,5.*R-2*(i-10)*R); //down-tube
      }
    }
    else {
      if(i < 9){
	c.SetXYZ(0.,R,7.*R-2*(i-1)*R); //up-tube
      }
      else {
	c.SetXYZ(0.,-R,7.*R-2*(i-9)*R); //down-tube
      }
    }  
    if(stagCentral==1 && (i==5 || i==13)){
      c.SetXYZ(c.X(),c.Y()-dstagCentral,c.Z());
    }
    
    if(ROTATE==1){
      c.SetXYZ(c.X(),c.Z()*sin(ROTANGLE) + c.Y()*cos(ROTANGLE),c.Z()*cos(ROTANGLE) - c.Y()*sin(ROTANGLE));
    }

    wire_posX[i] = c.X();
    wire_posY[i] = c.Y();
    wire_posZ[i] = c.Z();
    
  }

  ifstream cabling("cabling.dat");

  int wire;

  cout << "CABLING" << endl;
  for(Int_t i=0;i<17;i++){
    tdc[wire] = -1;
  }

  for(Int_t i=0;i<17;i++){

    cabling >> wire ;
    cabling >> tdc[wire];

    if(wire != 0)
      cout << "Wire " << wire << " ---- TDC channel " << tdc[wire] << endl;
    else
      cout << "Trigger ---- TDC channel " << tdc[0] << endl;

  }

  cabling.close();

  ifstream offsets("time_offsets.dat");

  for(Int_t i=0;i<17;i++){

    offsets >> wire ;
    offsets >> offset[wire];

  }

}

Geometry::~Geometry()
{

}

Int_t Geometry::GetWireIndex(Int_t tdc_ch)
{

  for(Int_t i=0;i<17;i++)
    if(tdc_ch == tdc[i]) return i;

  return -1;

}

Geometry *Geometry::GetInstance()
{

  if(Instance) return Instance;

  Instance = new Geometry();
  return Instance; 

}

