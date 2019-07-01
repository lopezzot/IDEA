#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <TVector3.h> 
#include <TObject.h> 
#include <TGraph.h> 

class Geometry : public TObject {

 protected:

  Double_t wire_posX[17];
  Double_t wire_posY[17];
  Double_t wire_posZ[17];
  int tdc[17];
  Double_t offset[17];

  Geometry();
  virtual ~Geometry();

 private:

  static Geometry *Instance;

 public:

  static Geometry *GetInstance();
  Double_t Getwire_posX(Int_t wire) { return wire_posX[wire]; }
  Double_t Getwire_posY(Int_t wire) { return wire_posY[wire]; }
  Double_t Getwire_posZ(Int_t wire) { return wire_posZ[wire]; }
  TVector3 Getwire_pos(Int_t wire) { return TVector3(wire_posX[wire],wire_posY[wire],wire_posZ[wire]); }
  Int_t GetTDCChannel(Int_t wire) { return tdc[wire]; };
  Int_t GetWireIndex(Int_t tdc_ch);
  Double_t GetTimeOffset(Int_t wire) { return offset[wire]; }
  
  ClassDef(Geometry,1) //Geometry class
  
};

#endif 
