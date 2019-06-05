//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: IDEADetectorConstruction.hh 75215 2013-10-29 16:07:06Z gcosmo $
// 
/// \file IDEADetectorConstruction.hh
/// \brief Definition of the IDEADetectorConstruction class

#ifndef IDEADetectorConstruction_h
#define IDEADetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4Material.hh"


class G4Material;
class G4VPhysicalVolume;
class G4GlobalMagFieldMessenger;
class G4LogicalVolume;

/// Detector construction class to define materials and geometry.
/// The calorimeter is a box made of a given number of layers. A layer consists
/// of an absorber plate and of a detection gap. The layer is replicated.
///
/// Four parameters define the geometry of the calorimeter :
///
/// - the thickness of an absorber plate,
/// - the thickness of a gap,
/// - the number of layers,
/// - the transverse size of the calorimeter (the input face is a square).
///
/// In ConstructSDandField() sensitive detectors of IDEACalorimeterSD type
/// are created and associated with the Absorber and Gap volumes.
/// In addition a transverse uniform magnetic field is defined 
/// via G4GlobalMagFieldMessenger class.

class IDEADetectorConstruction : public G4VUserDetectorConstruction
{
public:
  IDEADetectorConstruction();
  virtual ~IDEADetectorConstruction();
  
public:
  virtual G4VPhysicalVolume* Construct();
  virtual void ConstructSDandField();
  
  void SetNewVolumeOn(G4bool b);
  G4bool GetNewVolumeOn(){return fIDEANewVolumeOn;}
  
private:
  
  // methods
  //
  void DefineMaterials();
  G4VPhysicalVolume* DefineVolumes();
  
  // data members
  //
  static G4ThreadLocal G4GlobalMagFieldMessenger*  fMagFieldMessenger; 
  // magnetic field messenger
  
  G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps
  G4int   fNofLayers;     // number of layers
  G4bool fIDEANewVolumeOn;

  //MATERIALS
  G4Material* AIR;
  G4Material* COP;
  G4Material* KAP;
  G4Material* SiO2;
  G4Material* B2O3;
  G4Material* Al2O3;
  G4Material* CaO;
  G4Material *fglas;
  G4Material *epoxy;
  G4Material *ARGON;
  G4Material *CO2;
  G4Material* VET;
  G4Material* dlc;
  G4Material *CF4;
  G4Material *GAS;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

