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
// $Id: IDEADetectorConstruction.cc 101905 2016-12-07 11:34:39Z gunter $
// 
/// \file IDEADetectorConstruction.cc
/// \brief Implementation of the IDEADetectorConstruction class

#include "IDEADetectorConstruction.hh"
#include "IDEAPreShowerSD.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4SDManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

// include a New Volume
#include "IDEANewVolume.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal 
G4GlobalMagFieldMessenger* IDEADetectorConstruction::fMagFieldMessenger = 0; 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

IDEADetectorConstruction::IDEADetectorConstruction()
 : G4VUserDetectorConstruction(),
   fCheckOverlaps(true),
   fNofLayers(-1)
{
// if you want a new volume
   fIDEANewVolumeOn=true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

IDEADetectorConstruction::~IDEADetectorConstruction()
{ 
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* IDEADetectorConstruction::Construct()
{
  // Define materials 
  DefineMaterials();
  
  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void IDEADetectorConstruction::DefineMaterials()
{ 
  // Lead material defined using NIST Manager
  auto nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_Pb");
  
  // Liquid argon material
  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;  
  G4double density; 
  new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
         // The argon by NIST Manager is a gas with a different density

  // Vacuum
  new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                  kStateGas, 2.73*kelvin, 3.e-18*pascal);

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* IDEADetectorConstruction::DefineVolumes()
{
    
  // Get materials
  auto defaultMaterial = G4Material::GetMaterial("Galactic");
  auto absorberMaterial = G4Material::GetMaterial("G4_Pb");
  auto gapMaterial = G4Material::GetMaterial("liquidArgon");
  
  if ( ! defaultMaterial || ! absorberMaterial || ! gapMaterial ) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined."; 
    G4Exception("IDEADetectorConstruction::DefineVolumes()",
      "MyCode0001", FatalException, msg);
  }  
   
  // Sizes

  G4double worldR  = 2500.0*mm;
  //
  // World
  //
  auto solidW = new G4Tubs("World",0.,worldR,2.0*worldR,0.,twopi);
  auto fLogicWorld = new G4LogicalVolume( solidW,defaultMaterial,"World");
  auto world = new G4PVPlacement(0,G4ThreeVector(),
                                       fLogicWorld,"World",0,false,0);


  // -------------------------------------------------------------------------
  // Cylinder
  // -------------------------------------------------------------------------
  G4double innerRadius = 1050.*mm;
  G4double outerRadius = 1250.*mm;
  G4double hz = 1250.*mm;
  G4double startAngle = 0.*deg;
  G4double spanningAngle = 360.*deg;
  G4ThreeVector positionPreshower = G4ThreeVector (0., 0., 0.);
  
  auto solidPreshower = new G4Tubs ("solidPreshower", innerRadius, outerRadius, hz, startAngle, spanningAngle);

  auto logicPreshower = new G4LogicalVolume (solidPreshower,defaultMaterial, "logicPreshower", 0, 0, 0);

  auto physPreshower =new G4PVPlacement(0,            
					positionPreshower,  
					logicPreshower,    
					"physPreshower",  
					fLogicWorld,     
					false,          
					0);             


  // -------------------------------------------------------------------------
  // mu-RWELL
  // -------------------------------------------------------------------------
  G4double rwell_x = (500./2.)*mm;
  G4double rwell_y = (20./2.)*mm;
  G4double rwell_z = (500./2.)*mm;
  G4double px, py, pz;
  // 13 sector with 5 chambers along Z-axis:
  // we define the angle in the plane XY occupied by each chamber
  // taking into account a minimum separation between each other (given by the rest divided per 13!)
  G4double rwell_rotAngle = (0.46748636173 + 0.01583558497)*rad;
  G4int rwell_nCopyXY = 13;
  G4int rwell_nCopyZ = 5;
  fNofLayers = rwell_nCopyZ*rwell_nCopyXY;
  
  auto solidRwell = new G4Box ("solidRwell", rwell_x, rwell_y, rwell_z);

  auto logicRwell = new G4LogicalVolume (solidRwell, absorberMaterial, "logicRwell", 0, 0, 0);

  G4RotationMatrix *zRot[13];
  
  for (G4int i=0; i < rwell_nCopyXY; i++)
    {
      zRot[i] = new G4RotationMatrix;
      zRot[i] -> rotateZ (i * rwell_rotAngle);
      px = (innerRadius + rwell_y) * cos((twopi/4.) - i * rwell_rotAngle);   
      py = (innerRadius + rwell_y) * sin((twopi/4.) - i * rwell_rotAngle);   

      for (G4int z=0; z < rwell_nCopyZ; z++)
	{
	  pz = -hz + (2*z+1) * rwell_z;
	  new G4PVPlacement (zRot[i], G4ThreeVector (px, py, pz), logicRwell, "physRwell", logicPreshower, false, i*rwell_nCopyZ + z);
	}
    }

  //------------------
 
  //Place IDEANewVolume
  if(fIDEANewVolumeOn){
    new IDEANewVolume(0,G4ThreeVector(0.,0.,(hz + 200.0)*mm),fLogicWorld,false,0,this);
  }//end of Place IDEANewVolume

  
 

  // colors
  G4VisAttributes zero = G4VisAttributes::GetInvisible();
  fLogicWorld->SetVisAttributes(&zero);
 
  //G4VisAttributes* verdecolor = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0));//verde
  G4VisAttributes* redcolor = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));//red
  //G4VisAttributes* cyancolor = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0));//cyan
  G4VisAttributes* blucolor = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));//blu 
  //G4VisAttributes* giallocolor = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));// giallo 
  //G4VisAttributes* magentacolor = new G4VisAttributes(G4Colour(1.0, 0.0, 1.0));//magenta
  //G4VisAttributes* grigiocolor = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));//grigio

  logicPreshower->SetVisAttributes(blucolor);//era grigio
  logicRwell->SetVisAttributes(redcolor);//era grigio

  
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  return world;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void IDEADetectorConstruction::ConstructSDandField()
{
  // G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

  // 
  // Sensitive detectors
  //
  auto absoSD 
    = new IDEAPreShowerSD("AbsorberSD", "AbsorberHitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(absoSD);
  SetSensitiveDetector("logicRwell",absoSD);

 
  auto gapSD 
    = new IDEAPreShowerSD("GapSD", "GapHitsCollection", 1);
  G4SDManager::GetSDMpointer()->AddNewDetector(gapSD);
  SetSensitiveDetector("logicPreshower",gapSD);
 
  // 
  // Magnetic field
  //
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue;
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);
  
  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
