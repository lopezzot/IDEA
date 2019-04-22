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
//
//
#include "IDEANewVolume.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4NistManager.hh"

G4LogicalVolume* IDEANewVolume::IDEANewVolume_log=NULL;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

IDEANewVolume::IDEANewVolume(G4RotationMatrix *pRot,
                       const G4ThreeVector &tlate,
                       G4LogicalVolume *pMotherLogical,
                       G4bool pMany,
                       G4int pCopyNo,
                       IDEADetectorConstruction* c)
  :G4PVPlacement(pRot,tlate,
                 new G4LogicalVolume(new G4Box("temp",1,1,1),
                                     G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic"),
                                     "temp",0,0,0),
                 "IDEANewVolumePhysi",pMotherLogical,pMany,pCopyNo),fConstructor(c)
{
 
   
  G4Box* IDEANewVolume_box = new G4Box("IDEANewVolume_box",10.0*cm,10.0*cm,10.0*cm);
 
  IDEANewVolume_log
      = new G4LogicalVolume(IDEANewVolume_box,
                            G4NistManager::Instance()->FindOrBuildMaterial("G4_Al"),
                            "IDEANewVolume_log",0,0,0);
 
 SetLogicalVolume(IDEANewVolume_log); 
}
