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
// $Id: GMCG4DetectorConstruction.cc 87359 2014-12-01 16:04:27Z gcosmo $
//
/// \file GMCG4DetectorConstruction.cc
/// \brief Implementation of the GMCG4DetectorConstruction class

// just to remove a Warning due to a similar definition bettwen a ROOT variable and a CLHEP one
// this warnign came out after splitting the code
#pragma GCC diagnostic ignored "-Wshadow"
 
// GMC includes
#include <ConstructMaterials.hh>
#include <GeomService.hh>
#include <findMaterialOrThrow.hh>
#include "GMCG4DetectorConstruction.hh"
#include "GMCG4RomaProtoBuilder.hh"
#include "GMCG4LecceProtoBuilder.hh"
#include <GMCG4ModupixTrkrBuilder.hh>
#include "GMCG4DetectorMessenger.hh"
#include "GMCG4FieldSetup.hh"

// Beam Pipe
#include "GMCG4BeamPipeBuilder.hh"

// CDCH includes
#include "CDCHMaker.hh"
//#include "CDCHtracker.hh"
#include <CDCHtrackerBuilder.hh>
#include <SensitiveDetectorName.hh>

// SVX includes
#include "SVXMaker.hh"
#include "SVXtracker.hh"
#include <SVXtrackerBuilder.hh>

// PSHW includes
#include "PSHWMaker.hh"
#include "PSHWtracker.hh"
#include <PSHWBuilder.hh>
#include "PSHWRadiator.hh"

// MEG Target
#include "GMCG4MegTAR.hh"

// G4 includes
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVDivision.hh"
//#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4GDMLParser.hh"
#include "G4SystemOfUnits.hh"

#include "G4ios.hh"

//c++ includes
#include <vector>
#include <string>
#include <algorithm>
#include "GMCG4TrackerSD.hh"
#include "RootIO.hh"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
//G4ThreadLocal
//G4GlobalMagFieldMessenger* GMCG4DetectorConstruction::fMagFieldMessenger = 0;

GMCG4DetectorConstruction::GMCG4DetectorConstruction(G4String fGeomConfName, int runNum, G4String dataOutFold)
: G4VUserDetectorConstruction()
 ,fStepLimit(NULL)
 ,fTheWorld(NULL)
 ,fCheckOverlaps(false)
 ,fMegTar(NULL)
 {
 
  fMessenger = new GMCG4DetectorMessenger(this);
  std::cout<<"Messenger "<<fMessenger<<std::endl;
//  GeomService *geoms = GeomService::Instance("config/geom_01.txt");
//  /*const crd::SimpleConfig &*/cRd = &geoms->getConfig();

  ReadConfigGeomFile(fGeomConfName);

  //Load data saver
  RootIO::GetInstance(runNum,dataOutFold);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
GMCG4DetectorConstruction::~GMCG4DetectorConstruction() {

  if (fMegTar!=NULL) delete fMegTar;

  delete fStepLimit;
  delete fMessenger;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 G4VPhysicalVolume* GMCG4DetectorConstruction::Construct() {

  // Define materials
//  DefineMaterials();
  gmc::ConstructMaterials matConst;
  matConst.construct();

  //Create the world
  CreateWorld();

  // Define volumes
  DefineVolumes();

  // Visualization attributes
  G4VisAttributes* boxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  fTheWorld->GetLogicalVolume()->SetVisAttributes(boxVisAtt);

  // Set step limit from configuration
  SetMaxStep(0.);
  
  return GetWorld();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GMCG4DetectorConstruction::CreateWorld() {

  G4Material* air  = gmc::findMaterialOrThrow(cRd->getString("world.materialName"));

  //G4double worldLength = cRd->getDouble("world.Length") ;
  std::vector<double> worldDims;
  cRd->getVectorDouble("world.dims",worldDims);

  // Definitions of Solids, Logical Volumes, Physical Volumes

  // World
  G4GeometryManager::GetInstance()->SetWorldMaximumExtent(*std::max_element(worldDims.begin(),worldDims.end())/*worldLength*/);
  if (fCheckOverlaps) {
    G4cout << "Computed tolerance = "
        << G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/mm
        << " mm" << G4endl;
  }

  G4Box* worldS
    = new G4Box("world",                                    //its name
                  worldDims[0]/2,worldDims[1]/2,worldDims[2]/2); //its size
//                worldLength/2,worldLength/2,worldLength/2); //its size

  G4LogicalVolume* worldLV
    = new G4LogicalVolume(
                 worldS,   //its solid
                 air,      //its material
                 "World"); //its name
  
  //  Must place the World Physical volume unrotated at (0,0,0).
  // 

  fTheWorld = new G4PVPlacement(
                 0,               // no rotation
                 G4ThreeVector(), // at (0,0,0)
                 worldLV,         // its logical volume
                 "World",         // its name
                 0,               // its mother  volume
                 false,           // no boolean operations
                 0,               // copy number
                 fCheckOverlaps); // checking overlaps 

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GMCG4DetectorConstruction::DefineVolumes() {

  ConstructBeamPipe();

  ConstructMEGTarget();

  ConstructCentralTracker();

  ConstructVertexTracker();

  ConstructPreShower();

  //export geometry in GDML file
  if (cRd->getBool("writeGDML",false)) {
    G4GDMLParser parser;
    parser.Write(cRd->getString("GDMLFileName"), fTheWorld);
  }
  
  // Always return the physical world
  //  return TheWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GMCG4DetectorConstruction::ConstructBeamPipe() {

  if (cRd->getBool("hasBeamPipe",false)) GMCG4BeamPipeBuilder::constructPipe(fTheWorld->GetLogicalVolume());

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GMCG4DetectorConstruction::ConstructMEGTarget() {

  if (cRd->getBool("hasMEGTar",false)) {
    fMegTar = new GMCG4MegTAR();
    fMegTar->construct(fTheWorld->GetLogicalVolume());
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GMCG4DetectorConstruction::ConstructCentralTracker() {

  if (cRd->getBool("hasRMproto",false)) GMCG4RomaProtoBuilder::constructTracker(fTheWorld->GetLogicalVolume());
  else if (cRd->getBool("hasLEproto",false)) GMCG4LecceProtoBuilder::constructTracker(fTheWorld->GetLogicalVolume());

  if (cRd->getBool("hasCDCH",false)) {

    RootIO::GetInstance()->CreateMCStepBranches(SensitiveDetectorName::TrackerGas(),"HitsStepCh");
    if (cRd->getBool("cdch.ActiveWiresSD",false)) {
      RootIO::GetInstance()->CreateMCStepBranches(SensitiveDetectorName::TrackerSWires(),"HitsStepChSW");
      if(cRd->getBool("cdch.ActiveFWiresSD",false)) {
        RootIO::GetInstance()->CreateMCStepBranches(SensitiveDetectorName::TrackerFWires(),"HitsStepChFW");
      }
    }

    cdch::CDCHMaker cdchtm( *cRd );
    GeomService::Instance()->addDetector( cdchtm.getCDCHtrackerPtr() );

    cdch::CDCHtrackerBuilder::instantiateSensitiveDetectors("DriftTrackerHitsCollection");
    cdch::CDCHtrackerBuilder::constructTracker(fTheWorld->GetLogicalVolume());
  }


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GMCG4DetectorConstruction::ConstructVertexTracker() {

  if (cRd->getBool("hasMdPvertex",false)) GMCG4ModupixTrkrBuilder::constructVrtxTracker(fTheWorld->GetLogicalVolume());

  if (cRd->getBool("hasSVX",false)) {

    RootIO::GetInstance()->CreateMCStepBranches(SensitiveDetectorName::SVXTrackerRO(),"SVXHitsStepCh");

    svx::SVXMaker svxtm( *cRd );
    GeomService::Instance()->addDetector( svxtm.getSVXtrackerPtr() );

    svx::SVXtrackerBuilder::instantiateSensitiveDetectors("SVXTrackerHitsCollection");
    svx::SVXtrackerBuilder::constructTracker(fTheWorld->GetLogicalVolume());
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GMCG4DetectorConstruction::ConstructPreShower() {

  if (cRd->getBool("hasPSHW",false)) {

    RootIO::GetInstance()->CreateMCStepBranches(SensitiveDetectorName::PSHWTrackerRO(),"PSHWHitsStepCh");

    pshw::PSHWMaker pshwtm( *cRd );
    GeomService::Instance()->addDetector( pshwtm.getPSHWtrackerPtr() );

    pshw::PSHWBuilder::instantiateSensitiveDetectors("PSHWTrackerHitsCollection");
    VolumeInfo pshwvolinf = pshw::PSHWBuilder::constructTracker( fTheWorld->GetLogicalVolume() );

    GeomService::Instance()->addDetector( pshwtm.getPSHWradiatorPtr() );
    pshw::PSHWBuilder::constructRadiator( pshwvolinf.logical );

  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GMCG4DetectorConstruction::ConstructSDandField() {
  // Sensitive detectors

  if ( cRd->getBool("hasRMproto",false) || cRd->getBool("hasLEproto",false) ) {
    RootIO::GetInstance()->CreateMCStepBranches("TrackerChamberSD","HitsStepCh");
    G4String trackerChamberSDname = "GMCG4/TrackerChamberSD";
    GMCG4TrackerSD* aDriftTracker = new GMCG4TrackerSD(trackerChamberSDname,
        "DriftTrackerHitsCollection");
    aDriftTracker->SetDetectorVersion(cRd->getBool("hasLEproto",false));

    // Setting aTrackerSD to all logical volumes with the same name
    // of "Chamber_LV".
    SetSensitiveDetector("tracker_log", aDriftTracker, true);
  }
//  else if (cRd->getBool("hasCDCH",false)) {
//    cdch::CDCHtrackerBuilder::instantiateSensitiveDetectors("DriftTrackerHitsCollection");
//  }

  if ( cRd->getBool("hasMdPvertex",false) ) {
    RootIO::GetInstance()->CreateMCStepBranches("PixelTrackerSD","HitsStepPx");
    G4String pixeltrackerSDname = "GMCG4/PixelTrackerSD";
    GMCG4TrackerSD* aPixelTracker = new GMCG4TrackerSD(pixeltrackerSDname,
        "PixelTrackerHitsCollection");
    SetSensitiveDetector("modupix_log", aPixelTracker, true);
  }

  // Construct the field creator - this will register the field it creates
  if (!fEmFieldSetup.Get()) {
    GMCG4FieldSetup* fieldSetup
      = new GMCG4FieldSetup(G4ThreeVector()/* 3.3*tesla, 0.0, 0.0 )*/ );
    G4AutoDelete::Register(fieldSetup); // Kernel will delete the F01FieldSetup
    fEmFieldSetup.Put(fieldSetup);
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GMCG4DetectorConstruction::ReadConfigGeomFile(G4String fGeomConfName) {

  GeomService *geoms = GeomService::Instance(fGeomConfName);
  cRd = &geoms->getConfig();

  fCheckOverlaps = cRd->getBool("g4.doSurfaceCheck",false);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GMCG4DetectorConstruction::SetMaxStep(G4double maxStep)
{
  if ((fStepLimit==NULL)&&(maxStep>0.)) fStepLimit = new G4UserLimits();
  if ((fStepLimit)&&(maxStep>0.)) {
    fStepLimit->SetMaxAllowedStep(maxStep);
    fTheWorld->GetLogicalVolume()->SetUserLimits(fStepLimit);
  }
  if (cRd->getBool("hasCDCH",false)) { cdch::CDCHtrackerBuilder::constructStepLimiters(); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GMCG4DetectorConstruction::SetCheckOverlaps(G4bool checkOverlaps)
{
  fCheckOverlaps = checkOverlaps;
}  
