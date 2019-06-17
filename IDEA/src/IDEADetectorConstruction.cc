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

    // ---------- Material defined by NIST Manager
  AIR =  nistManager->FindOrBuildMaterial("G4_AIR");
  COP =  nistManager->FindOrBuildMaterial("G4_Cu");
  KAP =  nistManager->FindOrBuildMaterial("G4_KAPTON");
  SiO2 = nistManager->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
  B2O3 = nistManager->FindOrBuildMaterial("G4_BORON_OXIDE");
  Al2O3 = nistManager->FindOrBuildMaterial("G4_ALUMINUM_OXIDE");
  CaO   = nistManager->FindOrBuildMaterial("G4_CALCIUM_OXIDE");
  ARGON = nistManager->FindOrBuildMaterial("G4_Ar"); //0.00178 * g/cm3
  CO2 = nistManager->FindOrBuildMaterial("G4_CARBON_DIOXIDE"); //0.00198 * g/cm3

  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;
  G4double density;
  
  // VACUUM
  new G4Material("Galactic", z=1., a=1.01*g/mole,density=universe_mean_density, kStateGas, 2.73*kelvin, 3.e-18*pascal);
  // LIQUID ARGON
  new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3); // The argon by NIST Manager is a gas with a different density


  // ---------- Epoxy, from CGEM-IT code
  G4double epoxy_density = 1.25*g/cm3;
  epoxy =  new G4Material ("epoxy", epoxy_density, 3);
  // H
  z = 1., a = 1.01*g/mole;
  G4Element *H = new G4Element("H", "H", z, a);
  // O
  z=8.; a = 16.0*g/mole;
  G4Element *O = new G4Element("O", "O", z, a);
  // C
  z=6.; a = 12.01*g/mole;
  G4Element* C = new G4Element("C", "C", z, a);
  epoxy->AddElement(C, 18);
  epoxy->AddElement(H, 31);
  epoxy->AddElement(O, 3);
  
  // ---------- Fiberglas, a mean value from tab. 1 in  https://www.asminternational.org/documents/10192/1849770/06781G_p27-34.pdf
  //            mass fraction %:
  //            SiO2  60%
  //            B2O3   5%
  //            Al2O3 13%
  //            CaO   22%
  G4double fglas_density = 1.99*g/cm3;
  fglas =  new G4Material("fglas", fglas_density, 4);
  fglas->AddMaterial(SiO2, 0.6);
  fglas->AddMaterial(B2O3, 0.05);
  fglas->AddMaterial(Al2O3, 0.13);
  fglas->AddMaterial(CaO, 0.22);
  
  // ---------- FR4 vetronite (fiberglas 60%  + epoxy 40%)
  //            Simulated as permaglas (density=1.97*g/cm3) with FR4 density
  G4double vet_density = 1.85*g/cm3;
  VET =  new G4Material("vetronite", vet_density, 2);
  VET->AddMaterial(fglas, 0.6);
  VET->AddMaterial(epoxy, 0.4);
  
  // ---------- DLC - Diamond like carbon
  // ---------- Pre-preg - We assumed the same density for the pre-preg material
  G4double dlc_density = 2.*g/cm3;
  dlc =  new G4Material("DLC", dlc_density, 1);
  dlc->AddElement(C, 1);
  
  // ---------- GAS, Ar:CO2:CF4 45:15:40
  G4double gas_density= 0.002536 * g/cm3; 
  
  G4double cf4_density = 0.00378 * g/cm3;
  a = 12.011*g/mole;
  G4Element* elC = new G4Element("Carbon", "C", z=6., a);  
  a = 18.998*g/mole;
  G4Element* elF = new G4Element("Fluorine", "F", z=9., a);  
  CF4 = new G4Material("G4_CF4", cf4_density, 2);
  CF4->AddElement(elC, 1);
  CF4->AddElement(elF, 4);

  GAS =  new G4Material("ArCO2CF4", gas_density, 3);
  GAS->AddMaterial(ARGON, 0.295); 
  GAS->AddMaterial(CO2, 0.109);   
  GAS->AddMaterial(CF4, 0.596);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* IDEADetectorConstruction::DefineVolumes()
{
  // Get MATERIALS
  auto defaultMaterial = G4Material::GetMaterial("Galactic");
  auto absorberMaterial = G4Material::GetMaterial("G4_Pb");
  auto gapMaterial = G4Material::GetMaterial("liquidArgon");

  if ( ! defaultMaterial || ! absorberMaterial || ! gapMaterial )
    {
      G4ExceptionDescription msg;
      msg << "Cannot retrieve materials already defined."; 
      G4Exception("IDEADetectorConstruction::DefineVolumes()",
      "MyCode0001", FatalException, msg);
    }  

  

  // -------------------------------------------------------------------------
  // World
  // -------------------------------------------------------------------------
  G4double worldR  = 2500.0 * CLHEP::mm;
  auto solidW = new G4Tubs ("World", 0., worldR, 2.0*worldR, 0., twopi);
  auto fLogicWorld = new G4LogicalVolume ( solidW, defaultMaterial, "World");
  auto world = new G4PVPlacement (0, G4ThreeVector(), fLogicWorld, "World", 0, false, 0);


  // -------------------------------------------------------------------------
  // Cylinder
  // -------------------------------------------------------------------------
  G4double innerRadius = 1050. * CLHEP::mm;
  G4double outerRadius = 1250. * CLHEP::mm;
  G4double hz = 1250. * CLHEP::mm;
  G4double startAngle = 0.*deg;
  G4double spanningAngle = 360.*deg;
  G4ThreeVector positionPreshower = G4ThreeVector (0., 0., 0.);
  
  auto solidPreshower = new G4Tubs ("solidPreshower", innerRadius, outerRadius, hz, startAngle, spanningAngle);

  auto logicPreshower = new G4LogicalVolume (solidPreshower, defaultMaterial, "logicPreshower", 0, 0, 0);

  auto physPreshower =new G4PVPlacement(0,            
					positionPreshower,  
					logicPreshower,    
					"physPreshower",  
					fLogicWorld,     
					false,          
					0);             


  // -------------------------------------------------------------------------
  // micro-RWELL
  // -------------------------------------------------------------------------

  // CATHODE
  G4double rwell_cat_fr4_thick  = 1.6 * CLHEP::mm;
  G4double rwell_cat_cop_thick  = 35 * CLHEP::micrometer;
  G4cout << "Cathode thicknesses: " << G4endl;
  G4cout << "vetronite " << rwell_cat_fr4_thick << G4endl;
  G4cout << "copper    " << rwell_cat_cop_thick << G4endl;

  // GAS GAP
  G4double rwell_drift_gap = 6 * CLHEP::mm;

  // micro-RWELL + READOUT
  G4double rwell_pcb_cop_hole_thick  = 5 * CLHEP::micrometer; // w/ hole
  G4double rwell_pcb_kap_thick       = 50 * CLHEP::micrometer; // w/ hole
  G4double rwell_pcb_dlc_thick       = 0.1 * CLHEP::micrometer;
  G4double rwell_pcb_grid_thick      = 35 * CLHEP::micrometer;
  G4double rwell_pcb_prepreg_thick   = 100 * CLHEP::micrometer;
  G4double rwell_pcb_cop_thick       = 35 * CLHEP::micrometer; // w/o hole
  G4double rwell_pcb_fr4_thick       = 1.6 * CLHEP::mm;
  G4cout << "micro-RWELL + readout thicknesses: " << G4endl;
  G4cout << "copper hole  " << rwell_pcb_cop_hole_thick << G4endl;
  G4cout << "kapton hole  " << rwell_pcb_kap_thick << G4endl;
  G4cout << "dlc          " << rwell_pcb_dlc_thick << G4endl;
  G4cout << "grid         " << rwell_pcb_grid_thick << G4endl;
  G4cout << "prepreg      " << rwell_pcb_prepreg_thick << G4endl;
  G4cout << "copper strip " << rwell_pcb_cop_thick << G4endl;
  G4cout << "vetronite    " << rwell_pcb_fr4_thick << G4endl;

  // micro-RWELL CHAMBER COMPONENTS
  G4double rwell_cat_thick   = rwell_cat_fr4_thick + rwell_cat_cop_thick; // CATHODE (vetronite+copper)
  G4double rwell_pcb_thick   = rwell_pcb_cop_hole_thick + rwell_pcb_kap_thick + rwell_pcb_dlc_thick + rwell_pcb_grid_thick + rwell_pcb_prepreg_thick + rwell_pcb_cop_thick + rwell_pcb_fr4_thick; //PCB
  G4double rwell_ch_thick    = rwell_cat_thick + rwell_pcb_thick + rwell_drift_gap; // FULL micro-RWELL CHAMBER

  // CHAMBER DIMENSIONS
  G4double rwell_x_width = 500 * CLHEP::mm;
  G4double rwell_y_width = rwell_ch_thick;
  G4double rwell_z_width = 500 * CLHEP::mm;
  G4double px, py, pz;

  G4cout << "micro-RWELL thicknesses: " << G4endl;
  G4cout << "CATHODE                  " << rwell_cat_thick << G4endl;
  G4cout << "micro-RWELL + readout    " << rwell_pcb_thick << G4endl;
  G4cout << "FULL micro-RWELL chamber " << rwell_ch_thick << '\n' << G4endl;

 
  // ---------- Take into account holes in copper and kapton:
  //            consider one element of ANSYS geometry ===> p=pitch, area=p X sqrt(3)*p, 2 holes
  //            cylinder in copper / double cone in kapton
  //   copper:
  //   5 micron      thickness
  //   70 micron     diameter
  //   140 micron    pitch
  //   kapton:
  //   50 micron     thickness
  //   50-70 micron  diameter
  //   25-35 micron  r-R
  //   140 micron    pitch
  //
  // cylinder in copper/double cone in kapton
  G4double hole_p = 140 * CLHEP::micrometer;
  G4double hole_r = 25 * CLHEP::micrometer;
  G4double hole_R = 35 * CLHEP::micrometer;

  // COPPER: only hole
  G4double V_cop_full = hole_p * hole_p * sqrt(3.) * rwell_pcb_cop_hole_thick;
  G4double V_cop_hole = rwell_pcb_cop_hole_thick * M_PI * hole_R * hole_R;
  G4double V_cop = V_cop_full - 2*V_cop_hole;
  G4double f_cop = V_cop/V_cop_full;
  G4double cop_hole_density = f_cop * COP->GetDensity();
  G4Material* COP_hole = G4NistManager::Instance()->BuildMaterialWithNewDensity("COP_hole", "G4_Cu", cop_hole_density,
		    COP->GetTemperature(), COP->GetPressure());
  G4cout << "--Take into account THE COPPER holes--" << G4endl;
  G4cout << "COPPER density (std) " << COP->GetDensity()/g/cm3 << "; COPPER density (in micro-RWELL) " << COP_hole->GetDensity()/g/cm3 << '\n' << G4endl;
  
  // KAPTON: only hole
  G4double V_kap_full = hole_p * hole_p * sqrt(3.) * rwell_pcb_kap_thick;
  G4double V_kap_hole = (1./3) * rwell_pcb_kap_thick * M_PI * (hole_r*hole_r + hole_R*hole_R + hole_r*hole_R);
  G4double V_kap = V_kap_full - 2*V_kap_hole;
  G4double f_kap = V_kap/V_kap_full;
  G4double kap_hole_density = f_kap * KAP->GetDensity();
  G4Material* KAP_hole = G4NistManager::Instance()->BuildMaterialWithNewDensity("KAP_hole", "G4_KAPTON", kap_hole_density,
                    KAP->GetTemperature(), KAP->GetPressure());
  G4cout << "Volume FULL " << V_kap_full << ", HOLE " << V_kap_hole << " and FINAL " << V_kap << G4endl;
  G4cout << "KAPTON density (std) " << KAP->GetDensity()/g/cm3 << "; KAPTON density (in micro-RWELL)  " << KAP_hole->GetDensity()/g/cm3 << '\n' << G4endl;
  
  // ---------- Take into account dead surfaces in copper and kapton due to the inclusion of the silver grid (HR layout - SG2++):
  //   copper and kapton:
  //   12 mm             pitch
  //   0.6 mm            dead surface
  //   (12 - 0.6)/12     active surface => 95%
  //   (1 - 0.95)        full (=dead) surface => 5% 
  G4double grid_p = 12 * CLHEP::mm;
  G4double dead_s = 0.6 * CLHEP::mm;
  G4double f_dead = (grid_p - dead_s)/grid_p;
  G4double cop_dead_density = f_dead * COP_hole->GetDensity() + (1 - f_dead) * COP->GetDensity();
  G4double kap_dead_density = f_dead * KAP_hole->GetDensity() + (1 - f_dead) * KAP->GetDensity();

  G4Material* COP_deadS = G4NistManager::Instance()->BuildMaterialWithNewDensity("COP_deadS", "G4_Cu", cop_dead_density,
		    COP->GetTemperature(), COP->GetPressure());
  G4Material* KAP_deadS = G4NistManager::Instance()->BuildMaterialWithNewDensity("KAP_deadS", "G4_KAPTON", kap_dead_density,
		    KAP->GetTemperature(), KAP->GetPressure());
  G4cout << "--Take into account the dead surface in the copper and kapton layers due to the presence of the silver grid--" << G4endl;
  G4cout << "COPPER density (std) " <<COP->GetDensity()/g/cm3<< "; COPPER density (holes) " <<COP_hole->GetDensity()/g/cm3<< "; COPPER density (grid) " <<COP_deadS->GetDensity()/g/cm3<<'\n'<< G4endl;
  G4cout << "KAPTON density (std) " <<KAP->GetDensity()/g/cm3<< "; KAPTON density (holes) " <<KAP_hole->GetDensity()/g/cm3<< "; KAPTON density (grid) " <<KAP_deadS->GetDensity()/g/cm3<<'\n'<< G4endl;

  // ---------- Take into account grid grounded between DLC and pre-preg:
  //   copper:
  //   35 micron     thickness
  //   100 micron    size
  //   12 mm         pitch
  G4double grid_s = 100 * CLHEP::micrometer;
  G4double f_grid_cop = grid_s/grid_p;
  G4double pcb_cop_grid_density = f_grid_cop * COP->GetDensity();
  G4Material* COP_grid = G4NistManager::Instance()->BuildMaterialWithNewDensity("COP_grid", "G4_Cu", pcb_cop_grid_density,
		     COP->GetTemperature(), COP->GetPressure());
  G4cout << "--Take into account THE SILVER GRID for the micro-RWELL HR layout--" << G4endl;
  G4cout << "COPPER density (std) " << COP->GetDensity()/g/cm3 << "; COPPER density (grid)  " << COP_grid->GetDensity()/g/cm3 << '\n' << G4endl;

  // ---------- Take into account strips in the electrode cathode:
  //   copper:
  //   35 micron     thickness
  //   250 micron    size
  //   400 micron    pitch
  G4double strip_s = 250 * CLHEP::micrometer;
  G4double strip_p = 400 * CLHEP::micrometer;
  G4double f_pcb_cop = strip_s/strip_p;
  G4double pcb_cop_strip_density = f_pcb_cop * COP->GetDensity();
  G4Material* COP_strip = G4NistManager::Instance()->BuildMaterialWithNewDensity("COP_strip", "G4_Cu", pcb_cop_strip_density,
		     COP->GetTemperature(), COP->GetPressure());
  G4cout << "--Take into account THE ANODE STRIPS in the PCB electrode--" << G4endl;
  G4cout << "COPPER density (std) " << COP->GetDensity()/g/cm3 << "; COPPER density (strips)  " << COP_strip->GetDensity()/g/cm3 << '\n' << G4endl;

  
  // ---------- CHAMBER ----------
  // solid
  auto rwell_ch_solid = new G4Box("rwell_ch_solid", rwell_x_width * 0.5,  rwell_y_width * 0.5, rwell_z_width * 0.5);
  // logic
  auto rwell_ch_logic = new G4LogicalVolume(rwell_ch_solid, GAS, "rwell_ch_logic", 0, 0, 0);

  // ---------- CATHODE ----------
  // solid
  auto rwell_cat_fr4_solid = new G4Box("rwell_cat_fr4_solid", rwell_x_width * 0.5,  (rwell_cat_fr4_thick + rwell_cat_cop_thick) * 0.5, rwell_z_width * 0.5);
  auto rwell_cat_cop_solid = new G4Box("rwell_cat_cop_solid", rwell_x_width * 0.5, rwell_cat_cop_thick * 0.5, rwell_z_width * 0.5);
  // logic
  auto rwell_cat_logic = new G4LogicalVolume(rwell_cat_fr4_solid, VET, "rwell_cat_logic", 0, 0, 0);
  auto rwell_cat_cop_logic = new G4LogicalVolume(rwell_cat_cop_solid, COP, "rwell_cat_cop_logic", 0, 0, 0);

  // ------------------ micro-RWELL + READOUT ------------------
  //   5     micron rame        --> hole den
  //   50    micron kapton      --> hole den
  //   0.1   micron dlc         --> full                 
  //   35    micron grid        --> cop grid den
  //   100   micron prepreg     --> full              
  //   35    micron copper      --> strip den
  //   1.6   mm vetronite       --> full
  // -----------------------------------------------------------
  // solid
  auto rwell_pcb_cop_hole_solid = new G4Box("rwell_pcb_cop_hole_solid", rwell_x_width * 0.5, rwell_pcb_cop_hole_thick * 0.5, rwell_z_width * 0.5); 
  auto rwell_pcb_kap_solid      = new G4Box("rwell_pcb_kap_solid", rwell_x_width * 0.5, rwell_pcb_kap_thick * 0.5, rwell_z_width * 0.5);
  auto rwell_pcb_dlc_solid      = new G4Box("rwell_pcb_dlc_solid", rwell_x_width * 0.5, rwell_pcb_dlc_thick * 0.5, rwell_z_width * 0.5); 
  auto rwell_pcb_grid_solid     = new G4Box("rwell_pcb_grid_solid", rwell_x_width * 0.5, rwell_pcb_grid_thick * 0.5, rwell_z_width * 0.5); 
  auto rwell_pcb_prepreg_solid  = new G4Box("rwell_pcb_prepreg_solid", rwell_x_width * 0.5, rwell_pcb_prepreg_thick * 0.5, rwell_z_width * 0.5); 
  auto rwell_pcb_cop_solid      = new G4Box("rwell_pcb_cop_solid", rwell_x_width * 0.5, rwell_pcb_cop_thick * 0.5, rwell_z_width * 0.5); 
  auto rwell_pcb_fr4_solid      = new G4Box("rwell_pcb_fr4_solid", rwell_x_width * 0.5, rwell_pcb_thick * 0.5, rwell_z_width * 0.5);
  // logic
  auto rwell_pcb_cop_hole_logic = new G4LogicalVolume(rwell_pcb_cop_hole_solid, COP_deadS, "rwell_pcb_cop_hole_logic", 0, 0, 0);
  auto rwell_pcb_kap_logic      = new G4LogicalVolume(rwell_pcb_kap_solid, KAP_deadS, "rwell_pcb_kap_logic", 0, 0, 0);
  auto rwell_pcb_dlc_logic      = new G4LogicalVolume(rwell_pcb_dlc_solid, dlc, "rwell_pcb_dlc_logic", 0, 0, 0);
  auto rwell_pcb_grid_logic     = new G4LogicalVolume(rwell_pcb_grid_solid, COP_grid, "rwell_pcb_grid_logic", 0, 0, 0);
  auto rwell_pcb_prepreg_logic  = new G4LogicalVolume(rwell_pcb_prepreg_solid, dlc, "rwell_pcb_prepreg_logic", 0, 0, 0);
  auto rwell_pcb_cop_logic      = new G4LogicalVolume(rwell_pcb_cop_solid, COP_strip, "rwell_pcb_cop_logic", 0, 0, 0);
  auto rwell_pcb_logic          = new G4LogicalVolume(rwell_pcb_fr4_solid, VET, "rwell_pcb_logic", 0, 0, 0);


  // ---------- BUILD CATHODE ----------
  // It looks at the interaction point:
  // consequently, we place the copper layer up and the fr4 layer down
  double rwell_cat_cop_y = rwell_cat_thick * 0.5 - rwell_cat_cop_thick * 0.5;
  auto rwell_cat_cop_phys = new G4PVPlacement(0, G4ThreeVector(0, rwell_cat_cop_y, 0),
                rwell_cat_cop_logic, "rwell_cat_cop_phys", rwell_cat_logic, false, 0);
  G4cout << "============================" << G4endl;
  G4cout << "Cathode y position:\t" << G4endl;
  G4cout << "copper\t" << rwell_cat_cop_y << G4endl;
   
  // ------------------ BUILD micro-RWELL + readout -----------------  
  //   5     micron rame    --> hole den + "dead strips" --------
  //   50    micron kapton  --> hole den + "dead strips" --------
  //   0.1   micron dlc     --> full                     --------
  //   35    micron grid    --> cop grid den             --------
  //   100   micron prepreg --> full                     --------
  //   35    micron copper  --> full                     --------
  //   1.6   mm vetronite   --> full                     --------
  // ----------------------------------------------------------------
  double rwell_pcb_cop_hole_y = - rwell_pcb_thick * 0.5 + rwell_pcb_cop_hole_thick * 0.5; 
  double rwell_pcb_kap_y      = rwell_pcb_cop_hole_y + rwell_pcb_cop_hole_thick * 0.5 + rwell_pcb_kap_thick * 0.5; 
  double rwell_pcb_dlc_y      = rwell_pcb_kap_y + rwell_pcb_kap_thick * 0.5 + rwell_pcb_dlc_thick * 0.5;
  double rwell_pcb_grid_y     = rwell_pcb_dlc_y + rwell_pcb_dlc_thick * 0.5 + rwell_pcb_grid_thick * 0.5;
  double rwell_pcb_prepreg_y  = rwell_pcb_grid_y + rwell_pcb_grid_thick * 0.5 + rwell_pcb_prepreg_thick * 0.5;
  double rwell_pcb_cop_y      = rwell_pcb_prepreg_y + rwell_pcb_prepreg_thick * 0.5 + rwell_pcb_cop_thick * 0.5;

  G4cout << "============================" << G4endl;
  G4cout << "PCB y position:\t " << G4endl;
  G4cout << "copper\t" << rwell_pcb_cop_y << G4endl;
  G4cout << "prepreg\t" << rwell_pcb_prepreg_y << G4endl;
  G4cout << "grid\t" << rwell_pcb_grid_y << G4endl;
  G4cout << "dlc\t" << rwell_pcb_dlc_y << G4endl;
  G4cout << "kapton\t" << rwell_pcb_kap_y << G4endl;
  G4cout << "copper\t" << rwell_pcb_cop_hole_y << G4endl;
  G4cout << "============================\n" << G4endl;

  auto rwell_pcb_cop_hole_phys = new G4PVPlacement(0, G4ThreeVector(0, rwell_pcb_cop_hole_y, 0),
                 rwell_pcb_cop_hole_logic, "rwell_pcb_cop_hole_phys", rwell_pcb_logic, false, 0);
  auto rwell_pcb_kap_phys = new G4PVPlacement(0, G4ThreeVector(0, rwell_pcb_kap_y, 0),
                 rwell_pcb_kap_logic, "rwell_pcb_kap_phys", rwell_pcb_logic, false, 0);
  auto rwell_pcb_dlc_phys = new G4PVPlacement(0, G4ThreeVector(0, rwell_pcb_dlc_y, 0),
                 rwell_pcb_dlc_logic, "rwell_pcb_dlc_phys", rwell_pcb_logic, false, 0);
  auto rwell_pcb_grid_phys = new G4PVPlacement(0, G4ThreeVector(0, rwell_pcb_grid_y, 0),
                 rwell_pcb_grid_logic, "rwell_pcb_grid_phys", rwell_pcb_logic, false, 0);
  auto rwell_pcb_prepreg_phys = new G4PVPlacement(0, G4ThreeVector(0, rwell_pcb_prepreg_y, 0),
                 rwell_pcb_prepreg_logic, "rwell_pcb_prepreg_phys", rwell_pcb_logic, false, 0);
  auto rwell_pcb_cop_phys = new G4PVPlacement(0, G4ThreeVector(0, rwell_pcb_cop_y, 0),
                 rwell_pcb_cop_logic, "rwell_pcb_cop_phys", rwell_pcb_logic, false, 0);
  
  // -----------------------------------
  // -- PLACEMENT WITHIN THE CHAMBER  -- 
  // -----------------------------------	

  // ---------- PLACE CATHODE ----------
  // It looks at the interaction point
  double rwell_cat_y = - rwell_ch_thick * 0.5 + rwell_cat_thick * 0.5;
  auto rwell_cat_phys = new G4PVPlacement(0, G4ThreeVector(0, rwell_cat_y, 0),
                rwell_cat_logic, "rwell_cat_phys", rwell_ch_logic, false, 0);

  // ------------ PLACE micro-RWELL + readout ------------  
  double rwell_pcb_y = rwell_ch_thick * 0.5 - rwell_pcb_thick * 0.5;
  auto rwell_pcb_phys = new G4PVPlacement(0, G4ThreeVector(0, rwell_pcb_y, 0),
                 rwell_pcb_logic, "rwell_pcb_phys", rwell_ch_logic, false, 0);

  G4cout << "===== PLACEMENT =====" << G4endl;
  G4cout << "Chamber size: \t x = " << rwell_x_width << "\t y = " << rwell_y_width << "\t z = " << rwell_z_width << G4endl;
  G4cout << "CAT: \t y = " << rwell_cat_y << "\n\t copper_y = " << rwell_cat_cop_y << G4endl; 
  G4cout << "PCB: \t y = " << rwell_pcb_y << " \n\t top cop = " << rwell_pcb_cop_hole_y << "\n\t kapton = " 
      << rwell_pcb_kap_y << "\n\t dlc = " << rwell_pcb_dlc_y << "\n\t grid = " << rwell_pcb_grid_y << "\n\t prepreg = " << rwell_pcb_prepreg_y << "\n\t cop = " << rwell_pcb_cop_y << G4endl;

  G4cout << "DISTANCES" << G4endl;
  double rwell_cat_to_pcb = rwell_pcb_y - rwell_cat_y - (rwell_cat_thick + rwell_pcb_thick) * 0.5;
  G4cout << "Drift gap " << rwell_cat_to_pcb << '\n' << G4endl;


  // -----------------------------------
  // ---- PLACEMENT OF THE CHAMBERS ---- 
  // -----------------------------------	
  
  // 13 sector with 5 chambers along Z-axis:
  // We define the angle in the plane XY occupied by each chamber
  // taking into account a minimum separation between each other (given by the rest divided per 13!)
  G4double rwell_rotAngle = (0.46748636173 + 0.01583558497)*rad;
  G4int rwell_nCopyXY = 13;
  G4int rwell_nCopyZ = 5;
  fNofLayers = rwell_nCopyZ*rwell_nCopyXY;
  
  //auto solidRwell = new G4Box ("solidRwell", rwell_x_width * 0.5, rwell_y_width * 0.5, rwell_z_width * 0.5);
  //auto logicRwell = new G4LogicalVolume (solidRwell, absorberMaterial, "logicRwell", 0, 0, 0);
  G4RotationMatrix *zRot[13];
  
  for (G4int i=0; i < rwell_nCopyXY; i++)
    {
      zRot[i] = new G4RotationMatrix;
      zRot[i] -> rotateZ (i * rwell_rotAngle);
      px = (innerRadius + rwell_y_width * 0.5) * cos((twopi/4.) - i * rwell_rotAngle);   
      py = (innerRadius + rwell_y_width * 0.5) * sin((twopi/4.) - i * rwell_rotAngle);   

      for (G4int j=0; j < rwell_nCopyZ; j++)
	{
	  pz = -hz + (2*j+1) * rwell_z_width * 0.5;
	  new G4PVPlacement (zRot[i], G4ThreeVector (px, py, pz), rwell_ch_logic, "rwell_ch_phys", logicPreshower, false, i*rwell_nCopyZ + j);
	}
    }

  //------------------
 
  //Place IDEANewVolume
  if(fIDEANewVolumeOn)
    {
      new IDEANewVolume(0,G4ThreeVector(0.,0.,(hz + 200.0)*mm),fLogicWorld,false,0,this);
    }//end of Place IDEANewVolume

  // ----------------
  // ---- COLORS ---- 
  // ----------------	  
  
  fLogicWorld->SetVisAttributes(G4VisAttributes::Invisible); 
  
  G4VisAttributes* green   = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0)); 
  G4VisAttributes* red     = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0)); 
  G4VisAttributes* cyan    = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0)); 
  G4VisAttributes* blue    = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));  
  G4VisAttributes* yellow  = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));  
  G4VisAttributes* grey    = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5)); 
  G4VisAttributes* magenta = new G4VisAttributes(G4Colour(1.0, 0.0, 1.0)); 

  //  vis fr4 cat
  rwell_cat_logic->SetVisAttributes(grey);
  //  vis cop cat
  rwell_cat_logic->SetVisAttributes(red);
  
  // vis gas
  rwell_ch_logic->SetVisAttributes(cyan);

  // vis cop w/ hole
  rwell_pcb_cop_hole_logic->SetVisAttributes(red);
  // vis kap
  rwell_pcb_kap_logic->SetVisAttributes(yellow);
  // vis dlc
  rwell_pcb_dlc_logic->SetVisAttributes(green);
  // vis grid
  rwell_pcb_grid_logic->SetVisAttributes(red);
  // vis prepreg
  rwell_pcb_prepreg_logic->SetVisAttributes(green);
  // vis cop 
  rwell_pcb_cop_logic->SetVisAttributes(red);
  // vis fr4 pcb
  rwell_pcb_logic->SetVisAttributes(grey);
 
  logicPreshower->SetVisAttributes(blue);
  
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
  SetSensitiveDetector("rwell_ch_logic",absoSD);

 
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
