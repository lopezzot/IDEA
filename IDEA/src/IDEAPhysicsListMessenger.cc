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
/// \file hadronic/Hadr01/src/PhysicsListMessenger.cc
/// \brief Implementation of the PhysicsListMessenger class
//
//
// $Id: PhysicsListMessenger.cc 70761 2013-06-05 12:30:51Z gcosmo $
//
//
/////////////////////////////////////////////////////////////////////////
//
// PhysicsListMessenger
//
// Created: 31.01.2006 V.Ivanchenko
//
// Modified:
// 04.06.2006 Adoptation of Hadr01 (V.Ivanchenko)
//
////////////////////////////////////////////////////////////////////////
//
// 

#include "IDEAPhysicsListMessenger.hh"

#include "IDEAPhysicsList.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UImanager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

IDEAPhysicsListMessenger::IDEAPhysicsListMessenger(IDEAPhysicsList* pPhys)
:G4UImessenger(), fIDEAPhysicsList(pPhys),
 fGammaCutCmd(0), fElectCutCmd(0), fPosCutCmd(0), fCutCmd(0), fAllCutCmd(0),
 fPListCmd(0), fListCmd(0) 
{   
  fGammaCutCmd = new G4UIcmdWithADoubleAndUnit("/IDEA/CutGamma",this);  
  fGammaCutCmd->SetGuidance("Set gamma cut.");
  fGammaCutCmd->SetParameterName("Gcut",false);
  fGammaCutCmd->SetUnitCategory("Length");
  fGammaCutCmd->SetRange("Gcut>=0.0");
  fGammaCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fElectCutCmd = new G4UIcmdWithADoubleAndUnit("/IDEA/CutEl",this);  
  fElectCutCmd->SetGuidance("Set electron cut.");
  fElectCutCmd->SetParameterName("Ecut",false);
  fElectCutCmd->SetUnitCategory("Length");
  fElectCutCmd->SetRange("Ecut>=0.0");
  fElectCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fPosCutCmd = new G4UIcmdWithADoubleAndUnit("/IDEA/CutPos",this);
  fPosCutCmd->SetGuidance("Set positron cut.");
  fPosCutCmd->SetParameterName("Pcut",false);
  fPosCutCmd->SetUnitCategory("Length");
  fPosCutCmd->SetRange("Pcut>=0.0");
  fPosCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fCutCmd = new G4UIcmdWithADoubleAndUnit("/IDEA/CutProt",this);
  fCutCmd->SetGuidance("Set proton cut.");
  fCutCmd->SetParameterName("ProtCut",false);
  fCutCmd->SetUnitCategory("Length");
  fCutCmd->SetRange("ProtCut>=0.0");
  fCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fAllCutCmd = new G4UIcmdWithADoubleAndUnit("/IDEA/CutsAll",this);
  fAllCutCmd->SetGuidance("Set cut for all.");
  fAllCutCmd->SetParameterName("cut",false);
  fAllCutCmd->SetUnitCategory("Length");
  fAllCutCmd->SetRange("cut>=0.0");
  fAllCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fPListCmd = new G4UIcmdWithAString("/IDEA/Physics",this);
  fPListCmd->SetGuidance("Add modula physics list.");
  fPListCmd->SetParameterName("PList",false);
  fPListCmd->AvailableForStates(G4State_PreInit);

  fListCmd = new G4UIcmdWithoutParameter("/IDEA/ListPhysics",this);
  fListCmd->SetGuidance("Available Physics Lists");
  fListCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

IDEAPhysicsListMessenger::~IDEAPhysicsListMessenger()
{
  delete fGammaCutCmd;
  delete fElectCutCmd;
  delete fPosCutCmd;
  delete fCutCmd;
  delete fAllCutCmd;
  delete fPListCmd;
  delete fListCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void IDEAPhysicsListMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if( command == fGammaCutCmd ) {
    if(fIDEAPhysicsList) {
      fIDEAPhysicsList->SetCutForGamma(fGammaCutCmd->GetNewDoubleValue(newValue));
    } else {
      UI->ApplyCommand("/run/setCutForAGivenParticle gamma " + newValue);
    }

  } else if( command == fElectCutCmd ) {
    if(fIDEAPhysicsList) {
      fIDEAPhysicsList->SetCutForElectron(
        fElectCutCmd->GetNewDoubleValue(newValue));
    } else {
      UI->ApplyCommand("/run/setCutForAGivenParticle e- " + newValue);
    }

  } else if( command == fPosCutCmd ) {
    if(fIDEAPhysicsList) {
      fIDEAPhysicsList->SetCutForPositron(fPosCutCmd->GetNewDoubleValue(newValue));
    } else {
      UI->ApplyCommand("/run/setCutForAGivenParticle e+ " + newValue);
    }

  } else if( command == fCutCmd ) {
    if(fIDEAPhysicsList) {
      fIDEAPhysicsList->SetCutForProton(fCutCmd->GetNewDoubleValue(newValue));
    } else {
      UI->ApplyCommand("/run/setCutForAGivenParticle proton " + newValue);
    }

  } else if( command == fAllCutCmd ) {

    if(fIDEAPhysicsList) {
      G4double cut = fAllCutCmd->GetNewDoubleValue(newValue);
      fIDEAPhysicsList->SetCutForGamma(cut);
      fIDEAPhysicsList->SetCutForElectron(cut);
      fIDEAPhysicsList->SetCutForPositron(cut);
      fIDEAPhysicsList->SetCutForProton(cut);
    } else {
      UI->ApplyCommand("/run/setCut " + newValue);
    }

  } else if( command == fPListCmd ) {
    if(fIDEAPhysicsList) {
      G4String name = newValue;
      if(name == "PHYSLIST") {
        char* path = getenv(name);
        if (path) name = G4String(path);
        else {
          G4cout << "### IDEAPhysicsListMessenger WARNING: "
                 << " environment variable PHYSLIST is not defined"
                 << G4endl;
          return; 
        }
      }
      fIDEAPhysicsList->AddIDEAPhysicsList(name);
    } else {
      G4cout << "### IDEAPhysicsListMessenger WARNING: "
             << " /IDEA/Physics UI command is not available "
             << "for reference Physics List" << G4endl;
    }

  } else if( command == fListCmd ) {
    if(fIDEAPhysicsList) {
      fIDEAPhysicsList->List();
    } else { 
      G4cout << "### IDEAPhysicsListMessenger WARNING: "
             << " /IDEA/ListPhysics UI command is not available "
             << "for reference Physics List" << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
