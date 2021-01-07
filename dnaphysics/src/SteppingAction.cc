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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software 
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// $Id$
//
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class

#include "Analysis.hh"
#include "G4ios.hh"
#include "SteppingAction.hh"
#include "B1RunAction.hh"
#include "G4RunManager.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4SteppingManager.hh"
#include "G4TransportationManager.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4Gamma.hh"
#include "G4Alpha.hh"
#include "G4DNAGenericIonsManager.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"
#include "B1EventAction.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(B1EventAction* eventAction)
: G4UserSteppingAction(),
  fScoringVolume(0),
  fEventAction(eventAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{ 
  G4double flagParticle=-1.;
  G4double flagProcess=-1.;
  G4double x,y,z,xp,yp,zp;

  // Process sub-types are listed in G4PhysicsListHelper.cc

  /*
 // The following method avoids the usage of string comparison 

 if (step->GetTrack()->GetDynamicParticle()->GetDefinition()
     == G4Electron::ElectronDefinition())
    flagParticle = 1; 

 if (step->GetTrack()->GetDynamicParticle()->GetDefinition()
 == G4Proton::ProtonDefinition())
    flagParticle = 2; 

 if (step->GetTrack()->GetDynamicParticle()->GetDefinition()
 == G4Alpha::AlphaDefinition())
    flagParticle = 4; 

    G4DNAGenericIonsManager *instance;
    instance = G4DNAGenericIonsManager::Instance();
    G4ParticleDefinition* protonDef = G4Proton::ProtonDefinition();
    G4ParticleDefinition* hydrogenDef = instance->GetIon("hydrogen");
    G4ParticleDefinition* alphaPlusPlusDef = instance->GetIon("alpha++");
    G4ParticleDefinition* alphaPlusDef = instance->GetIon("alpha+");
    G4ParticleDefinition* heliumDef = instance->GetIon("helium");

 if (step->GetTrack()->GetDynamicParticle()->GetDefinition()
 == instance->GetIon("hydrogen"))
    flagParticle = 3; 

 if (step->GetTrack()->GetDynamicParticle()->GetDefinition()
 == instance->GetIon("alpha+"))
    flagParticle = 5; 

 if (step->GetTrack()->GetDynamicParticle()->GetDefinition()
 == instance->GetIon("helium"))
    flagParticle = 6; 
   */

  const G4String& particleName = step->GetTrack()->GetDynamicParticle()->
      GetDefinition()->GetParticleName();
  const G4String& processName = step->GetPostStepPoint()->
      GetProcessDefinedStep()->GetProcessName();

  // Particle
  if (particleName == "e-")       flagParticle = 1;
  else if (particleName == "proton")   flagParticle = 2;
  else if (particleName == "hydrogen") flagParticle = 3;
  else if (particleName == "alpha")    flagParticle = 4;
  else if (particleName == "alpha+")   flagParticle = 5;
  else if (particleName == "helium")   flagParticle = 6;

  // Processes
  if (processName=="e-_G4DNAElastic")    flagProcess =11;
  else if (processName=="e-_G4DNAExcitation")    flagProcess =12;
  else if (processName=="e-_G4DNAIonisation")    flagProcess =13;
  else if (processName=="e-_G4DNAAttachment")    flagProcess =14;
  else if (processName=="e-_G4DNAVibExcitation")  flagProcess =15;

  else if (processName=="proton_G4DNAExcitation")  flagProcess =17;
  else if (processName=="proton_G4DNAIonisation")  flagProcess =18;
  else if (processName=="proton_G4DNAChargeDecrease")  flagProcess =19;

  else if (processName=="hydrogen_G4DNAExcitation")    flagProcess =20;
  else if (processName=="hydrogen_G4DNAIonisation")    flagProcess =21;
  else if (processName=="hydrogen_G4DNAChargeIncrease")  flagProcess =22;

  else if (processName=="alpha_G4DNAExcitation")    flagProcess =23;
  else if (processName=="alpha_G4DNAIonisation")    flagProcess =24;
  else if (processName=="alpha_G4DNAChargeDecrease")    flagProcess =25;

  else if (processName=="alpha+_G4DNAExcitation")  flagProcess =26;
  else if (processName=="alpha+_G4DNAIonisation")  flagProcess =27;
  else if (processName=="alpha+_G4DNAChargeDecrease")  flagProcess =28;
  else if (processName=="alpha+_G4DNAChargeIncrease")  flagProcess =29;

  else if (processName=="helium_G4DNAExcitation")  flagProcess =30;
  else if (processName=="helium_G4DNAIonisation")  flagProcess =31;
  else if (processName=="helium_G4DNAChargeIncrease")  flagProcess =32;

  if (processName!="Transportation")
  {
    x=step->GetPreStepPoint()->GetPosition().x()/nanometer;
    y=step->GetPreStepPoint()->GetPosition().y()/nanometer;
    z=step->GetPreStepPoint()->GetPosition().z()/nanometer;
    xp=step->GetPostStepPoint()->GetPosition().x()/nanometer;
    yp=step->GetPostStepPoint()->GetPosition().y()/nanometer;
    zp=step->GetPostStepPoint()->GetPosition().z()/nanometer;

// get analysis manager
//    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
/*
    // fill ntuple
    analysisManager->FillNtupleDColumn(0, flagParticle);
    analysisManager->FillNtupleDColumn(1, flagProcess);
    analysisManager->FillNtupleDColumn(2, x);
    analysisManager->FillNtupleDColumn(3, y);
    analysisManager->FillNtupleDColumn(4, z);
    analysisManager->FillNtupleDColumn(5, step->GetTotalEnergyDeposit()/eV);
    analysisManager->FillNtupleDColumn(6,
                                       std::sqrt((x-xp)*(x-xp)+
                                           (y-yp)*(y-yp)+(z-zp)*(z-zp))/nm);
    analysisManager->FillNtupleDColumn(7,
                                       (step->GetPreStepPoint()->
                                           GetKineticEnergy() -
                                           step->GetPostStepPoint()->
                                           GetKineticEnergy())/eV );
    analysisManager->FillNtupleIColumn(8, G4EventManager::GetEventManager()->
                                       GetConstCurrentEvent()->GetEventID());
    analysisManager->AddNtupleRow();*/
  }
G4StepPoint* PreStep = step->GetPreStepPoint();
G4StepPoint* PostStep = step->GetPostStepPoint();

//G4double PreStepX = PreStep->GetPosition().x();
//G4double PreStepY = PreStep->GetPosition().y();
//G4double PreStepZ = PreStep->GetPosition().z();
G4double parentID = step->GetTrack()->GetParentID();
G4double trackID = step->GetTrack()->GetTrackID();

G4double PostStepX = PostStep->GetPosition().x();
G4double PostStepY = PostStep->GetPosition().y();
G4double PostStepZ = PostStep->GetPosition().z();

// positions in the global coordinate system:
//G4ThreeVector posPreStep  = PreStep->GetPosition();
// G4ThreeVector posPostStep = PostStep->GetPosition();

G4TouchableHandle touchPreStep = PreStep->GetTouchableHandle();
G4TouchableHandle touchPostStep = PostStep->GetTouchableHandle();
//To get the current volume:
G4VPhysicalVolume* volumePre = touchPreStep->GetVolume();
G4VPhysicalVolume* volumePost =touchPostStep->GetVolume();

//To get its name:
G4String namePre = volumePre->GetName();
G4String namePost;

if (!fScoringVolume) { 
    const DetectorConstruction* detectorConstruction
      = static_cast<const DetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume = detectorConstruction->GetScoringVolume();   
  }

// get volume of the current step
  G4LogicalVolume* volume
    = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetLogicalVolume();

  // check if we are in scoring volume
  if (volume != fScoringVolume) return;

  // collect energy deposited in this step
  G4double edepStep = step->GetTotalEnergyDeposit();
  fEventAction->AddEdep(edepStep);
  //G4cout << "\n\n !!!!!!!ENERGY DEPOSIT:  !!!!!!\n"<< G4BestUnit(edepStep,"Energy") <<"\n" << G4endl;
  


if(volumePost){
    namePost = volumePost->GetName(); 
}


G4int eventNum = G4RunManager::GetRunManager() -> GetCurrentEvent() -> GetEventID(); 
//G4double preeKin = step -> GetPreStepPoint() -> GetKineticEnergy()/eV;
G4double posteKin = step -> GetPostStepPoint() -> GetKineticEnergy()/eV;
//G4double differeceeEKIN = posteKin - preeKin;
G4double PosX = step->GetTrack()->GetPosition().x()/nanometer;
G4double PosY = step->GetTrack()->GetPosition().y()/nanometer;
G4double PosZ = step->GetTrack()->GetPosition().z()/nanometer;
//G4double globalTime = step -> GetTrack()->GetGlobalTime()/ns;
//G4String material= step -> GetTrack() -> GetMaterial() -> GetName();
G4String Step_Volume = step->GetTrack()->GetVolume()->GetName();
G4Track* theTrack = step->GetTrack();
//G4double VerX = theTrack->GetVertexPosition().x()/nanometer;
//G4double VerY = theTrack->GetVertexPosition().y()/nanometer;
//G4double VerZ = theTrack->GetVertexPosition().z()/nanometer;
//G4double phi = atan2(PosY, PosZ);
//G4double RADIUS = sqrt(PosX*PosX + PosY*PosY + PosZ*PosZ);
//G4double logOfrad = log10(RADIUS);

G4double edep = step->GetTotalEnergyDeposit()/eV;
//G4int nofHits = fHitsCollection->entries();

//G4cout << "X:   "<<PosX << "    Y:      "<<PosY << "    Z:      "<< PosZ << G4endl;


if((step->GetTrack()->GetDefinition()->GetParticleName() == "e-")  && (namePost != "physical-gold-sphere" && namePre == "physical-gold-sphere")  && (parentID == 0)){
//if(volume == "physical-gold-sphere" && step->GetTrack()->GetDefinition()->GetParticleName() == "e-"  && namePost =="physical-blue-sphere1"  && parentID == 0){
    	
	theTrack -> SetTrackStatus(fStopAndKill);
  //G4cout << "Particle KILLED" << G4endl; test is succesful
}

if((namePost == "physical-gold-sphere" && namePre == "physical-gold-sphere") && (step->GetTrack()->GetDefinition()->GetParticleName() == "e-") && (parentID > 0) && (posteKin < 10*eV)){
	
	theTrack -> SetTrackStatus(fStopAndKill);
  //G4cout << " Particle KIlled inside NP"<< G4endl;
}
/*
if((step->GetTrack()->GetDefinition()->GetParticleName() == "e-") && namePost == "physical-blue-sphere1" && (parentID > 0)){
		std::ofstream WriteDataIn("Secondaries.txt", std::ios::app);
		WriteDataIn		<<   edep       << '\t'
                    	<<   phi		<< '\t'
                    	<<   logOfrad   << '\t'
                    	<<   G4endl;
}*/
/*if((step->GetTrack()->GetDefinition()->GetParticleName() == "e-") && namePost == "physical-blue-sphere1" && (parentID > 0)){

    std::ofstream outbinary("SecondaryElectrons.dat",std::ios::app | std::ios::binary);
    //outbinary.write((char*)&parentID,sizeof(G4double));
    outbinary.write((char*)&preeKin,sizeof(G4double));
    outbinary.write((char*)&PosX,sizeof(G4double));
    outbinary.write((char*)&PosY,sizeof(G4double));
    outbinary.write((char*)&PosZ,sizeof(G4double));
    outbinary.write((char*)&edep,sizeof(G4double));
    outbinary.write((char*)&phi,sizeof(G4double));
    outbinary.write((char*)&logOfrad,sizeof(G4double));
    outbinary.close();

}*/

/*if((step->GetTrack()->GetDefinition()->GetParticleName() == "e-") || (step->GetTrack()->GetDefinition()->GetParticleName() == "gamma")){
		std::ofstream WriteDataIn("Secondaries.txt", std::ios::app);
		WriteDataIn		<<   edep       << '\t'
                    	<<   phi		<< '\t'
                    	<<   logOfrad   << '\t'
                    	<<   G4endl;
}*/

}
