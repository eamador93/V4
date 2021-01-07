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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DetectorConstruction::DetectorConstruction() 
: G4VUserDetectorConstruction(), 
  //fpWaterMaterial(0), 
  //fUserLimits(0),
  fScoringVolume(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DetectorConstruction::~DetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* DetectorConstruction::Construct()

{
  DefineMaterials();
  return ConstructDetector();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DetectorConstruction::DefineMaterials()
{
  /*
  // Water is defined from NIST material database
  G4NistManager * man = G4NistManager::Instance();

  G4Material * H2O = man->FindOrBuildMaterial("G4_WATER");
  G4Material * Gold = man->FindOrBuildMaterial("G4_Au");
  //G4Material * Air = man->FindOrBuildMaterial("G4_AIR");
  
   If one wishes to test other density value for water material,
   one should use instead:
   G4Material * H2O = man->BuildMaterialWithNewDensity("G4_WATER_MODIFIED",
   "G4_WATER",1.100*g/cm3);

   Note: any string for "G4_WATER_MODIFIED" parameter is accepted
   and "G4_WATER" parameter should not be changed
   Both materials are created and can be selected from dna.mac
   
  fpWaterMaterial = H2O;
  fpNanoMaterial = Gold;
  //G4cout << "-> Density of water material (g/cm3)="
  // << fpWaterMaterial->GetDensity()/(g/cm/cm/cm) << G4endl;

  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4VPhysicalVolume* DetectorConstruction::ConstructDetector()
{

  G4NistManager * man = G4NistManager::Instance();

  G4Material * H2O = man->FindOrBuildMaterial("G4_WATER");
  G4Material * Gold = man->FindOrBuildMaterial("G4_Au");


  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

  // WORLD VOLUME DIMENSIONS

  G4double worldSizeX = 2 * m;
  G4double worldSizeY = worldSizeX;
  G4double worldSizeZ = worldSizeX;

  //Shapes
  
  G4Box* solidWorld = new G4Box("World", worldSizeX / 2, worldSizeY / 2, worldSizeZ / 2); 

  G4Orb* gold_solid = new G4Orb("gold-sphere", 50 *nm);

  G4Sphere* solidBp1 = new G4Sphere("blue-sphere1", 50*nm,  72.6466*nm, 0, 2*CLHEP::pi, 0, CLHEP::pi); 
  
  //72.6466nm original outer radius
  //G4Sphere* solidBp1 = new G4Sphere("blue-sphere1", 50*nm,  2*cm, 0, 2*CLHEP::pi, 0, CLHEP::pi);
 
  //make smaller outer radius 
  //check exmaple  B2b 
  
  //Logical Volumes

  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, H2O, "World"); 

  G4LogicalVolume* gold_logic = new G4LogicalVolume(gold_solid, Gold,"logic-gold-sphere");
  
  G4LogicalVolume* logicBp1 = new G4LogicalVolume(solidBp1, H2O,"logic-blue-sphere1");
   
   //gold_logic->SetUserLimits(fUserLimits);

  //Physical Volumes

  G4PVPlacement* physiWorld = new G4PVPlacement(0, G4ThreeVector(),"World", logicWorld, 0, false, 0);

  G4PVPlacement* gold_physical = new G4PVPlacement(0, G4ThreeVector(0,0,0), gold_logic, "physical-gold-sphere", logicWorld, false, 0);

  G4PVPlacement* physiBp1 = new G4PVPlacement(0, G4ThreeVector(0,0,0), logicBp1, "physical-blue-sphere1", logicWorld, false, 0);

  // Visualization attributes
  G4VisAttributes* worldVisAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0)); //
  worldVisAtt->SetVisibility(true);
  logicWorld->SetVisAttributes(worldVisAtt);
 
  auto GoldAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0)); //gold
  GoldAtt->SetVisibility(true);
  gold_logic->SetVisAttributes(GoldAtt);

  auto WaterAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0)); //blue
  WaterAtt->SetVisibility(true);
  logicBp1->SetVisAttributes(WaterAtt);

  fScoringVolume = logicBp1;
  // 
  // Shows how to introduce a 20 eV tracking cut
  //
  //logicWorld->SetUserLimits(new G4UserLimits(DBL_MAX,DBL_MAX,DBL_MAX,20*eV));

  return physiWorld;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*
void DetectorConstruction::SetMaterial(G4String materialChoice)
{
  // Search the material by its name   
  G4Material* pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(
      materialChoice);

  if (pttoMaterial)
  {
    fpWaterMaterial = pttoMaterial;
    G4LogicalVolume* logicWorld =
        G4LogicalVolumeStore::GetInstance()->GetVolume("World");
    logicWorld->SetMaterial(fpWaterMaterial);
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
  }
}*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructDetector());
}


void DetectorConstruction::SetMaxStepSize(G4double val)
{
  // set the maximum allowed step size
  //
  if (val <= DBL_MIN)
    { G4cout << "\n --->warning from SetMaxStepSize: maxStep "
             << val  << " out of range. Command refused" << G4endl;
      return;
    }
  //fUserLimits->SetMaxAllowedStep(1*nm);
}
