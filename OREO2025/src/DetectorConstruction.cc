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
/// \file B4/B4c/src/DetectorConstruction.cc
/// \brief Implementation of the B4c::DetectorConstruction class

#include "DetectorConstruction.hh"

#include "CalorimeterSD.hh"

#include "G4AutoDelete.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4PhysicalConstants.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"

namespace B4c
{

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  G4ThreadLocal G4GlobalMagFieldMessenger *DetectorConstruction::fMagFieldMessenger = nullptr;

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  G4VPhysicalVolume *DetectorConstruction::Construct()
  {
    // Define materials
    DefineMaterials();

    // Define volumes
    return DefineVolumes();
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  void DetectorConstruction::DefineMaterials()
  {
    // Lead material defined using NIST Manager
    auto nistManager = G4NistManager::Instance();
    nistManager->FindOrBuildMaterial("G4_Pb");
    nistManager->FindOrBuildMaterial("G4_PbWO4");
    nistManager->FindOrBuildMaterial("G4_GLASS_LEAD");
    nistManager->FindOrBuildMaterial("G4_AIR");
    nistManager->FindOrBuildMaterial("G4_Si");

    // Liquid argon material
    G4double a; // mass of a mole;
    G4double z; // z=mean number of protons;
    G4double density;
    new G4Material("liquidArgon", z = 18., a = 39.95 * g / mole, density = 1.390 * g / cm3);
    // The argon by NIST Manager is a gas with a different density

    // Vacuum
    new G4Material("Galactic", z = 1., a = 1.01 * g / mole, density = universe_mean_density,
                   kStateGas, 2.73 * kelvin, 3.e-18 * pascal);

    // Print materials
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  G4VPhysicalVolume *DetectorConstruction::DefineVolumes()
  {
    // Geometry parameters
    G4int nLayers = 2;
    G4int nCrystals = 9;
    G4ThreeVector crystalSizeUp(2.5 * cm, 2.5 * cm, 4.5 * cm);
    G4ThreeVector crystalSizeDown(2.5 * cm, 2.5 * cm, 10.0 * cm);
    G4double gapLayers = 0.01 * cm;
    G4ThreeVector teleSize(9.3 * cm, 9.3 * cm, 0.041 * cm);

    auto OREOsize = crystalSizeUp + crystalSizeDown + G4ThreeVector(0, 0, gapLayers);
    auto worldSizeXY = 5 * m;
    auto worldSizeZ = 10 * m; // 4.332 * cm + 244.5 * cm + 4.333 * cm + 30.3 * cm + OREOsize.z() + 7.6 * cm + 30.0 * cm;

    // Get materials
    auto defaultMaterial = G4Material::GetMaterial("Galactic");
    auto absorberMaterial = G4Material::GetMaterial("G4_PbWO4");
    auto gapMaterial = G4Material::GetMaterial("G4_AIR");
    auto teleMaterial = G4Material::GetMaterial("G4_Si");

    if (!defaultMaterial || !absorberMaterial || !gapMaterial || !teleMaterial)
    {
      G4ExceptionDescription msg;
      msg << "Cannot retrieve materials already defined.";
      G4Exception("DetectorConstruction::DefineVolumes()", "MyCode0001", FatalException, msg);
    }

    //
    // World
    //
    auto worldS = new G4Box("World",                                           // its name
                            worldSizeXY / 2, worldSizeXY / 2, worldSizeZ / 2); // its size

    auto worldLV = new G4LogicalVolume(worldS,          // its solid
                                       defaultMaterial, // its material
                                       "World");        // its name

    auto worldPV = new G4PVPlacement(nullptr,         // no rotation
                                     G4ThreeVector(), // at (0,0,0)
                                     worldLV,         // its logical volume
                                     "World",         // its name
                                     nullptr,         // its mother  volume
                                     false,           // no boolean operation
                                     100,               // copy number
                                     fCheckOverlaps); // checking overlaps


    //tele
    auto teleS = new G4Box("tele", teleSize.x() / 2, teleSize.y() / 2, teleSize.z() / 2);

    auto teleLV = new G4LogicalVolume(teleS,         // its solid
                                      teleMaterial,   // its material
                                      "tele");        // its name

    auto telePV1 = new G4PVPlacement(nullptr,         // no rotation
                                    G4ThreeVector(0, 0 , 10. * cm),
                                    teleLV,         // its logical volume
                                    "tele1",         // its name
                                    worldLV,         // its mother  volume
                                    false,           // no boolean operation
                                    200,               // copy number
                                    fCheckOverlaps); // checking overlaps

    auto telePV2 = new G4PVPlacement(nullptr,         // no rotation
                                    G4ThreeVector(0, 0 , 254.5 * cm),
                                    teleLV,         // its logical volume
                                    "tele2",         // its name
                                    worldLV,         // its mother  volume
                                    false,           // no boolean operation
                                    201,               // copy number
                                    fCheckOverlaps); // checking overlaps
                                        
    //
    // Calorimeter
    //

    auto PWO_crystal_up = new G4Box("PWO_crystal", // its name
                                    crystalSizeUp.x() / 2, crystalSizeUp.y() / 2,
                                    crystalSizeUp.z() / 2);
    auto PWO_crystal_down = new G4Box("PWO_crystal", // its name
                                      crystalSizeDown.x() / 2, crystalSizeDown.y() / 2,
                                      crystalSizeDown.z() / 2);
    auto PWO_LV_up = new G4LogicalVolume(PWO_crystal_up,     // its solid
                                         absorberMaterial,   // its material
                                         "PWO_LV_up");       // its name
    auto PWO_LV_down = new G4LogicalVolume(PWO_crystal_down, // its solid
                                           absorberMaterial, // its material
                                           "PWO_LV_down");   // its name

    // layer up
    G4double zOreo = 3 * m;
    const G4int nRows = 3, nCols = 3;
    const G4double pitchX = crystalSizeUp.x();
    const G4double pitchY = crystalSizeUp.y();

    for (G4int i = 0; i < nRows * nCols; ++i)
    {
      G4int row = i / nCols; // 0..2
      G4int col = i % nCols; // 0..2

      // centro la matrice: col,row ∈ {0,1,2} → { -1, 0, +1 }
      G4double x = (col - (nCols - 1) / 2.0) * pitchX;
      G4double y = (row - (nRows - 1) / 2.0) * pitchY;
      G4double z1 = zOreo + 0.5 * crystalSizeUp.z();                                   // come nel tuo codice
      G4double z2 = zOreo + crystalSizeUp.z() + 0.5 * crystalSizeDown.z() + gapLayers; // come nel tuo codice

      new G4PVPlacement(nullptr,
                        G4ThreeVector(x, y, z1),
                        PWO_LV_up,
                        "PWO",
                        worldLV,
                        false, i, fCheckOverlaps);

      new G4PVPlacement(nullptr,
                        G4ThreeVector(x, y, z2),
                        PWO_LV_down,
                        "PWO",
                        worldLV,
                        false, i+9, fCheckOverlaps);
    }

    // layer down

    //
    // Visualization attributes
    //
    worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());
    PWO_LV_up->SetVisAttributes(G4VisAttributes(G4Colour::Blue()));
    PWO_LV_down->SetVisAttributes(G4VisAttributes(G4Colour::Red()));
    teleLV->SetVisAttributes(G4VisAttributes(G4Colour::Green()));

    //
    // Always return the physical World
    //
    return worldPV;
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  void DetectorConstruction::ConstructSDandField()
  {
    G4SDManager::GetSDMpointer()->SetVerboseLevel(1);
    G4int ncells = 18;

    // Sensitive detectors

    auto pwoSD = new CalorimeterSD("PWO_SD", "PWO_HitsCollection", ncells);
    G4SDManager::GetSDMpointer()->AddNewDetector(pwoSD);
    SetSensitiveDetector("PWO_LV_up", pwoSD);
    SetSensitiveDetector("PWO_LV_down", pwoSD);
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

} // namespace B4c
