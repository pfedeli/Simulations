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
/// \file B4/B4c/src/RunAction.cc
/// \brief Implementation of the B4::RunAction class

#include "RunAction.hh"

#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "globals.hh"

namespace B4
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
{
  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);

  // Create analysis manager
  // The choice of the output format is done via the specified
  // file extension.
  auto analysisManager = G4AnalysisManager::Instance();

  // Create directories
  // analysisManager->SetHistoDirectoryName("histograms");
  // analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetNtupleMerging(true);
  // Note: merging ntuples is available only with Root output

  // Book histograms, ntuple
  //

  // Creating histograms
  analysisManager->CreateH1("Eup", "Edep in layer up", 200, 0., 200 * GeV);
  analysisManager->CreateH1("Edown", "Edep in layer down", 200, 0., 200 * GeV);
  analysisManager->CreateH1("Etot", "Edep in total", 200, 0., 200 * GeV);

  // Creating ntuple
  //
  analysisManager->CreateNtuple("OREO", "Edep");
  analysisManager->CreateNtupleDColumn("Cry0");
  analysisManager->CreateNtupleDColumn("Cry1");
  analysisManager->CreateNtupleDColumn("Cry2");
  analysisManager->CreateNtupleDColumn("Cry3");
  analysisManager->CreateNtupleDColumn("Cry4");
  analysisManager->CreateNtupleDColumn("Cry5");
  analysisManager->CreateNtupleDColumn("Cry6");
  analysisManager->CreateNtupleDColumn("Cry7");
  analysisManager->CreateNtupleDColumn("Cry8");
  analysisManager->CreateNtupleDColumn("Cry9");
  analysisManager->CreateNtupleDColumn("Cry10");
  analysisManager->CreateNtupleDColumn("Cry11");
  analysisManager->CreateNtupleDColumn("Cry12");
  analysisManager->CreateNtupleDColumn("Cry13");
  analysisManager->CreateNtupleDColumn("Cry14");
  analysisManager->CreateNtupleDColumn("Cry15");
  analysisManager->CreateNtupleDColumn("Cry16");
  analysisManager->CreateNtupleDColumn("Cry17");
  analysisManager->FinishNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* /*run*/)
{
  // inform the runManager to save random number seed
  // G4RunManager::GetRunManager()->SetRandomNumberStore(true);

  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // Open an output file
  //
  //G4String fileName = "todelete.root";
  // Other supported output types:
  // G4String fileName = "B4.csv";
  // G4String fileName = "B4.hdf5";
  // G4String fileName = "B4.xml";
  analysisManager->OpenFile(); //  analysisManager->OpenFile(fileName); se da c++
  G4cout << "Using " << analysisManager->GetType() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* /*run*/)
{
  // print histogram statistics
  //
  auto analysisManager = G4AnalysisManager::Instance();
  if (analysisManager->GetH1(1)) {
    G4cout << G4endl << " ----> print histograms statistic ";
    if (isMaster) {
      G4cout << "for the entire run " << G4endl << G4endl;
    }
    else {
      G4cout << "for the local thread " << G4endl << G4endl;
    }

    G4cout << " Eup : mean = " << G4BestUnit(analysisManager->GetH1(0)->mean(), "Energy")
           << " rms = " << G4BestUnit(analysisManager->GetH1(0)->rms(), "Energy") << G4endl;

    G4cout << " Edown : mean = " << G4BestUnit(analysisManager->GetH1(1)->mean(), "Energy")
           << " rms = " << G4BestUnit(analysisManager->GetH1(1)->rms(), "Energy") << G4endl;

    G4cout << " Etot : mean = " << G4BestUnit(analysisManager->GetH1(2)->mean(), "Energy")
           << " rms = " << G4BestUnit(analysisManager->GetH1(2)->rms(), "Energy") << G4endl;

  }

  // save histograms & ntuple
  //
  analysisManager->Write();
  analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace B4
