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
/// \file B4/B4c/src/EventAction.cc
/// \brief Implementation of the B4c::EventAction class

#include "EventAction.hh"

#include "CalorHit.hh"

#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4UnitsTable.hh"

#include <iomanip>

namespace B4c
{

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  CalorHitsCollection *EventAction::GetHitsCollection(G4int hcID, const G4Event *event) const
  {
    auto hitsCollection = static_cast<CalorHitsCollection *>(event->GetHCofThisEvent()->GetHC(hcID));

    if (!hitsCollection)
    {
      G4ExceptionDescription msg;
      msg << "Cannot access hitsCollection ID " << hcID;
      G4Exception("EventAction::GetHitsCollection()", "MyCode0003", FatalException, msg);
    }

    return hitsCollection;
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  void EventAction::PrintEventStatistics(G4double absoEdep, G4double absoTrackLength,
                                         G4double gapEdep, G4double gapTrackLength) const
  {
    // print event statistics
    G4cout << "   Absorber: total energy: " << std::setw(7) << G4BestUnit(absoEdep, "Energy")
           << "       total track length: " << std::setw(7) << G4BestUnit(absoTrackLength, "Length")
           << G4endl << "        Gap: total energy: " << std::setw(7) << G4BestUnit(gapEdep, "Energy")
           << "       total track length: " << std::setw(7) << G4BestUnit(gapTrackLength, "Length")
           << G4endl;
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  void EventAction::BeginOfEventAction(const G4Event * /*event*/) {}

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  void EventAction::EndOfEventAction(const G4Event *event)
  {
    // Get hits collections IDs (only once)
    if (fPWOhcID == -1)
    {
      fPWOhcID = G4SDManager::GetSDMpointer()->GetCollectionID("PWO_SD/PWO_HitsCollection");
    }

    // Get hits collections
    auto PWOhc = GetHitsCollection(fPWOhcID, event); // di tutti il calo, è un array con 19 indici uno per cristallo +1 per totali

    // Get hit with total values
    auto PWOHit = (*PWOhc)[PWOhc->entries() - 1];

    // Print per event (modulo n)
    //
    auto eventID = event->GetEventID();
    auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
    if ((printModulo > 0) && (eventID % printModulo == 0))
    {
      G4cout << "--> End of event: " << eventID << "\n"
             << G4endl;
    }

    // Fill histograms, ntuple
    //

    // get analysis manager
    auto analysisManager = G4AnalysisManager::Instance();

    // fill histograms
    analysisManager->FillH1(2, (*PWOhc)[18]->GetEdep());

    // fill ntuple
    G4double eneup = 0.;
    G4double enedown = 0.;
    for (G4int i = 0; i < 18; i++)
    {
      analysisManager->FillNtupleDColumn(i, (*PWOhc)[i]->GetEdep());
      if (i < 9)
        eneup += (*PWOhc)[i]->GetEdep();
      else
        enedown += (*PWOhc)[i]->GetEdep();
    }
    analysisManager->AddNtupleRow();
    analysisManager->FillH1(0, eneup);
    analysisManager->FillH1(1, enedown);
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

} // namespace B4c
