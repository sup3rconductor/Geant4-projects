#include "DetRunAction.hh"
#include "DetPrimaryGeneratorAction.hh"
#include "DetDetectorConstruction.hh"


#include "G4RunManager.hh"
#include "G4Run.hh"
//#include "G4AccumulableManager.hh"
//#include "G4LogicalVolumeStore.hh"
//#include "G4LogicalVolume.hh"
//#include "G4UnitsTable.hh"
//#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//extern G4int Nph;
//extern FILE *pt, *pos;
extern G4int Nevent;

DetRunAction::DetRunAction()
: G4UserRunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetRunAction::~DetRunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetRunAction::BeginOfRunAction(const G4Run*)
{ 
  // inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetRunAction::EndOfRunAction(const G4Run* run)
{
  Nevent = run->GetNumberOfEvent();
  if (Nevent==0) return;
  
  //fclose(pt);
  //fclose(pos);
}
