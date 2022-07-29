
#include "DetActionInitialization.hh"
#include "DetPrimaryGeneratorAction.hh"
#include "DetRunAction.hh"
#include "DetEventAction.hh"
#include "DetSteppingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetActionInitialization::DetActionInitialization()
 : G4VUserActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetActionInitialization::~DetActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetActionInitialization::BuildForMaster() const
{
  DetRunAction* runAction = new DetRunAction;
  SetUserAction(runAction);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetActionInitialization::Build() const
{
  SetUserAction(new DetPrimaryGeneratorAction);

  DetRunAction* runAction = new DetRunAction;
  SetUserAction(runAction);
  
  DetEventAction* eventAction = new DetEventAction(runAction);
  SetUserAction(eventAction);
  
  SetUserAction(new DetSteppingAction(eventAction));
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
