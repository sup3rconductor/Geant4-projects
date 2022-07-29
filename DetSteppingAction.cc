
#include "DetSteppingAction.hh"
#include "DetEventAction.hh"
#include "DetDetectorConstruction.hh"

#include "G4SystemOfUnits.hh"
#include "G4ParticleTypes.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

extern G4int Nevent;
extern G4double Z0const;
extern char fpartname[7];
G4int Nph;
G4double Tph[150000];
G4double Energy1, Eph[150000];
G4double scx, scy, scz, scux, scuy, scuz;

DetSteppingAction::DetSteppingAction(DetEventAction* eventAction)
: G4UserSteppingAction(), fEventAction(eventAction), fScoringVolume(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetSteppingAction::~DetSteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetSteppingAction::UserSteppingAction(const G4Step* step)
{
  G4Track* theTrack = step->GetTrack();
  G4String Name = theTrack->GetDefinition()->GetParticleName();
  G4ParticleDefinition* particleType = theTrack->GetDefinition();
  G4StepPoint* prePoint = step->GetPreStepPoint();
  G4String vname = prePoint->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetName();

  if(vname=="sc_l")
  {
    if(Name==fpartname)
    {
      //G4cout<<Name<<G4endl;
      G4ThreeVector worldPosition = theTrack->GetPosition();
      scx = worldPosition.x();
      scy = worldPosition.y();
      scz = worldPosition.z();
      G4ThreeVector direction = theTrack->GetMomentumDirection();
      scux = direction.x();
      scuy = direction.y();
      scuz = direction.z();
    }
  }

  if(vname=="alPMT_l") // PMT volume
  {
    if(particleType==G4OpticalPhoton::OpticalPhotonDefinition())
    {
    // collect the time
    //  Tph[Nph] = theTrack->GetGlobalTime()/ns;
    // collect the energy 
      Energy1 = theTrack->GetKineticEnergy()/eV;
      if(Energy1>=2.3 && Energy1<=3.2 && Nph<150000)
      {
        Eph[Nph] = Energy1;
        Tph[Nph] = theTrack->GetGlobalTime()/ns;
        Nph++;
      }
    // collect optical photons 
      theTrack->SetTrackStatus(fStopAndKill);
    }
    if(Name=="mu"||Name=="e"||Name=="proton")
    {
      theTrack->SetTrackStatus(fStopAndKill);
    }
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

