
#include "DetEventAction.hh"
#include "DetRunAction.hh"
#include "DetSteppingAction.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
extern G4int Nph;
extern G4int Nevent;
extern FILE *pt, *pos;

extern G4double Tph[150000];
extern G4double Eph[150000], E0;
extern G4double x0, yy0, z0, teta, phi, x, y, z, ux, uy, uz;
extern G4double scx, scy, scz, scux, scuy, scuz;
extern char fpartname[7];

G4int j;

DetEventAction::DetEventAction(DetRunAction* runAction)
: G4UserEventAction(), fRunAction(runAction)
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetEventAction::~DetEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetEventAction::BeginOfEventAction(const G4Event*)
{    
 Nph = 0;
 for(j=0; j<150000; j++)
 {
   Tph[j] = 0;
   Eph[j] = 0;
 }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetEventAction::EndOfEventAction(const G4Event*)
{   
 if(pt)
 {
   //G4cout<<"Nevent="<<Nevent<<"\tE0="<<E0<<"\tx0="<<x0<<"\ty0="<<yy0<<"\tz0="<<z0<<"\tteta="<<teta<<"\tphi="<<phi<<"\tNph="<<Nph<<G4endl;
   //G4cout<<"RandX="<<x<<"\tRandY="<<y<<"\tRandZ="<<z<<G4endl;
   //G4cout<<"Particle type "<<fpartname<<G4endl;
   //G4cout<<"Event number="<<Nevent<<"\tParticle type: "<<fpartname<<"\tTeta="<<teta<<G4endl;
   fprintf(pt, "%d\t%s\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", Nevent, fpartname, E0, x0, yy0, z0, teta/pi*180.0, phi/pi*180.0, ux, uy, uz, Nph, scx, scy, scz, scux, scuy, scuz);
   fprintf(pos, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", scx, scy, scz, scux, scuy, scuz);
   for(j=0; j<Nph; j++)
   {
     fprintf(pt, "%lf\t%lf\n", Eph[j], Tph[j]);
   }
   fprintf(pt, "\n");
 }
 
 Nevent++;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
