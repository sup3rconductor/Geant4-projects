#include "DetPrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"

G4double x0, yy0, z0, teta, phi, ux, uy, uz, E0, x, y, z;
char fpartname[7];
G4int fEvent, fpartnum;
G4double fteta, fphi, fEkin;

extern G4double Z0const;
extern FILE* rdata;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetPrimaryGeneratorAction::DetPrimaryGeneratorAction()
	: G4VUserPrimaryGeneratorAction(),
	fParticleGun(0)
{
	G4int n_particle = 1;
	fParticleGun = new G4ParticleGun(n_particle);

	// default particle kinematic
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName;
	G4ParticleDefinition* particle = particleTable->FindParticle(particleName = "mu+");
	fParticleGun->SetParticleDefinition(particle);
	fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., -1.));
	fParticleGun->SetParticleEnergy(4. * GeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetPrimaryGeneratorAction::~DetPrimaryGeneratorAction()
{
	delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	//G4double x0 = 5*mm;//*(G4UniformRand()-0.5);
	//G4double y0 = 5*mm;//*(G4UniformRand()-0.5);
	//G4double z0 = 600*mm;
	x = -400 + 800 * G4UniformRand();//1/4
	y = 400 * G4UniformRand();//the same
	z = Z0const;
	//G4cout<<Z0const<<G4endl;

	//G4double costeta = pow(G4UniformRand(), 1/4.2);
	//G4cout<<costeta<<G4endl;

	//teta = acos(pow((1+G4UniformRand()*(pow(cos(15*twopi/360),4.2)-1)), 1/4.2));
	//phi = twopi*G4UniformRand();
	//fscanf(rdata, "%d\t%s\t%d\t%lf\t%lf\t%lf\n", &fEvent, &fpartname, &fpartnum, &fteta, &fphi, &fEkin);

	teta = fteta * pi / 180.0;
	phi = fphi * pi / 180.0;

	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4ParticleDefinition* particle = particleTable->FindParticle("mu+");
	fParticleGun->SetParticleDefinition(particle);

	/* x0 = x + 1500 * sin(teta) * cos(phi);
	yy0 = y + 1500 * sin(teta) * sin(phi);
	z0 = z + 1500 * cos(teta); */

	fParticleGun->SetParticlePosition(G4ThreeVector(0 * mm, 0 * mm, 30 * mm));

	ux = 0;
	uy = 0;
	uz = -1;
	/* ux = -sin(teta) * cos(phi);
	uy = -sin(teta) * sin(phi);
	uz = -cos(teta); */

	fParticleGun->SetParticleMomentumDirection(G4ThreeVector(ux, uy, uz));

	//E0 = 100./pow(G4UniformRand(),1/1.7);
	E0 = 1;

	fParticleGun->SetParticleEnergy(E0 * GeV);

	fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

