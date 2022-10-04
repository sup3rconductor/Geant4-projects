#include "DetDetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SubtractionSolid.hh"
#include "G4RotationMatrix.hh"

#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpBoundaryProcess.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double Z0const, X0const, Y0const;

DetDetectorConstruction::DetDetectorConstruction()
	: G4VUserDetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetDetectorConstruction::~DetDetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetDetectorConstruction::Construct()
{
	/*	MATERIALS	*/

	G4double a, z;  //Atomic mass, atomic number
	G4double density, fractionmass;
	G4int ncomponents, nelements;


	//Chemical elements
	G4Element* elH = new G4Element("Hydrogen", "H", z = 1., a = 1.01 * g / mole);
	G4Element* elC = new G4Element("Carbon", "C", z = 6., a = 12.01 * g / mole);
	G4Element* elN = new G4Element("Nitrogen", "N", z = 7., a = 14.01 * g / mole);
	G4Element* elO = new G4Element("Oxygen", "O", z = 8., a = 16.00 * g / mole);
	G4Element* elSi = new G4Element("Silicium", "Si", z = 14., a = 28.09 * g / mole);
	G4Element* elAl = new G4Element("Aluminium", "Al", z = 13., a = 26.98 * g / mole);
	G4Element* elB = new G4Element("Boron", "B", z = 5., a = 10.812 * g / mole);
	G4Element* elFe = new G4Element("Ferrum", "Fe", z = 26., a = 55.85 * g / mole);


	//Air
	G4Material* Air = new G4Material("MAir", density = 1.290 * mg / cm3, ncomponents = 2);
	Air->AddElement(elN, fractionmass = 0.8);
	Air->AddElement(elO, fractionmass = 0.2);

	//Aluminium
	G4Material* AlMaterial = new G4Material("MAluminium", z = 13., a = 26.98 * g / mole, density = 2.8 * g / cm3);

	//Ferrum
	G4Material* FeMaterial = new G4Material("MFerrum", z = 26., a = 55.85 * g / mole, density = 7.9 * g / cm3);

	//Scintillator material
	G4Material* Scint = new G4Material("MScint", density = 1.032 * g / cm3, ncomponents = 2);
	Scint->AddElement(elC, fractionmass = 0.476);
	Scint->AddElement(elH, fractionmass = 0.524);

	//Outer cover of fiberglass
	G4Material* FP = new G4Material("MFP", density = 1.43 * g / cm3, ncomponents = 3);
	FP->AddElement(elC, nelements = 5);
	FP->AddElement(elH, nelements = 8);
	FP->AddElement(elO, nelements = 2);

	//Inner cover of fiberglass 
	G4Material* PMMA = new G4Material("MPMMA", density = 1.19 * g / cm3, ncomponents = 3);
	PMMA->AddElement(elC, nelements = 5);
	PMMA->AddElement(elH, nelements = 8);
	PMMA->AddElement(elO, nelements = 2);

	//Core of fiberglass
	G4Material* PS = new G4Material("MPS", density = 1.05 * g / cm3, ncomponents = 2);
	PS->AddElement(elC, nelements = 8);
	PS->AddElement(elH, nelements = 8);

	//Photocatode (borosilicate glass)
	G4Material* SiO2 = new G4Material("MSiO2", density = 2.1 * g / cm3, ncomponents = 2);
	SiO2->AddElement(elSi, nelements = 1);
	SiO2->AddElement(elO, nelements = 2);
	G4Material* B2O3 = new G4Material("MB2O3", density = 2.26 * g / cm3, ncomponents = 2);
	B2O3->AddElement(elB, nelements = 2);
	B2O3->AddElement(elO, nelements = 3);
	G4Material* Al2O3 = new G4Material("MAl2O3", density = 3.99 * g / cm3, ncomponents = 2);
	Al2O3->AddElement(elAl, nelements = 2);
	Al2O3->AddElement(elO, nelements = 3);
	G4Material* Na2O = new G4Material("MNa2O", density = 2.27 * g / cm3, ncomponents = 2);
	Na2O->AddElement(elAl, nelements = 2);
	Na2O->AddElement(elO, nelements = 1);

	G4Material* PhotCat = new G4Material("MPhotCat", density = 2.5 * g / cm3, ncomponents = 4);
	PhotCat->AddMaterial(SiO2, fractionmass = 80. * perCent);
	PhotCat->AddMaterial(B2O3, fractionmass = 14. * perCent);
	PhotCat->AddMaterial(Al2O3, fractionmass = 4. * perCent);
	PhotCat->AddMaterial(Na2O, fractionmass = 2. * perCent);




	/*	DETECTOR	*/


	G4bool checkOverlaps = true;

	//World
	G4double world_sizeX = 5 * m;
	G4double world_sizeY = 5 * m;
	G4double world_sizeZ = 5 * m;

	G4Box* solidWorld = new G4Box("World_s", 0.5 * world_sizeX, 0.5 * world_sizeY, 0.5 * world_sizeZ);
	G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, Air, "World_l");
	G4VPhysicalVolume* physWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "World", 0, false, 0, checkOverlaps);

	//Scintillator
	G4double ScrLength = 200 * mm;
	G4double ScrWidth = 200 * mm;
	G4double ScrHeight = 5 * mm;

	//Holes
	G4double HoleLength = ScrLength;
	G4double HoleWidth = 1. * mm;
	G4double HoleHeight = 2. * mm;
	G4RotationMatrix* HoleRot = new G4RotationMatrix;

	//Optical fiber
	G4double OptRad = 0.5 * mm;
	G4double OptHeight = 100 * mm;
	G4double CovThickness = 0.02 * mm;
	G4RotationMatrix* OptRot = new G4RotationMatrix;
	OptRot->rotateY(90. * deg);

	//Photomultiplier
	G4double GlassRad = 1.5 * mm;
	G4double GlassHeight = 1 * mm;

	G4double PhotRad = GlassRad;
	G4double PhotHeight = 0.1 * mm;

	G4double BodyOuterRad = GlassRad + 0.1 * mm;
	G4double BodyHeight = GlassHeight + PhotHeight;

	//Detector steel shell and air hollow
	G4double GapH = 0.1 * mm;		//Gap between plates on same level
	G4double GapV = 0.1 * mm;		//Gap between levels 
	G4double GapFP = 0.1 * mm;		//Gap between frame and plate

	G4double HollowLength = 5 * ScrLength + 6 * GapH + BodyHeight;
	G4double HollowWidth = 5 * ScrWidth + 6 * GapH;
	G4double HollowHeight = 8 * ScrHeight + 9 * GapV;

	G4double ShellThickness = 2 * mm;

	G4double ShellLength = HollowLength + 2 * ShellThickness;
	G4double ShellWidth = HollowWidth + 2 * ShellThickness;
	G4double ShellHeight = HollowHeight + 2 * ShellThickness;

	//Start position
	G4double X0 = 0.5 * (HollowLength - ScrLength) - GapFP - BodyHeight;
	G4double Y0 = -0.5 * (HollowWidth - ScrWidth) + GapFP;
	G4double Z0 = 0.5 * (HollowHeight - ScrHeight) - GapFP;
	G4double StartPosHole = 47. * mm;
	G4double StartPosOpt = StartPosHole;
	G4double distance = 35. * mm; //distance between holes

	const G4int NLvls = 8, NRows = 5, NCols = 5;		//Scintillation plates
	const G4int NOpt = 20;								//Num of optical fiber lines in one row
	const G4int NPhot = 20;								//Num of photomultipliers in one row
	const G4int NHoles = 4;								//Num of holes in one scintillation plate
	const G4int Nhls = 20;

	G4double Glass_zpos = 0.5 * BodyHeight - 0.5 * GlassHeight;
	G4double Phot_zpos = -0.5 * BodyHeight + 0.5 * PhotHeight;

	//Positions of plate, optical fiber and photomultier
	G4double XPlate = X0, YPlate = Y0, ZPlate = Z0;
	G4double XPhot = X0, YPhot = Y0, ZPhot = Z0;

	G4Box* solidScintplate[NLvls][NCols][NRows] = { NULL }, * solidHole[NLvls][NCols][Nhls] = { NULL };
	G4Tubs* solidCore[NLvls][NCols][NOpt] = { NULL }, * solidInCov[NLvls][NCols][NOpt] = { NULL }, * solidOutCov[NLvls][NCols][NOpt] = { NULL },
		* solidBody[NLvls][NPhot] = { NULL }, * solidGlass[NLvls][NPhot] = { NULL }, * solidPhot[NLvls][NPhot] = { NULL };
	G4LogicalVolume* logicScintplate[NLvls][NCols][NRows] = { NULL }, * logicHole[NLvls][NCols][Nhls] = { NULL }, * logicCore[NLvls][NCols][NOpt] = { NULL }, * logicInCov[NLvls][NCols][NOpt] = { NULL },
		* logicOutCov[NLvls][NCols][NOpt] = { NULL }, * logicBody[NLvls][NPhot] = { NULL }, * logicGlass[NLvls][NPhot] = { NULL }, * logicPhot[NLvls][NPhot] = { NULL };
	G4VPhysicalVolume* physScintplate[NLvls][NCols][NRows] = { NULL }, *physHole[NLvls][NCols][Nhls] = { NULL }, *physCore[NLvls][NCols][NOpt] = { NULL }, *physInCov[NLvls][NCols][NOpt] = { NULL },
		*physOutCov[NLvls][NCols][NOpt] = { NULL }, *physBody[NLvls][NPhot] = { NULL }, *physGlass[NLvls][NPhot] = { NULL }, *physPhot[NLvls][NPhot] = { NULL };

	//Steel shell
	G4Box* solidShell = new G4Box("shell_s", 0.5 * ShellLength, 0.5 * ShellWidth, 0.5 * ShellHeight);
	G4LogicalVolume* logicShell = new G4LogicalVolume(solidShell, FeMaterial, "shell_l");
	G4VPhysicalVolume* physShell = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicShell, "Shell", logicWorld, false, 0, checkOverlaps);

	//Air hollow inside shell
	G4Box* solidHollow = new G4Box("hollow_s", 0.5 * HollowLength, 0.5 * HollowWidth, 0.5 * HollowHeight);
	G4LogicalVolume* logicHollow = new G4LogicalVolume(solidHollow, Air, "hollow_l");
	G4VPhysicalVolume* physHollow = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicHollow, "Hollow", logicShell, false, 0, checkOverlaps);

	//Variables for creating copies
	G4int row, column, level;
	G4int opt;
	G4int PlateNCopy = 0, HoleNCopy = 0, OptNCopy = 0, PhotNCopy = 0;
	G4int HoleCount = 0, OptCount = 0, PhotCount = 0;


	/* SOLID AND LOGICAL VOLUMES OF DETECTOR ELEMENTS */

	//Hole position
	G4double H_x[NHoles], H_y[NHoles], H_z[NHoles];
	H_x[0] = 0 * mm;
	H_y[0] = 0.5 * ScrWidth - StartPosHole;
	H_z[0] = 0.5 * ScrHeight - 0.5 * HoleHeight;

	G4int h;
	for (h = 1; h < NHoles; h++)
	{
		H_x[h] = H_x[h - 1];
		H_y[h] = H_y[h - 1] - distance;
		H_z[h] = H_z[h - 1];
	}


	/* physScintplate[0][0][0] = new G4PVPlacement(0, G4ThreeVector(2 * ScrLength, -2 * ScrWidth, 0), logicScintplate, "scintplate", logicHollow, true, 0); */

	//Optical fiber position
	G4double XOpt[NHoles], YOpt[NHoles], ZOpt[NHoles];
	for (h = 0; h < NHoles; h++)
	{
		XOpt[h] = H_x[h];
		YOpt[h] = H_y[h];
		ZOpt[h] = H_z[h] - OptRad;
	}

	//Constructing detector
	for (level = 0; level < NLvls; level++)
	{
		for (column = 0; column < NCols; column++)
		{
			for (row = 0; row < NRows; row++)
			{
				solidScintplate[level][column][row] = new G4Box("scintplate_s", 0.5 * ScrLength, 0.5 * ScrWidth, 0.5 * ScrHeight);
				logicScintplate[level][column][row] = new G4LogicalVolume(solidScintplate[level][column][row], Scint, "scintplate_l");
				physScintplate[level][column][row] = new G4PVPlacement(0, G4ThreeVector(XPlate, YPlate, ZPlate), logicScintplate[level][column][row], "scintplate", logicHollow, false, PlateNCopy, checkOverlaps);

				for (opt = 0; opt < 4; opt++)
				{
					//Hole for optical fiber
					solidHole[level][column][HoleCount] = new G4Box("hole_s", 0.5 * HoleLength, 0.5 * HoleWidth, 0.5 * HoleHeight);
					logicHole[level][column][HoleCount] = new G4LogicalVolume(solidHole[level][column][HoleCount], Air, "hole_l");
					physHole[level][column][HoleCount] = new G4PVPlacement(0, G4ThreeVector(H_x[opt], H_y[opt], H_z[opt]), logicHole[level][column][HoleCount], "scintplate", logicScintplate[level][column][row], false, HoleNCopy, checkOverlaps);

					//Outer cover
					solidOutCov[level][column][OptCount] = new G4Tubs("OutCov_s", 0, OptRad, OptHeight, 0. * deg, 360. * deg);
					logicOutCov[level][column][OptCount] = new G4LogicalVolume(solidOutCov[level][column][OptCount], FP, "OutCov_l");
					physOutCov[level][column][OptCount] = new G4PVPlacement(OptRot, G4ThreeVector(XOpt[opt], YOpt[opt], ZOpt[opt]), logicOutCov[level][column][OptCount], "OUTER COVER", logicScintplate[level][column][row], false, OptNCopy, checkOverlaps);

					//Inner cover
					solidInCov[level][column][OptCount] = new G4Tubs("InCov_s", 0, OptRad - CovThickness, OptHeight, 0. * deg, 360. * deg);
					logicInCov[level][column][OptCount] = new G4LogicalVolume(solidInCov[level][column][OptCount], PMMA, "InCov_l");
					physInCov[level][column][OptCount] = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicInCov[level][column][OptCount], "INNER COVER", logicOutCov[level][column][OptCount], false, OptNCopy, checkOverlaps);

					//Core
					solidCore[level][column][OptCount] = new G4Tubs("core_s", 0, OptRad - 2 * CovThickness, OptHeight, 0. * deg, 360. * deg);
					logicCore[level][column][OptCount] = new G4LogicalVolume(solidCore[level][column][OptCount], PS, "core_l");
					physCore[level][column][OptCount] = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicCore[level][column][OptCount], "CORE", logicInCov[level][column][OptCount], false, OptNCopy, checkOverlaps);

					HoleCount++;
					HoleNCopy++;
					OptNCopy++;
					OptCount++;
				}

				YPlate += (ScrWidth + GapH);
				PlateNCopy++;
				HoleCount = 0;
				OptCount = 0;
			}
			XPlate -= (ScrLength + GapH);
			YPlate = Y0;
		}
		ZPlate -= (ScrHeight + GapV);
		XPlate = X0;
		YPlate = Y0;
	}




	/*	OPTICAL PROPERTIES	*/


	//Scintillator optical properties
	const G4int nEntries = 60;
	G4double PhotonEnergy[nEntries] = { 2.3, 2.31525, 2.33051, 2.34576, 2.36102, 2.37627, 2.39153, 2.40678, 2.42203, 2.43729, 2.45254, 2.4678, 2.48305, 2.49831, 2.51356,
		 2.52881, 2.54407, 2.55932, 2.57458, 2.58983, 2.60508, 2.62034, 2.63559, 2.65085, 2.6661, 2.68136, 2.69661, 2.71186, 2.72712, 2.74237,
		 2.75763, 2.77288, 2.78814, 2.80339, 2.81864, 2.8339, 2.84915, 2.86441, 2.87966, 2.89492, 2.91017, 2.92542, 2.94068, 2.95593, 2.97119,
		 2.98644, 3.00169, 3.01695, 3.0322, 3.04746, 3.06271, 3.07797, 3.09322, 3.10847, 3.12373, 3.13898, 3.15424, 3.16949, 3.18475, 3.2 };
	G4double RefractiveScin[nEntries];
	G4double AbsLengthScin[nEntries];
	G4double SpIzlStr[nEntries] = { 0, 0, 0.04304, 0.09311, 0.14318, 0.19325, 0.24331, 0.29338, 0.34345, 0.39352, 0.44359, 0.49365, 0.54372, 0.59379, 0.65703,
		 0.72516, 0.7829, 0.85487, 0.93619, 1.0156, 1.10002, 1.19322, 1.29936, 1.41172, 1.53233, 1.65876, 1.79893, 1.98186, 2.18771, 2.4366,
		 2.78324, 3.0698, 3.27276, 3.39218, 3.46918, 3.4941, 3.52619, 3.60856, 3.88683, 4.28688, 4.71702, 4.93565, 4.80817, 4.56821, 4.23367,
		 3.56117, 2.30136, 1.47323, 1.10353, 0.84005, 0.61903, 0.46259, 0.35545, 0.2483, 0.14115, 0.034, 0, 0, 0, 0 };

	G4int j;

	for (j = 0; j < nEntries; j++)
	{
		RefractiveScin[j] = 1.58;
		AbsLengthScin[j] = 1. * m;
		PhotonEnergy[j] = PhotonEnergy[j] * eV;
	}

	G4MaterialPropertiesTable* ScintillatorProperties = new G4MaterialPropertiesTable();
	ScintillatorProperties->AddProperty("RINDEX", PhotonEnergy, RefractiveScin, nEntries);
	ScintillatorProperties->AddProperty("ABSLENGTH", PhotonEnergy, AbsLengthScin, nEntries);
	ScintillatorProperties->AddProperty("SCINTILLATIONCOMPONENT1", PhotonEnergy, SpIzlStr, nEntries);
	ScintillatorProperties->AddProperty("SCINTILLATIONCOMPONENT2", PhotonEnergy, SpIzlStr, nEntries);
	ScintillatorProperties->AddConstProperty("RESOLUTIONSCALE", 1.0);
	ScintillatorProperties->AddConstProperty("SCINTILLATIONYIELD", 12000 / MeV); // 12000
	ScintillatorProperties->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 2.4 * ns);
	ScintillatorProperties->AddConstProperty("SCINTILLATIONTIMECONSTANT2", 5 * ns);
	ScintillatorProperties->AddConstProperty("SCINTILLATIONYIELD1", 1.0);
	ScintillatorProperties->AddConstProperty("SCINTILLATIONYIELD2", 0.0);
	Scint->SetMaterialPropertiesTable(ScintillatorProperties);
	Scint->GetIonisation()->SetBirksConstant(0.126 * mm / MeV);


	G4double EnergyOpt[10] = { 1.9 * eV, 2.2 * eV, 2.3 * eV, 2.4 * eV, 2.56 * eV, 2.66 * eV, 2.68 * eV, 3.69 * eV, 3.7 * eV, 4.0 * eV };
	G4double AbsLenOpt[10] = { 5.0 * m, 5.0 * m, 5.0 * m, 5.0 * m, 5.0 * m, 5.0 * m, 0.1 * mm, 0.1 * mm, 5.0 * m, 5.0 * m };
	G4double SpIzlOpt[10] = { 0.001, 0.05, 0.25, 0.7, 1., 1., 0., 0., 0., 0. };
	G4double RindexOptCore[10] = { 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59 };
	G4double RindexOptInCov[10] = { 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49 };
	G4double RindexOptOutCov[10] = { 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42 };

	//Core optical properties
	G4MaterialPropertiesTable* OptCore = new G4MaterialPropertiesTable();
	OptCore->AddProperty("RINDEX", EnergyOpt, RindexOptCore, 10);
	OptCore->AddProperty("WLSABSLENGTH", EnergyOpt, AbsLenOpt, 10);
	OptCore->AddProperty("WLSCOMPONENT", EnergyOpt, SpIzlOpt, 10);
	OptCore->AddConstProperty("WLSTIMECONSTANT", 0.5 * ns);
	PS->SetMaterialPropertiesTable(OptCore);

	//Inner cover optical properties
	G4MaterialPropertiesTable* OptInCov = new G4MaterialPropertiesTable();
	OptInCov->AddProperty("RINDEX", EnergyOpt, RindexOptInCov, 10);
	OptInCov->AddProperty("WLSABSLENGTH", EnergyOpt, AbsLenOpt, 10);
	OptInCov->AddProperty("WLSCOMPONENT", EnergyOpt, SpIzlOpt, 10);
	OptInCov->AddConstProperty("WLSTIMECONSTANT", 0.5 * ns);
	PMMA->SetMaterialPropertiesTable(OptInCov);

	//Outer cover optical properties
	G4MaterialPropertiesTable* OptOutCov = new G4MaterialPropertiesTable();
	OptOutCov->AddProperty("RINDEX", EnergyOpt, RindexOptOutCov, 10);
	OptOutCov->AddProperty("WLSABSLENGTH", EnergyOpt, AbsLenOpt, 10);
	OptOutCov->AddProperty("WLSCOMPONENT", EnergyOpt, SpIzlOpt, 10);
	OptOutCov->AddConstProperty("WLSTIMECONSTANT", 0.5 * ns);
	FP->SetMaterialPropertiesTable(OptOutCov);

	//Air optical properties
	G4double EnergyAir[2] = { 1.9 * eV, 4.0 * eV };
	G4double AbsLenAir[2] = { 5.0 * m,  5.0 * m };
	G4double RindAir[2] = { 1.0002926, 1.0002926 };

	G4MaterialPropertiesTable* AirPT = new G4MaterialPropertiesTable();
	AirPT->AddProperty("RINDEX", EnergyAir, RindAir, 2);
	AirPT->AddProperty("ABSLENGTH", EnergyAir, AbsLenAir, 2);
	Air->SetMaterialPropertiesTable(AirPT);

	//Photomultiplier
	G4double EnergyPhotCat[2] = { 1.9 * eV, 4.0 * eV };
	G4double AbsLenPhotCat[2] = { 5.0 * m,  5.0 * m };
	G4double RindPhotCat[2] = { 1.5, 1.5 };

	G4MaterialPropertiesTable* PhotCatPT = new G4MaterialPropertiesTable();
	PhotCatPT->AddProperty("RINDEX", EnergyPhotCat, RindPhotCat, 2);
	PhotCatPT->AddProperty("ABSLENGTH", EnergyPhotCat, AbsLenPhotCat, 2);
	PhotCat->SetMaterialPropertiesTable(PhotCatPT);

	//Border borosilicate glass - aluminium: mirror reflection
	G4double reflectivity[2] = { 0.8, 0.8 };
	G4double PhotonEnergyBord[2] = { 1.9 * eV, 4.0 * eV };

	G4OpticalSurface* OptPovPhot = new G4OpticalSurface("PovPhotocathode");
	OptPovPhot->SetType(dielectric_metal);
	OptPovPhot->SetFinish(ground);
	OptPovPhot->SetModel(unified);

	G4MaterialPropertiesTable* PovPhotCatPT = new G4MaterialPropertiesTable();
	PovPhotCatPT->AddProperty("REFLECTIVITY", PhotonEnergyBord, reflectivity, 2);
	OptPovPhot->SetMaterialPropertiesTable(PovPhotCatPT);

	/* G4LogicalBorderSurface* PhotCatSurface = new G4LogicalBorderSurface("PhotoCathodeSurface", physPhotGlass, physPhot, OptPovPhot); */







	/*	VISUAL PROPERTIES	*/

	//Making world invisible
	auto UniverseVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
	UniverseVisAtt->SetVisibility(true);
	UniverseVisAtt->SetForceWireframe(true);
	logicWorld->SetVisAttributes(UniverseVisAtt);
	logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());

	//Making holes invisible
	G4int hle;
	for (level = 0; level < NLvls; level++)
	{
		for (column = 0; column < NCols; column++)
		{
			for (hle = 0; hle < Nhls; hle++)
			{
				logicHole[level][column][HoleCount]->SetVisAttributes(UniverseVisAtt);
				logicHole[level][column][HoleCount]->SetVisAttributes(G4VisAttributes::GetInvisible());
			}
		}
	}

	return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

