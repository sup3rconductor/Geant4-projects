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
	

	//Air
	G4Material* Air = new G4Material("MAir", density = 1.290 * mg / cm3, ncomponents = 2);
	Air->AddElement(elN, fractionmass = 0.8);
	Air->AddElement(elO, fractionmass = 0.2);

	//Aluminium
	G4Material* AlMaterial = new G4Material("MAluminium", z = 13., a = 26.98 * g / mole, density = 2.8 * g / cm3);

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
	G4double Dlinascr = 200 * mm;
	G4double Shirinascr = 200 * mm;
	G4double Visotascr = 5 * mm;


	//Hole parameters
	G4double Shirinaotv = 1 * mm;
	G4double Visotaotv =  4 * mm;

	G4double h1_xpos = 0 * mm;
	G4double h1_ypos = 0.5 * Shirinascr - 47 * mm;
	G4double h1_zpos = 0.5 * Visotascr;


	G4Box* solidScr = new G4Box("sc_s", 0.5 * Dlinascr, 0.5 * Shirinascr, 0.5 * Visotascr);
	G4Box* overstie = new G4Box("otv", 0.5 * Dlinascr + 0.1 * mm, 0.5 * Shirinaotv, 0.5 * Visotaotv);
	
	//Hole 1
	G4RotationMatrix* Rot = new G4RotationMatrix;
	G4ThreeVector trans1(h1_xpos, h1_ypos, h1_zpos);
	G4SubtractionSolid* solidScintplate1 = new G4SubtractionSolid("scintplate_s1", solidScr, overstie, Rot, trans1);

	G4double h2_ypos = h1_ypos - 35 * mm;
	G4ThreeVector trans2(h1_xpos, h2_ypos, h1_zpos);

	//Hole 2
	G4SubtractionSolid* solidScintplate2 = new G4SubtractionSolid("scintplate_s2", solidScintplate1, overstie, Rot, trans2);
	G4double h3_ypos = h2_ypos - 35 * mm;
	G4ThreeVector trans3(h1_xpos, h3_ypos, h1_zpos);

	//Hole 3
	G4SubtractionSolid* solidScintplate3 = new G4SubtractionSolid("scintplate_s3", solidScintplate2, overstie, Rot, trans3);
	G4double h4_ypos = h3_ypos - 35 * mm;
	G4ThreeVector trans4(h1_xpos, h4_ypos, h1_zpos);

	//Hole 4
	G4SubtractionSolid* solidScintplate4 = new G4SubtractionSolid("scintplate_s4", solidScintplate3, overstie, Rot, trans4);
	
	G4LogicalVolume* logicScintplate = new G4LogicalVolume(solidScintplate4, Scint, "scintplate_l");
	G4ThreeVector transplate(0., 0., 0.);

	G4VPhysicalVolume* physScintplate = new G4PVPlacement(0, transplate, logicScintplate, "scintplate", logicWorld, true, 0);


	//Optical fiber
	G4double OptRad = 0.5 * mm;
	G4double OptHeight = 100 * mm;
	G4double CovThickness = 0.03 * mm;
	G4RotationMatrix* OptRot = new G4RotationMatrix;
	OptRot->rotateY(90. * deg);
	Z0const = 0.5 * Visotascr - 0.5 * Visotaotv + OptRad;
	Y0const = h1_ypos;
	X0const = h1_xpos;

	//Core
	G4Tubs* solidCore = new G4Tubs("core_s", 0, OptRad - 2 * CovThickness, OptHeight, 0. * deg, 360. * deg);
	G4LogicalVolume* logicCore = new G4LogicalVolume(solidCore, PS, "core_l");
	G4VPhysicalVolume* physCore1 = new G4PVPlacement(OptRot, G4ThreeVector(h1_xpos, h1_ypos, 0.5 * Visotascr - 0.5 * Visotaotv + OptRad), logicCore, "CORE_1", logicWorld, false, 0, checkOverlaps);
	G4VPhysicalVolume* physCore2 = new G4PVPlacement(OptRot, G4ThreeVector(h1_xpos, h2_ypos, 0.5 * Visotascr - 0.5 * Visotaotv + OptRad), logicCore, "CORE_2", logicWorld, false, 0, checkOverlaps);
	G4VPhysicalVolume* physCore3 = new G4PVPlacement(OptRot, G4ThreeVector(h1_xpos, h3_ypos, 0.5 * Visotascr - 0.5 * Visotaotv + OptRad), logicCore, "CORE_3", logicWorld, false, 0, checkOverlaps);
	G4VPhysicalVolume* physCore4 = new G4PVPlacement(OptRot, G4ThreeVector(h1_xpos, h4_ypos, 0.5 * Visotascr - 0.5 * Visotaotv + OptRad), logicCore, "CORE_4", logicWorld, false, 0, checkOverlaps);

	//Inner cover
	G4Tubs* solidInCov = new G4Tubs("InCov_s", OptRad - 2 * CovThickness, OptRad - CovThickness, OptHeight, 0. * deg, 360. * deg);
	G4LogicalVolume* logicInCov = new G4LogicalVolume(solidInCov, PMMA, "InCov_l");
	G4VPhysicalVolume* physInCov1 = new G4PVPlacement(OptRot, G4ThreeVector(h1_xpos, h1_ypos, 0.5 * Visotascr - 0.5 * Visotaotv + OptRad), logicInCov, "INNER_COVER_1", logicWorld, false, 0, checkOverlaps);
	G4VPhysicalVolume* physInCov2 = new G4PVPlacement(OptRot, G4ThreeVector(h1_xpos, h2_ypos, 0.5 * Visotascr - 0.5 * Visotaotv + OptRad), logicInCov, "INNER_COVER_2", logicWorld, false, 0, checkOverlaps);
	G4VPhysicalVolume* physInCov3 = new G4PVPlacement(OptRot, G4ThreeVector(h1_xpos, h3_ypos, 0.5 * Visotascr - 0.5 * Visotaotv + OptRad), logicInCov, "INNER_COVER_3", logicWorld, false, 0, checkOverlaps);
	G4VPhysicalVolume* physInCov4 = new G4PVPlacement(OptRot, G4ThreeVector(h1_xpos, h4_ypos, 0.5 * Visotascr - 0.5 * Visotaotv + OptRad), logicInCov, "INNER_COVER_4", logicWorld, false, 0, checkOverlaps);

	//Outer cover
	G4Tubs* solidOutCov = new G4Tubs("OutCov_s", OptRad - CovThickness, OptRad, OptHeight, 0. * deg, 360. * deg);
	G4LogicalVolume* logicOutCov = new G4LogicalVolume(solidOutCov, FP, "OutCov_l");
	G4VPhysicalVolume* physOutCov2 = new G4PVPlacement(OptRot, G4ThreeVector(h1_xpos, h2_ypos, 0.5 * Visotascr - 0.5 * Visotaotv + OptRad), logicOutCov, "OUTER_COVER_2", logicWorld, false, 0, checkOverlaps);
	G4VPhysicalVolume* physOutCov1 = new G4PVPlacement(OptRot, G4ThreeVector(h1_xpos, h1_ypos, 0.5 * Visotascr - 0.5 * Visotaotv + OptRad), logicOutCov, "OUTER_COVER_1", logicWorld, false, 0, checkOverlaps);
	G4VPhysicalVolume* physOutCov3 = new G4PVPlacement(OptRot, G4ThreeVector(h1_xpos, h3_ypos, 0.5 * Visotascr - 0.5 * Visotaotv + OptRad), logicOutCov, "OUTER_COVER_3", logicWorld, false, 0, checkOverlaps);
	G4VPhysicalVolume* physOutCov4 = new G4PVPlacement(OptRot, G4ThreeVector(h1_xpos, h4_ypos, 0.5 * Visotascr - 0.5 * Visotaotv + OptRad), logicOutCov, "OUTER_COVER_4", logicWorld, false, 0, checkOverlaps);

	/* Photomultiplier */

	//Glass
	G4double GlassRad = 10 * mm;
	G4double GlassHeight = 1 * mm;
	G4double Glass_xpos = 0.5 * (Dlinascr + GlassHeight);

	G4Tubs* solidPhotGlass = new G4Tubs("PhotGlass_s", 0, GlassRad, 0.5 * GlassHeight, 0. * deg, 360. * deg);
	G4LogicalVolume* logicPhotGlass = new G4LogicalVolume(solidPhotGlass, PhotCat, "PhotGlass_s");
	G4VPhysicalVolume* physPhotGlass = new G4PVPlacement(OptRot, G4ThreeVector(Glass_xpos, h1_ypos, 0.5 * Visotascr - 0.5 * Visotaotv + OptRad), logicPhotGlass, "GLASS", logicWorld, false, 0, checkOverlaps);

	//PhotoCathode
	G4double PhotRad = 10 * mm;
	G4double PhotHeight = 0.1 * mm;
	G4double Phot_xpos = 0.5 * (Dlinascr + PhotHeight) + GlassHeight;

	G4Tubs* solidPhot = new G4Tubs("phot_s", 0, PhotRad, 0.5 * PhotHeight, 0. * deg, 360. * deg);
	G4LogicalVolume* logicPhot = new G4LogicalVolume(solidPhot, AlMaterial, "phot_l");
	G4VPhysicalVolume* physPhot = new G4PVPlacement(OptRot, G4ThreeVector(Phot_xpos, h1_ypos, 0.5 * Visotascr - 0.5 * Visotaotv + OptRad), logicPhot, "PHOTOCATHODE", logicWorld, false, 0, checkOverlaps);

	//Body
	G4double BodyRad = 11 * mm;
	G4double BodyHeight = 1.3 * mm;
	G4double HollowRad = 10 * mm;
	G4double HollowHeight = 1.1 * mm;
	G4double Hollow_zpos = 0.1 * mm;
	G4double Body_xpos = 0.5 * (Dlinascr + BodyHeight);

	G4Tubs* body = new G4Tubs("bdy", 0, BodyRad, 0.5 * BodyHeight, 0. * deg, 360. * deg);
	G4Tubs* hollow = new G4Tubs("hllw", 0, HollowRad, 0.5 * HollowHeight, 0. * deg, 360. * deg);
	G4ThreeVector transfer(0, 0, Hollow_zpos);
	G4SubtractionSolid* solidBody = new G4SubtractionSolid("body_s", body, hollow, Rot, transfer);
	G4LogicalVolume* logicBody = new G4LogicalVolume(solidBody, AlMaterial, "body_l");
	G4VPhysicalVolume* physBody = new G4PVPlacement(OptRot, G4ThreeVector(Body_xpos, h1_ypos, 0.5 * Visotascr - 0.5 * Visotaotv + OptRad), logicBody, "BODY", logicWorld, true, 0);

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
	/*Scint->SetMaterialPropertiesTable(ScintillatorProperties);
	Scint->GetIonisation()->SetBirksConstant(0.126 * mm / MeV); */

	
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
	G4double reflectivity[2] = {0.8, 0.8};
	G4double PhotonEnergyBord[2] = {1.9 * eV, 4.0 * eV};

	G4OpticalSurface* OptPovPhot = new G4OpticalSurface("PovPhotocathode");
	OptPovPhot->SetType(dielectric_metal);
	OptPovPhot->SetFinish(ground);
	OptPovPhot->SetModel(unified);

	G4MaterialPropertiesTable* PovPhotCatPT = new G4MaterialPropertiesTable();
	PovPhotCatPT->AddProperty("REFLECTIVITY", PhotonEnergyBord, reflectivity, 2);
	OptPovPhot->SetMaterialPropertiesTable(PovPhotCatPT);


	/*	VISUAL PROPERTIES	*/


	//Making world invisible
	auto UniverseVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
	UniverseVisAtt->SetVisibility(true);
	UniverseVisAtt->SetForceWireframe(true);
	logicWorld->SetVisAttributes(UniverseVisAtt);
	logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());

	return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
