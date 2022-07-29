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

#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpBoundaryProcess.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double Z0const;

DetDetectorConstruction::DetDetectorConstruction()
: G4VUserDetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetDetectorConstruction::~DetDetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetDetectorConstruction::Construct()
{
  G4double a, z;  //Atomic mass, atomic number
  G4double density, fractionmass;
  G4int ncomponents, nelements;

  G4double ShirSteel = 1*mm; //Stainless Steel Width

  G4Element* elH  = new G4Element("Hydrogen", "H", z=1.,  a = 1.01*g/mole);
  G4Element* elC  = new G4Element("Carbon", "C", z=6.,  a = 12.01*g/mole);
  G4Element* elN  = new G4Element("Nitrogen", "N", z=7.,  a = 14.01*g/mole);
  G4Element* elO  = new G4Element("Oxygen", "O", z=8.,  a = 16.00*g/mole);
  G4Element* elNa = new G4Element("Natrium", "Na", z=11., a = 22.99*g/mole);
  G4Element* elMg= new G4Element("Magnesium", "Mg", z=12., a = 24.31*g/mole);
  G4Element* elSi = new G4Element("Silicium", "Si", z=14., a = 28.09*g/mole);
  G4Element* elCa = new G4Element("Calcium", "Ca", z=20., a = 40.08*g/mole);
  G4Element* elCr = new G4Element("Chromium", "Cr", z=24.,  a = 52.00*g/mole);
  G4Element* elFe = new G4Element("Ferrum", "Fe", z=26.,  a = 55.85*g/mole);


  //Lime glass (soda-lime-glass) PMT
  G4Material* SiO2 = new G4Material("MSiO2", density=2.2*g/cm3, nelements=2);
  SiO2->AddElement(elSi,  nelements=1);
  SiO2->AddElement(elO, nelements=2);
  G4Material* Na2O = new G4Material("MNa2O", density=2.3*g/cm3, nelements=2);
  Na2O->AddElement(elNa,  nelements=2);
  Na2O->AddElement(elO, nelements=1);
  G4Material* CaO = new G4Material("MCaO", density=3.34*g/cm3, nelements=2);
  CaO->AddElement(elCa,  nelements=1);
  CaO->AddElement(elO, nelements=1);
  G4Material* MgO = new G4Material("MMgO", density=3.6*g/cm3, nelements=2);
  MgO->AddElement(elMg,  nelements=1);
  MgO->AddElement(elO, nelements=1);

  G4Material* LimeGlass = new G4Material("MLimeGlass", density=2.53*g/cm3, ncomponents=4);
  LimeGlass->AddMaterial(SiO2, fractionmass=73.*perCent);
  LimeGlass->AddMaterial(Na2O, fractionmass=14.*perCent);
  LimeGlass->AddMaterial(CaO, fractionmass=9.*perCent);
  LimeGlass->AddMaterial(MgO, fractionmass=4.*perCent);

  //Aluminium
  G4Material* AlMaterial = new G4Material("MAluminium", z = 13., a = 26.98*g/mole, density = 2.7*g/cm3);

  //Stainless Steel
  G4Material* StSteel = new G4Material("MStSteel", density= 7.70*g/cm3, ncomponents=2);
  StSteel->AddElement(elFe, fractionmass=0.87);
  StSteel->AddElement(elCr, fractionmass=0.13);

  //Atmosphere (air)
  G4Material* Air = new G4Material("MAir", density= 1.290*mg/cm3, ncomponents=2);
  Air->AddElement(elN, fractionmass=0.7);
  Air->AddElement(elO, fractionmass=0.3);

  //Plexiglass
  G4Material* Plexiglass = new G4Material("MPlexiglass", density=1.19*g/cm3, nelements=3);
  Plexiglass->AddElement(elC, nelements=5);
  Plexiglass->AddElement(elO, nelements=2);
  Plexiglass->AddElement(elH, nelements=8);

  //Scintillator
  G4Material* Scint = new G4Material("MScint", density=1.032*g/cm3, ncomponents=2);
  Scint->AddElement(elC, fractionmass=0.476);
  Scint->AddElement(elH, fractionmass=0.524);


  //Optical

  //Optical properties Air
  G4double EnergyAir[2] = {2.3*eV, 3.2*eV};
  G4double RefractiveAir[2] = {1.0002926, 1.0002926};
  G4double AbsLengthAir[2]  = {5*m, 5*m};

  G4MaterialPropertiesTable* PT_Air = new G4MaterialPropertiesTable();
  PT_Air->AddProperty("RINDEX",   EnergyAir, RefractiveAir, 2);
  PT_Air->AddProperty("ABSLENGTH",EnergyAir, AbsLengthAir,  2);
  Air->SetMaterialPropertiesTable(PT_Air);

/*
//!!!!!!!!!!!!
  G4double EnergySc[2] = {2.3*eV, 3.2*eV};
  G4double RefractiveSc[2] = {1.58, 1.58};
  G4double AbsLengthSc[2]  = {1.6*m, 1.6*m};

  G4MaterialPropertiesTable* PT_Sc = new G4MaterialPropertiesTable();
  PT_Sc->AddProperty("RINDEX",   EnergySc, RefractiveSc, 2);
  PT_Sc->AddProperty("ABSLENGTH", EnergySc, AbsLengthSc,  2);
  Scint->SetMaterialPropertiesTable(PT_Sc);
//!!!!!!!!!!!!
*/

  //Optical properties Scintillator
  const G4int nEntries = 60;
  G4double PhotonEnergy[nEntries]={2.3, 2.31525, 2.33051, 2.34576, 2.36102, 2.37627, 2.39153, 2.40678, 2.42203, 2.43729, 2.45254, 2.4678, 2.48305, 2.49831, 2.51356,
       2.52881, 2.54407, 2.55932, 2.57458, 2.58983, 2.60508, 2.62034, 2.63559, 2.65085, 2.6661, 2.68136, 2.69661, 2.71186, 2.72712, 2.74237, 
       2.75763, 2.77288, 2.78814, 2.80339, 2.81864, 2.8339, 2.84915, 2.86441, 2.87966, 2.89492, 2.91017, 2.92542, 2.94068, 2.95593, 2.97119, 
       2.98644, 3.00169, 3.01695, 3.0322, 3.04746, 3.06271, 3.07797, 3.09322, 3.10847, 3.12373, 3.13898, 3.15424, 3.16949, 3.18475, 3.2};
  G4double RefractiveScin[nEntries];
  G4double AbsLengthScin[nEntries];
  G4double SpIzlStr[nEntries]={0, 0, 0.04304, 0.09311, 0.14318, 0.19325, 0.24331, 0.29338, 0.34345, 0.39352, 0.44359, 0.49365, 0.54372, 0.59379, 0.65703, 
       0.72516, 0.7829, 0.85487, 0.93619, 1.0156, 1.10002, 1.19322, 1.29936, 1.41172, 1.53233, 1.65876, 1.79893, 1.98186, 2.18771, 2.4366, 
       2.78324, 3.0698, 3.27276, 3.39218, 3.46918, 3.4941, 3.52619, 3.60856, 3.88683, 4.28688, 4.71702, 4.93565, 4.80817, 4.56821, 4.23367, 
       3.56117, 2.30136, 1.47323, 1.10353, 0.84005, 0.61903, 0.46259, 0.35545, 0.2483, 0.14115, 0.034, 0, 0, 0, 0};

  G4int j;

  for(j=0; j<nEntries; j++)
  {
    RefractiveScin[j] = 1.58;
    AbsLengthScin[j] = 1.*m;
    PhotonEnergy[j]=PhotonEnergy[j]*eV;
  }

  G4MaterialPropertiesTable* ScintillatorProperties = new G4MaterialPropertiesTable();
  ScintillatorProperties->AddProperty("RINDEX",        PhotonEnergy, RefractiveScin, nEntries);
  ScintillatorProperties->AddProperty("ABSLENGTH",     PhotonEnergy, AbsLengthScin,  nEntries);
  ScintillatorProperties->AddProperty("SCINTILLATIONCOMPONENT1", PhotonEnergy, SpIzlStr,       nEntries);
  ScintillatorProperties->AddProperty("SCINTILLATIONCOMPONENT2", PhotonEnergy, SpIzlStr,       nEntries);
  ScintillatorProperties->AddConstProperty("RESOLUTIONSCALE", 1.0);
  ScintillatorProperties->AddConstProperty("SCINTILLATIONYIELD", 12000/MeV); // 12000
  ScintillatorProperties->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 2.4*ns);
  ScintillatorProperties->AddConstProperty("SCINTILLATIONTIMECONSTANT2", 5*ns);
  ScintillatorProperties->AddConstProperty("SCINTILLATIONYIELD1", 1.0);
  ScintillatorProperties->AddConstProperty("SCINTILLATIONYIELD2", 0.0);
  Scint->SetMaterialPropertiesTable(ScintillatorProperties);
  Scint->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

  //Optical properties Plexiglass
  G4double EnergyPlgl[2] = {2.3*eV, 3.2*eV};
  G4double RefractivePlgl[2] = {1.49, 1.49};
  G4double AbsLengthPlgl[2]  = {5*m, 5*m};

  G4MaterialPropertiesTable* PT_Plgl = new G4MaterialPropertiesTable(); 
  PT_Plgl->AddProperty("RINDEX",   EnergyPlgl, RefractivePlgl, 2);
  PT_Plgl->AddProperty("ABSLENGTH",EnergyPlgl, AbsLengthPlgl,  2);
  Plexiglass->SetMaterialPropertiesTable(PT_Plgl);

  //Optical properties Limeglass PMT
  G4double EnergyLgl[2] = {2.3*eV, 3.2*eV};
  G4double RefractiveLgl[2] = {1.54, 1.54};
  G4double AbsLengthLgl[2]  = {5*m, 5*m};

  G4MaterialPropertiesTable* PT_Lgl = new G4MaterialPropertiesTable(); 
  PT_Lgl->AddProperty("RINDEX",   EnergyLgl, RefractiveLgl, 2);
  PT_Lgl->AddProperty("ABSLENGTH",EnergyLgl, AbsLengthLgl,  2);
  LimeGlass->SetMaterialPropertiesTable(PT_Lgl);


  //Option to switch on/off

  G4bool checkOverlaps = true;

  //World
  G4double world_sizeX = 5*m;
  G4double world_sizeY = 5*m;
  G4double world_sizeZ = 5*m;

  G4Box* solidWorld = new G4Box("World_g", 0.5*world_sizeX, 0.5*world_sizeY, 0.5*world_sizeZ);
  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, Air, "World_l");
  G4VPhysicalVolume* physWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "World", 0, false, 0, checkOverlaps);

  //Part 1 (little solid parallelepiped on bottom)
  G4double Dlinab1  = 245*mm;
  G4double Shirinab1 = 146*mm;
  G4double Visotab1 = 20*mm;

  G4double b1_xpos = 0*mm;
  G4double b1_ypos = 0*mm;
  G4double b1_zpos = 0.5*Visotab1;

  G4Box* solidb1 = new G4Box("b1_g", 0.5*Dlinab1 , 0.5*Shirinab1, 0.5*Visotab1);
  G4LogicalVolume* logicb1 =  new G4LogicalVolume(solidb1, StSteel, "b1_l");
  G4VPhysicalVolume* physb1 = new G4PVPlacement(0, G4ThreeVector(b1_xpos, b1_ypos, b1_zpos), logicb1, "b1", logicWorld, false, 0, checkOverlaps);

  //Part 2 (solid parallelepiped between part 1 and pyramid(trapezoid))
  G4double Dlinab2  = 207*mm;
  G4double Shirinab2 = 107*mm;
  G4double Visotab2 = 118*mm;

  G4double b2_xpos = b1_xpos;
  G4double b2_ypos = b1_ypos;
  G4double b2_zpos = (Visotab1+0.5*Visotab2);

  G4Box* solidb2 = new G4Box("b2_g", 0.5*Dlinab2 , 0.5*Shirinab2, 0.5*Visotab2);
  G4LogicalVolume* logicb2 =  new G4LogicalVolume(solidb2, StSteel, "b2_l");
  G4VPhysicalVolume* physb2 = new G4PVPlacement(0, G4ThreeVector(b2_xpos, b2_ypos, b2_zpos), logicb2, "b2", logicWorld, false, 0, checkOverlaps);

  //Air in part 2
  G4double Dlinap2  = Dlinab2-2*ShirSteel;
  G4double Shirinap2 = Shirinab2-2*ShirSteel;
  G4double Visotap2 = Visotab2;

  G4Box* solidp2 = new G4Box("p2_g", 0.5*Dlinap2 , 0.5*Shirinap2, 0.5*Visotap2);
  G4LogicalVolume* logicp2 =  new G4LogicalVolume(solidp2, Air, "p2_l");
  G4VPhysicalVolume* physp2 = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicp2, "p2", logicb2, false, 0, checkOverlaps);

  //Big part of air in part 1
  G4double Dlinap11  = Dlinab1-2*ShirSteel;
  G4double Shirinap11 = Shirinab1-2*ShirSteel;
  G4double Visotap11 = Visotab1-2*ShirSteel;

  G4Box* solidp11 = new G4Box("p11_g", 0.5*Dlinap11 , 0.5*Shirinap11, 0.5*Visotap11);
  G4LogicalVolume* logicp11 =  new G4LogicalVolume(solidp11, Air, "p11_l");
  G4VPhysicalVolume* physp11 = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicp11, "p11", logicb1, false, 0, checkOverlaps);

  //Small part of air in part 1
  G4double Dlinap12  = Dlinap2;
  G4double Shirinap12 = Shirinap2;
  G4double Visotap12 = ShirSteel;

  G4Box* solidp12 = new G4Box("p12_g", 0.5*Dlinap12, 0.5*Shirinap12, 0.5*Visotap12);
  G4LogicalVolume* logicp12 =  new G4LogicalVolume(solidp12, Air, "p12_l");
  G4VPhysicalVolume* physp12 = new G4PVPlacement(0, G4ThreeVector(0, 0, 0.5*(Visotab1-Visotap12)), logicp12, "p12", logicb1, false, 0, checkOverlaps);

  //Part 3 (trapezoid(pyramid))
  G4double Dlina1b3  = 200*mm;
  G4double Shirina1b3 = 100*mm;
  G4double Dlina2b3  = 851*mm;
  G4double Shirina2b3 = 851*mm;
  G4double Visotab3 = 233*mm;

  G4double b3_xpos = b2_xpos;
  G4double b3_ypos = b2_ypos;
  G4double b3_zpos = (Visotab1+Visotab2+0.5*Visotab3);

  G4Trd* solidb3 = new G4Trd("b3_g", 0.5*Dlina1b3, 0.5*Dlina2b3, 0.5*Shirina1b3, 0.5*Shirina2b3, 0.5*Visotab3);
  G4LogicalVolume* logicb3 =  new G4LogicalVolume(solidb3, StSteel, "b3_l");
  G4VPhysicalVolume* physb3 = new G4PVPlacement(0, G4ThreeVector(b3_xpos, b3_ypos, b3_zpos), logicb3, "b3", logicWorld, false, 0, checkOverlaps);

  //Air in part 3
  G4double Dlina1p3  = (Dlina1b3-ShirSteel*std::sqrt(Dlina1b3*Dlina1b3+4*Visotab3*Visotab3)/Visotab3);
  G4double Shirina1p3 = (Shirina1b3-ShirSteel*std::sqrt(Shirina1b3*Shirina1b3+4*Visotab3*Visotab3)/Visotab3);
  G4double Dlina2p3  = (Dlina2b3-ShirSteel*std::sqrt(Dlina2b3*Dlina2b3+4*Visotab3*Visotab3)/Visotab3);
  G4double Shirina2p3 = (Shirina2b3-ShirSteel*std::sqrt(Shirina2b3*Shirina2b3+4*Visotab3*Visotab3)/Visotab3);
  G4double Visotap3 = Visotab3;

  G4Trd* solidp3 = new G4Trd("p3_g", 0.5*Dlina1p3, 0.5*Dlina2p3, 0.5*Shirina1p3, 0.5*Shirina2p3, 0.5*Visotap3);
  G4LogicalVolume* logicp3 =  new G4LogicalVolume(solidp3, Air, "p3_l");
  G4VPhysicalVolume* physp3 = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicp3, "p3", logicb3, false, 0, checkOverlaps);

  //Part 4 (solid parallelepiped on top)
  G4double Dlinab4  = 858*mm;
  G4double Shirinab4 = 858*mm;
  G4double Visotab4 = 110*mm;

  G4double b4_xpos = b3_xpos;
  G4double b4_ypos = b3_ypos;
  G4double b4_zpos = (Visotab1+Visotab2+Visotab3+0.5*Visotab4);

  G4Box* solidb4 = new G4Box("b4_g", 0.5*Dlinab4 , 0.5*Shirinab4, 0.5*Visotab4);
  G4LogicalVolume* logicb4 =  new G4LogicalVolume(solidb4, StSteel, "b4_l");
  G4VPhysicalVolume* physb4 = new G4PVPlacement(0, G4ThreeVector(b4_xpos, b4_ypos, b4_zpos), logicb4, "b4", logicWorld, false, 0, checkOverlaps);

  //Big part of air in part 4
  G4double Dlinap41  = Dlinab4-2*ShirSteel;
  G4double Shirinap41 = Shirinab4-2*ShirSteel;
  G4double Visotap41 = Visotab4-2*ShirSteel;

  G4Box* solidp41 = new G4Box("p41_g", 0.5*Dlinap41 , 0.5*Shirinap41, 0.5*Visotap41);
  G4LogicalVolume* logicp41 =  new G4LogicalVolume(solidp41, Air, "p41_l");
  G4VPhysicalVolume* physp41 = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicp41, "p41", logicb4, false, 0, checkOverlaps);

  //Small part of air in part 4
  G4double Dlinap42  = Dlina2p3;
  G4double Shirinap42 = Shirina2p3;
  G4double Visotap42 = ShirSteel;

  G4Box* solidp42 = new G4Box("p42_g", 0.5*Dlinap42, 0.5*Shirinap42, 0.5*Visotap42);
  G4LogicalVolume* logicp42 =  new G4LogicalVolume(solidp42, Air, "p42_l");
  G4VPhysicalVolume* physp42 = new G4PVPlacement(0, G4ThreeVector(0, 0, -0.5*(Visotab4-Visotap42)), logicp42, "p42", logicb4, false, 0, checkOverlaps);

  //OrgBalka
  G4double Shirinaob = 20*mm;
  G4double Dlinaob = Shirinap41;
  G4double Visotaob = 50*mm;

  G4double otstup = 14*mm;//otstup org balki ot osnovania piramidi
  G4double rasst = 322*mm;//rasstoianie mezhdu org balkami

  G4Box* solidob = new G4Box("ob_g", 0.5*Dlinaob, 0.5*Shirinaob, 0.5*Visotaob);
  G4LogicalVolume* logicob =  new G4LogicalVolume(solidob, Plexiglass, "ob_l");
  G4VPhysicalVolume* physob1 = new G4PVPlacement(0, G4ThreeVector(0, -0.5*(rasst+Shirinaob), (-0.5*Visotap41+otstup+0.5*Visotaob)), logicob, "ob1", logicp41, false, 0, checkOverlaps);
  G4VPhysicalVolume* physob2 = new G4PVPlacement(0, G4ThreeVector(0, 0.5*(rasst+Shirinaob), (-0.5*Visotap41+otstup+0.5*Visotaob)), logicob, "ob2", logicp41, false, 0, checkOverlaps);

  //Scintillator
  G4double Dlinasc = 800*mm;
  G4double Shirinasc = 800*mm;
  G4double Visotasc = 40*mm;

  G4Box* solidsc = new G4Box("sc_g", 0.5*Dlinasc, 0.5*Shirinasc, 0.5*Visotasc);
  G4LogicalVolume* logicsc =  new G4LogicalVolume(solidsc, Scint, "sc_l");
  G4VPhysicalVolume* physsc = new G4PVPlacement(0, G4ThreeVector(0, 0, (-0.5*Visotap41+otstup+Visotaob+0.5*Visotasc)), logicsc, "sc", logicp41, false, 0, checkOverlaps);
  
  Z0const = Visotab1+Visotab2+Visotab3+0.5*(Visotab4-Visotap41)+otstup+Visotaob+0.5*Visotasc;
  
  //PMT glass
  G4double IRglPMT = 0*mm;
  G4double ORglPMT = 38*mm;
  G4double VisotaglPMT = 12.75*mm;
  G4double rsaPMT = 100*mm; //distance between standart and additional PMT

  G4Tubs* solidglPMT = new G4Tubs("glPMT_g", IRglPMT, ORglPMT, 0.5*VisotaglPMT, 0*deg, 360*deg);
  G4LogicalVolume* logicglPMT = new G4LogicalVolume(solidglPMT, LimeGlass, "glPMT_l");
  G4VPhysicalVolume* physglPMT = new G4PVPlacement(0, G4ThreeVector(-0.5*rsaPMT, 0, 0.5*(Visotap2-VisotaglPMT)), logicglPMT, "glsPMT", logicp2, false, 0, checkOverlaps); //standart PMT
  //G4PhysicalVolume* physglPMT = new G4PVPlacement(0, G4ThreeVector(0.5*rsaPMT, 0, 0.5*(Visotap2-VisotaglPMT)), logicglPMT, "glaPMT", logicp2, false, 0, checkOverlaps); //additional PMT


  //Photo cathode Al
  G4double IRalPMT = 74.5*mm;
  //G4double IRalPMT = 74.9*mm; //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!real
  G4double ORalPMT = 75*mm;
  G4double phi0alPMT = 0*deg;
  G4double dphialPMT = 360*deg;
  G4double teta0alPMT = 0*deg;
  G4double dtetaalPMT = 30*deg;

  G4double dh = 3*mm; //min depth PMT LimeGlass

  G4Sphere* solidalPMT = new G4Sphere("alPMT_g", IRalPMT, ORalPMT, phi0alPMT, dphialPMT, teta0alPMT, dtetaalPMT);
  G4LogicalVolume* logicalPMT = new G4LogicalVolume(solidalPMT, AlMaterial, "alPMT_l");
  G4VPhysicalVolume* physalPMT = new G4PVPlacement(0, G4ThreeVector(0, 0, (0.5*VisotaglPMT-dh-ORalPMT+0.8)), logicalPMT, "alPMT", logicglPMT, false, 0, checkOverlaps);//+0.8 - kostil v ramkah obiema

  //Korpus diffuse reflection
  G4double reflectivity_korp[2]= {0.9, 0.9};
  G4double PhotonEnergyPov[2] = {2.3*eV, 3.2*eV};

  G4OpticalSurface *OpticalPovKorpusa = new G4OpticalSurface("PovKorpusa");
  OpticalPovKorpusa->SetModel(unified);
  OpticalPovKorpusa->SetType(dielectric_dielectric);
  OpticalPovKorpusa->SetFinish(groundfrontpainted);

  G4MaterialPropertiesTable *PovKorpusaSurfacePT = new G4MaterialPropertiesTable();
  PovKorpusaSurfacePT->AddProperty("REFLECTIVITY", PhotonEnergyPov, reflectivity_korp, 2);
  OpticalPovKorpusa->SetMaterialPropertiesTable(PovKorpusaSurfacePT);

  G4LogicalBorderSurface *PovDif3 =  new G4LogicalBorderSurface("KorpuSurface3",    physp3,   physb3, OpticalPovKorpusa);
  G4LogicalBorderSurface *PovDif41 =  new G4LogicalBorderSurface("KorpuSurface41",    physp41,   physb4, OpticalPovKorpusa);
  G4LogicalBorderSurface *PovDif42 =  new G4LogicalBorderSurface("KorpuSurface42",    physp42,   physb4, OpticalPovKorpusa);


  //Mirror surface of photo cathode
  G4double reflectivity_al[2]= {0.98, 0.98};
  G4double PhotonEnergyAL[2] = {2.3*eV, 3.2*eV};

  G4OpticalSurface *OpticalPovAL = new G4OpticalSurface("PovCathoda");
  OpticalPovAL->SetModel(unified);
  OpticalPovAL->SetType(dielectric_metal);
  OpticalPovAL->SetFinish(polished);

  G4MaterialPropertiesTable *PovALSurfacePT = new G4MaterialPropertiesTable();
  PovALSurfacePT->AddProperty("REFLECTIVITY", PhotonEnergyAL, reflectivity_al, 2);
  OpticalPovAL->SetMaterialPropertiesTable(PovALSurfacePT);

  G4LogicalBorderSurface *PovAL =  new G4LogicalBorderSurface("PMTCathodeSurface",    physglPMT,   physalPMT, OpticalPovAL);

  logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());

/*
  G4Colour grey     (0.5, 0.5, 0.5);
  G4Colour black    (0.0, 0.0, 0.0);
  G4Colour red      (1.0, 0.0, 0.0);
  G4Colour green    (0.0, 1.0, 0.0);
  G4Colour blue     (0.0, 0.0, 1.0);
  G4Colour cyan     (0.0, 1.0, 1.0);
  G4Colour magenta  (1.0, 0.0, 1.0);
  G4Colour yellow   (1.0, 1.0, 0.0);

  G4VisAttributes* VisAtt_b4= new G4VisAttributes(grey);
  VisAtt_b4->SetForceSolid(true);
  logicb4->SetVisAttributes(VisAtt_b4);

  G4VisAttributes* VisAtt_b3= new G4VisAttributes(grey);
  VisAtt_b3->SetForceSolid(true);
  logicb3->SetVisAttributes(VisAtt_b3);

  G4VisAttributes* VisAtt_b2= new G4VisAttributes(grey);
  VisAtt_b2->SetForceSolid(true);
  logicb2->SetVisAttributes(VisAtt_b2);

  G4VisAttributes* VisAtt_b1= new G4VisAttributes(grey);
  VisAtt_b1->SetForceSolid(true);
  logicb1->SetVisAttributes(VisAtt_b1);

  G4VisAttributes* VisAtt_ob= new G4VisAttributes(cyan);
  VisAtt_ob->SetForceSolid(true);
  logicob->SetVisAttributes(VisAtt_ob);

  G4VisAttributes* VisAtt_sc= new G4VisAttributes(yellow);
  VisAtt_sc->SetForceSolid(true);
  logicsc->SetVisAttributes(VisAtt_sc);

  G4VisAttributes* VisAtt_PMTgl= new G4VisAttributes(cyan);
  VisAtt_PMTgl->SetForceSolid(true);
  logicglPMT->SetVisAttributes(VisAtt_PMTgl);

  G4VisAttributes* VisAtt_PMTal= new G4VisAttributes(magenta);
  VisAtt_PMTal->SetForceSolid(true);
  logicalPMT->SetVisAttributes(VisAtt_PMTal);
*/  
  fScoringVolume = logicalPMT; 

  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
