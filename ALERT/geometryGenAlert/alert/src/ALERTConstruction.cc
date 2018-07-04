#include "ALERTConstruction.hh"

#include "G4Isotope.hh"
#include "G4UnitsTable.hh"
#include "G4Material.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4CachedMagneticField.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4UserLimits.hh"
#include "G4ExtrudedSolid.hh"
#include "G4AutoDelete.hh"
#include "G4GDMLParser.hh"
#include "G4ChordFinder.hh"

ALERT::ALERT() : 
  G4VUserDetectorConstruction(), 
  world_x( 15.0*m ),
  world_y( 15.0*m ),
  world_z( 20.0*m )
{}

ALERT::~ALERT()
{
  delete Scinti;
  delete Silicon;
  delete Ar;
  delete Lead;
  delete Air;
  delete Kapton;
  delete He_Target;
  delete He_ClearS;
  delete Xe_varPT;
  delete O;
  delete N;
  delete C;
  delete H;
  delete elHe;
  delete elXe;
  delete fRecoilChamber;
  delete fRecoilHodo3;
}

void ALERT::DefineMaterials()
{
  G4NistManager* nist = G4NistManager::Instance();

  G4double a       = 0.0;
  G4double z       = 0.0;
  G4double density = 0.0;
  G4int    nel     = 0;
  
  // -----------------------------------------------
  // Elements 
  H = new G4Element("Hydrogen", "H", z=1., a=  1.01*g/mole);
  C = new G4Element("Carbon",   "C", z=6., a= 12.01*g/mole);
  N = new G4Element("Nitrogen", "N", z=7., a= 14.01*g/mole);
  O = new G4Element("Oxygen",   "O", z=8., a= 16.00*g/mole);
  elW = new G4Element("elW","W", z=74., a= 183.84*g/mole);
  elHe  = new G4Element("Helium","He" , 2., a = 4.003*g/mole);

  // -----------------------------------------------
  // Elemental Materials
  Lead     = new G4Material("Lead", z    = 82., a = 207.19*g/mole, density = 11.35*g/cm3);
  Ar       = new G4Material("ArgonGas",z = 18., a = 39.95*g/mole, density  = 1.782*mg/cm3);
  Silicon  = new G4Material("Silicon", z = 14., a = 28.09*g/mole, density  = 2.33*g/cm3);
  Tungsten = new G4Material("Tungsten",density = 19.3*g/cm3,1);
  Tungsten->AddElement(elW,1);
  elXe  = new G4Element("Xenon","Xe" , 54., a = 131.293*g/mole);

  // -----------------------------------------------
  // Air
  Air = new G4Material("Air", density= 1.29*mg/cm3, nel=2);
  Air->AddElement(N, 70.*perCent);
  Air->AddElement(O, 30.*perCent);

  // -----------------------------------------------
  // Mylar
  Mylar = nist->FindOrBuildMaterial("G4_MYLAR");

  // -----------------------------------------------
  // Scintillator plastic
  Scinti = new G4Material("Scintillator", density= 1.032*g/cm3, nel=2);
  Scinti->AddElement(C, 9);
  Scinti->AddElement(H, 10);

  // -----------------------------------------------
  // Helium for the target (3 atm)
  G4int    ncomponents;
  G4double fractionmass;
  G4double temperature, pressure;

  G4double a_noUnit     = 4.003;
  G4double pre_noUnit   = 3.0;
  G4double tempe_noUnit = 298.;
  pressure              = 3.0*atmosphere;
  temperature           = 298.*kelvin;
  density     = (a_noUnit*pre_noUnit)/(0.0821*tempe_noUnit)*kg/m3; //0.164*kg/m3 at 1 atm;
  He_Target = new G4Material("He", density, ncomponents=1, kStateGas,temperature,pressure);
  He_Target->AddElement(elHe, fractionmass=1.);

  // -----------------------------------------------
  // Helium for the clear space around the target (1 atm)
  pressure    = 1.0*atmosphere;
  pre_noUnit = 1.0;
  temperature = 298.*kelvin;
  tempe_noUnit = 298.;
  G4double d_He = (a_noUnit*pre_noUnit)/(0.0821*tempe_noUnit)*kg/m3; //0.164*kg/m3 at 1 atm;
  He_ClearS = new G4Material("He_ClearS", d_He, ncomponents=1, kStateGas,temperature,pressure);
  He_ClearS->AddElement(elHe, fractionmass=1.);

  // -----------------------------------------------
  // Xenon
  a_noUnit     = 131.293;
  pressure     = 0.3*atmosphere;
  pre_noUnit   = 0.3;
  temperature  = 298.*kelvin;
  tempe_noUnit = 298.;
  density      = (a_noUnit*pre_noUnit)/(0.0821*tempe_noUnit)*kg/m3; //5.37*kg/m3 at 1 atm;
  Xe_varPT = new G4Material("Xe", density, ncomponents=1,
			    kStateGas,temperature,pressure);
  Xe_varPT->AddElement(elXe, fractionmass=1.);
 

  // -----------------------------------------------
  // iso-Butane
  a = 58.122*g/mole;
  a_noUnit = 58.122;
  pressure    = 1.0*atmosphere;
  pre_noUnit = 1.0;        
  temperature = 298.*kelvin;
  tempe_noUnit = 298.;
  G4double d_iC4H10  = (a_noUnit*pre_noUnit)/(0.0821*tempe_noUnit)*kg/m3; 
  isobutane = new G4Material("isoC4H10",d_iC4H10,ncomponents=2,kStateGas,temperature,pressure);
  isobutane->AddElement(C,4);
  isobutane->AddElement(H,10);
  //G4cout << "Isobutane density is " << density/(kg/m3) << "." << G4endl;

  // -----------------------------------------------
  // Deuterium
  a = 2.0141*g/mole;
  a_noUnit = 2.0141;
  pressure = 7.0*atmosphere; 
  pre_noUnit = 7.0; 
  temperature = 298.*kelvin;
  tempe_noUnit = 298.;
  G4Isotope * iso_H2 = new G4Isotope("H2",1,2,a);
  ele_D = new G4Element("Deuterium Atom","D",1);
  ele_D->AddIsotope(iso_H2, 100.*perCent);
  density = (2.*a_noUnit*pre_noUnit)/(0.0821*tempe_noUnit)*kg/m3;
  Deuterium = new G4Material("DeuteriumGas", density, ncomponents=1, kStateGas, temperature, pressure); 
  Deuterium->AddElement(ele_D,100.*perCent);

  // -----------------------------------------------
  // LH2
  density   = 0.07085*g/cm3; //g/cm^3
  LH2 = new G4Material("LH2", density, ncomponents=1); 
  LH2->AddElement(H,100.*perCent);

  // -----------------------------------------------
  // Kapton
  Kapton = new G4Material("Kapton",density = 1.42*g/cm3, nel=4);
  Kapton->AddElement(H, fractionmass = 0.0273);
  Kapton->AddElement(C, fractionmass = 0.7213);
  Kapton->AddElement(N, fractionmass = 0.0765);
  Kapton->AddElement(O, fractionmass = 0.1749);

  // -----------------------------------------------
  // CO2
  a_noUnit = 44.01;
  pressure = 1.*atmosphere;
  pre_noUnit = 1.;
  temperature = 298.*kelvin;
  tempe_noUnit = 298.;
  G4double d_CO2 = (a_noUnit*pre_noUnit)/(0.0821*tempe_noUnit)*kg/m3;
  CO2 = new G4Material("CO2", d_CO2, nel=2,
		       kStateGas,temperature,pressure);
  CO2->AddElement(C,1);
  CO2->AddElement(O,2);

  
  // -----------------------------------------------
  // Gas mixtures

  // 90% He - 10% CO2
  density = (d_He*90./100. + d_CO2*10./100.);
  He10CO2 = new G4Material("He10CO2"  , density, ncomponents=2);
  He10CO2->AddElement(elHe, fractionmass = 0.9);
  He10CO2->AddMaterial(CO2,    fractionmass = 0.1);

  // 90% He - 10% iC4H10
  density = (d_He*90./100. + d_iC4H10*10./100.);
  HeiC4H10 = new G4Material("HeC4H10"  , density, ncomponents=2);
  HeiC4H10->AddElement(elHe,   fractionmass = 0.9);
  HeiC4H10->AddMaterial(isobutane, fractionmass = 0.1);
}

G4VPhysicalVolume* ALERT::Construct()
{  

  DefineMaterials();

  G4NistManager* nist = G4NistManager::Instance();

  bool    checkOverlaps    = false;
  double  red              = 0.0;
  double  green            = 0.0;
  double  blue             = 0.0;
  double  alpha            = 0.0;


  // ------------------------------------------------------------------------
  // World
  // ------------------------------------------------------------------------
  red         = 0.0/256.0;
  green       = 200.0/256.0;
  blue        = 0.0/256.0;
  alpha       = 0.4;

  world_mat   = nist->FindOrBuildMaterial("G4_AIR");
  if(!world_solid) world_solid = new G4Box( "World", 0.5*world_x, 0.5 * world_y, 0.5 * world_z );
  if(!world_log)   world_log   = new G4LogicalVolume( world_solid, world_mat, "world_log" );
  if(!world_phys)  world_phys  = new G4PVPlacement( 0, G4ThreeVector(), world_log, "world_phys", 0, false, 0, checkOverlaps );

  G4Colour            world_color {red, green, blue, alpha }; 
  G4VisAttributes   * world_vis   = new G4VisAttributes(world_color);
  world_vis->SetForceWireframe(true);
  world_log->SetVisAttributes(world_vis);

  // ------------------------------------------------------------------------
  // Target
  // ------------------------------------------------------------------------
  G4ThreeVector target_pos {target_posX, target_posY, target_posZ};
  G4Tubs* target = new G4Tubs("target_solid", innerRadiusOfTheTarget,outerRadiusOfTheTarget, fTargetLength/2.0, 0.0, 360.0*CLHEP::degree);
  G4Material * target_mat = He_Target; //LH2;//Deuterium;
  G4LogicalVolume* logicTarget = new G4LogicalVolume(target, target_mat,"target_log");
  new G4PVPlacement(0, target_pos,   logicTarget,  "Target_phys", world_log,    
		    false,           //no boolean operation
		    0,               //copy number
		    checkOverlaps);  //overlaps checking
  //Definition of visualisation attributes
  //Instantiation of a set of visualization attributes with cyan colour
  G4VisAttributes * TargetVisAtt = new G4VisAttributes(G4Colour(1.,1.,1.,0.3));
  // Set the forced wireframe style 
  TargetVisAtt->SetForceWireframe(true);
  // Assignment of the visualization attributes to the logical volume
  logicTarget->SetVisAttributes(TargetVisAtt);

  // ------------------------------------------------------------------------
  // Kapton foil around the target
  // ------------------------------------------------------------------------
  G4ThreeVector kapton_pos = G4ThreeVector(target_posX, target_posY, target_posZ);
  G4Tubs* kapton_cyl = new G4Tubs("KaptonCylinder", innerRadiusOfTheKapton, outerRadiusOfTheKapton, hightOfTheKapton, startAngleOfTheKapton, spanningAngleOfTheKapton);

  G4LogicalVolume* logicKapton =                         
    new G4LogicalVolume(kapton_cyl,          //its solid
			Kapton,              //its material
			"KaptonCylinder");   //its name

  new G4PVPlacement(0,                     //no rotation
		    kapton_pos,              //at position
		    logicKapton,             //its logical volume
		    "KaptonCylinder",        //its name
		    world_log,    //its mother  volume
		    false,                   //no boolean operation
		    0,                       //copy number
		    checkOverlaps);          //overlaps checking

  // Definition of visualisation attributes
  // Instantiation of a set of visualization attributes with cyan colour
  G4VisAttributes * KaptonVisAtt = new G4VisAttributes(G4Colour(1.,1.,0.));
  // Set the forced wireframe style 
  KaptonVisAtt->SetForceWireframe(true);
  // Assignment of the visualization attributes to the logical volume
  logicKapton->SetVisAttributes(KaptonVisAtt);


  // ------------------------------------------------------------------------
  // Clear Space - Parameters of space around the target
  // ------------------------------------------------------------------------
  G4ThreeVector around_pos = G4ThreeVector(target_posX, target_posY, target_posZ);
  G4Tubs* around = new G4Tubs("ClearSpace", innerRadiusOfAround, outerRadiusOfAround, hightOfAround, startAngleOfAround, spanningAngleOfAround);
  G4LogicalVolume* logicAround =                         
    new G4LogicalVolume(around,              //its solid
			He_ClearS,           //its material
			"ClearSpace");       //its name
  new G4PVPlacement(0,                     //no rotation
		    around_pos,              //at position
		    logicAround,             //its logical volume
		    "ClearSpace",            //its name
		    world_log,    //its mother  volume
		    false,                   //no boolean operation
		    0,                       //copy number
		    checkOverlaps);          //overlaps checking

  // Definition of visualisation attributes
  // Instantiation of a set of visualization attributes with cyan colour
  G4VisAttributes * AroundVisAtt = new G4VisAttributes(G4Colour(1.,0.,1.));
  // Set the forced wireframe style 
  AroundVisAtt->SetForceWireframe(true);
  // Assignment of the visualization attributes to the logical volume
  logicAround->SetVisAttributes(AroundVisAtt);


  // ------------------------------------------------------------------------
  // Mylar layer around clear space
  // ------------------------------------------------------------------------
  G4ThreeVector Oclkapton_pos = G4ThreeVector(target_posX, target_posY, target_posZ);

  G4Tubs* Oclkapton_cyl = new G4Tubs("OclKaptonCylinder",
				     innerRadiusOfTheOclKapton, outerRadiusOfTheOclKapton,
				     hightOfTheOclKapton, startAngleOfTheOclKapton,
				     spanningAngleOfTheOclKapton);

  G4LogicalVolume* logicOclKapton =                         
    new G4LogicalVolume(Oclkapton_cyl,          //its solid
			Mylar,                  //its material
			"OclKaptonCylinder");   //its name

  new G4PVPlacement(0,                     //no rotation
		    Oclkapton_pos,              //at position
		    logicOclKapton,             //its logical volume
		    "OclKaptonCylinder",        //its name
		    world_log,    //its mother  volume
		    false,                   //no boolean operation
		    0,                       //copy number
		    checkOverlaps);          //overlaps checking

  // Definition of visualisation attributes
  // Instantiation of a set of visualization attributes with cyan colour
  G4VisAttributes * OclKaptonVisAtt = new G4VisAttributes(G4Colour(1.,1.,0.));
  // Set the forced wireframe style 
  KaptonVisAtt->SetForceWireframe(true);
  // Assignment of the visualization attributes to the logical volume
  logicOclKapton->SetVisAttributes(OclKaptonVisAtt);
   
  // ------------------------------------------------------------------------
  // Outisde Mylar layer
  // ------------------------------------------------------------------------
  G4ThreeVector okapton_pos = G4ThreeVector(target_posX, target_posY, target_posZ);
  G4Tubs* okapton_cyl = new G4Tubs("OKaptonCylinder", 
				   innerRadiusOfTheOKapton, outerRadiusOfTheOKapton, 
				   hightOfTheOKapton, startAngleOfTheOKapton, spanningAngleOfTheOKapton);
  G4LogicalVolume* logicOKapton = new G4LogicalVolume(okapton_cyl, Mylar, "OKaptonCylinder");
  new G4PVPlacement(0,                     //no rotation
		    okapton_pos,              //at position
		    logicOKapton,             //its logical volume
		    "OKaptonCylinder",        //its name
		    world_log,    //its mother  volume
		    false,                   //no boolean operation
		    0,                       //copy number
		    checkOverlaps);          //overlaps checking

  // Definition of visualisation attributes
  // Instantiation of a set of visualization attributes with cyan colour
  G4VisAttributes * OKaptonVisAtt = new G4VisAttributes(G4Colour(1.,1.,0.));
  // Set the forced wireframe style 
  OKaptonVisAtt->SetForceWireframe(true);
  // Assignment of the visualization attributes to the logical volume
  logicOKapton->SetVisAttributes(OKaptonVisAtt);

  // Recoil Chamber
  fRecoilChamber = new RecoilChamberDetectorGeometry();
  fRecoilChamber->He10CO2 = He10CO2;
  fRecoilChamber->HeiC4H10 = HeiC4H10;
  fRecoilChamber->Tungsten = Tungsten; 
  fRecoilChamber->Mylar = Mylar;
  fRecoilChamber->PlaceParallelPhysicalVolume(world_log);
  
  // Recoil Hodoscope
  fRecoilHodo3 = new RecoilHodoDetectorGeometry3();
  fRecoilHodo3->BuildLogicalVolumes();
  fRecoilHodo3->PlacePhysicalVolume(world_log, world_phys);
  
  G4GDMLParser parser;
  G4String gdml_name = "ALERTdetector.gdml";
  parser.Write(gdml_name, world_phys);

  return world_phys;
}
