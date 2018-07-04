#ifndef __ALERT_HH_
#define __ARLERT_HH_

#include "RecoilChamberDetectorGeometry.hh"
#include "RecoilHodoDetectorGeometry3.hh"

#include "G4VUserDetectorConstruction.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4String.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4UserLimits.hh"
#include "G4VSolid.hh"

class ALERT : public G4VUserDetectorConstruction{

private:

  G4double world_x;
  G4double world_y;
  G4double world_z;

  //--Parameters of the magnetic field
  G4double   fieldValue = 5*tesla;
  G4double   acc_inter = 0.1*mm;

  //--Parameters of the target
  G4double   innerRadiusOfTheTarget = 0.0*mm;
  G4double   outerRadiusOfTheTarget = 6.*mm;
  G4double   fTargetLength = 40.0*cm;//5.0*cm;
  G4double   startAngleOfTheTarget = 0.*deg;
  G4double   spanningAngleOfTheTarget = 360.*deg;
  G4double   target_posX = 0.*mm;
  G4double   target_posY = 0.*mm;
  G4double   target_posZ = 0.*mm;

  //--Kapton foil around the target
  G4double      innerRadiusOfTheKapton = 6.0*mm;
  G4double      outerRadiusOfTheKapton = 6.028*mm;
  G4double      hightOfTheKapton = 200.*mm;
  G4double      startAngleOfTheKapton = 0.*deg;
  G4double      spanningAngleOfTheKapton = 360.*deg;

  //--Parameters of space around the target
  G4double  innerRadiusOfAround = 6.028*mm;
  G4double  outerRadiusOfAround = 29.990*mm;
  G4double  hightOfAround = 200.*mm;
  G4double  startAngleOfAround = 0.*deg;
  G4double spanningAngleOfAround = 360.*deg;

  //--Mylar foil around the clear space
  G4double   innerRadiusOfTheOclKapton = 29.996*mm;
  G4double   outerRadiusOfTheOclKapton = 30.000*mm;
  G4double   hightOfTheOclKapton = 200.*mm;
  G4double   startAngleOfTheOclKapton = 0.*deg;
  G4double   spanningAngleOfTheOclKapton = 360.*deg;
  //--


  //--Mylar foil around the gas detector
  G4double     innerRadiusOfTheOKapton = 79.995*mm;
  G4double     outerRadiusOfTheOKapton = 80.0*mm;
  G4double     hightOfTheOKapton = 200.*mm;
  G4double     startAngleOfTheOKapton = 0.*deg;
  G4double     spanningAngleOfTheOKapton = 360.*deg;
  //--

  //--Parameters of the scintillators
  G4double  innerRadiusOfTheSiDetector = 80.*mm;
  G4double  outerRadiusOfTheSiDetector = 140.*mm;
  G4double  hightOfTheSiDetector = 200.*mm;
  G4double  startAngleOfTheSiDetector = 0.*deg;
  G4double  spanningAngleOfTheSiDetector = 360.*deg;
  G4double  SiDetector_posX = 0.*mm;
  G4double  SiDetector_posY = 0.*mm;
  G4double  SiDetector_posZ = 0.*mm; 

  G4Material* Air;
  G4Material* Ar;
  G4Material* Silicon;
  G4Material* Scinti;
  G4Material* Lead;
  G4Material* Kapton;
  G4Material* He_Target;
  G4Material* He_ClearS;
  G4Material* Xe_varPT;
  G4Material* isobutane;
  G4Material* Deuterium;
  G4Material* CO2;
  G4Material* He10CO2;
  G4Material* HeiC4H10;
  G4Material* Tungsten; 
  G4Material* Mylar;
  G4Material* LH2;

  G4Element* H;
  G4Element* C;
  G4Element* N;
  G4Element* O;
  G4Element* elHe;
  G4Element* elXe;
  G4Element* ele_D;
  G4Element* elW;
  G4UserLimits* fStepLimit; // pointer to user step limits

  G4Material        * world_mat    = nullptr;
  G4VSolid          * world_solid  = nullptr;
  G4LogicalVolume   * world_log    = nullptr;
  G4VPhysicalVolume * world_phys   = nullptr;

  RecoilChamberDetectorGeometry* fRecoilChamber = nullptr;
  RecoilHodoDetectorGeometry3* fRecoilHodo3 = nullptr;
  
protected:

  G4LogicalVolume  * fScoringVolume;
  G4String    fCollimatorMatName;

  void DefineMaterials();

public:

  ALERT();

  virtual ~ALERT();

  virtual G4VPhysicalVolume* Construct();
};
#endif

  
