#include "RecoilChamberDetectorGeometry.hh"

#include "CLHEP/Units/PhysicalConstants.h"

#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4NistManager.hh"
#include "G4UserLimits.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4Isotope.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4GenericTrap.hh"
#include "G4IntersectionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"

#include <cmath>

RecoilChamberDetectorGeometry::RecoilChamberDetectorGeometry()
{
   using namespace clas12::wire;
   using namespace CLHEP;

   G4SDManager* SDMan = G4SDManager::GetSDMpointer();
}
//______________________________________________________________________________

RecoilChamberDetectorGeometry::~RecoilChamberDetectorGeometry(){}
//______________________________________________________________________________

G4VSolid * RecoilChamberDetectorGeometry::BuildWireSolid(int layer, int subcell)
{
   G4RotationMatrix aROT;
   aROT.rotateZ(5.0*CLHEP::degree);
   std::string name = "trap" + std::to_string(layer) + "_sub" + std::to_string(subcell);
   std::vector<G4TwoVector> trap_points = fRCGeometry.GetSubCellTrapPoints(layer,subcell);
   G4VSolid * trap  = new  G4GenericTrap(name.c_str(), fRCGeometry.WireLength/2.0, trap_points);
   G4VSolid * box   = new G4Box("container_box", 4.0*mm, 4.0*mm, fRCGeometry.WireLength/2.0);
   G4VSolid * trap0 = new G4IntersectionSolid("trap0",box,trap,fRCGeometry.GetFirstCellTransform(layer,subcell));
   return trap0;
}
//______________________________________________________________________________

void RecoilChamberDetectorGeometry::BuildUnitCells()
{
   for(int lay = 0; lay<fRCGeometry.NLayers; lay++){

      fWireVolume_solid[lay][0] = BuildWireSolid(lay,0);
      fWireVolume_solid[lay][1] = BuildWireSolid(lay,1);
      fWireVolume_solid[lay][2] = BuildWireSolid(lay,2);
      fWireVolume_solid[lay][3] = BuildWireSolid(lay,3);

      //--Step max
      G4double maxStep = 1.0*mm;
      G4UserLimits * fStepLimit = new G4UserLimits(maxStep);
      fStepLimit->SetMaxAllowedStep(maxStep);

      for(int i = 0; i<4; i++){

         std::string log_name = "wire_volume_log_" 
            + std::to_string(lay) + "_part" + std::to_string(i);

         fWireVolume_log[lay][i] = new G4LogicalVolume(
               fWireVolume_solid[lay][i], 
               HeiC4H10, 
               log_name);
         fWireVolume_log[lay][i]->SetUserLimits(fStepLimit);
      }
   }
}
//______________________________________________________________________________

void  RecoilChamberDetectorGeometry::PlaceCells(G4LogicalVolume * mother, int layer, double z_rotation, int wire_number ){

   bool check_overlaps = false;

   G4VisAttributes * vs = 0;
   if(layer%2==0) {
      vs = new G4VisAttributes(G4Colour(0.0,0.9,0.1));
   } else {
      vs = new G4VisAttributes(G4Colour(0.9,0.0,0.1));
   }
   vs->SetForceWireframe(true);

   for(int i = 0; i<4; i++) {
      std::string wire_name = "wire_volume_" + std::to_string(wire_number) + "_part" + std::to_string(i);
      fWireVolume_log[layer][i]->SetVisAttributes(vs);//G4VisAttributes::GetInvisible());

      new G4PVPlacement(
            G4Transform3D(
               fRCGeometry.GetSenseWireRotation(wire_number),
               fRCGeometry.GetSenseWirePosition(wire_number)),
            fWireVolume_log[layer][i], wire_name,
            mother,
            false,
            wire_number,
            check_overlaps);
   }
}
//______________________________________________________________________________


void RecoilChamberDetectorGeometry::BuildLogicalVolumes()
{
   using namespace CLHEP;
   using namespace clas12::geo;
   using namespace clas12::geo::DC;
}
//______________________________________________________________________________

G4VPhysicalVolume * RecoilChamberDetectorGeometry::PlaceParallelPhysicalVolume(G4LogicalVolume * mother )
{
   using namespace clas12::geo;

   bool check_overlaps = false;

   BuildUnitCells();

   //--Creating the geometry
   G4ThreeVector gasDetector_pos = G4ThreeVector(gasDetector_posX, gasDetector_posY, gasDetector_posZ);

   G4Tubs* gasDetector = new G4Tubs("GasDetectorPara", 
         innerRadiusOfTheGasDetector,outerRadiusOfTheGasDetector, 
         hightOfTheGasDetector, 
         startAngleOfTheGasDetector, spanningAngleOfTheGasDetector);

   G4LogicalVolume* logicGasDetector =                         
      new G4LogicalVolume(gasDetector, HeiC4H10, "GasDetectorPara_log");            

   G4VPhysicalVolume * gasDetector_phys = new G4PVPlacement(0, gasDetector_pos, logicGasDetector, 
         "gasDetectorPara_phys", mother, false, 0, check_overlaps);            

   G4VisAttributes * GasDetectorVisAtt = new G4VisAttributes(G4Colour(0.3,0.1,0.1));
   GasDetectorVisAtt->SetForceWireframe(true);
   logicGasDetector->SetVisAttributes(GasDetectorVisAtt);
 
   int wire_number = 0;

   for(G4int lay=0; lay<fRCGeometry.NLayers; lay++){

      bool placed_one = false;

      G4double PhiWire = fRCGeometry.fCellDeltaPhi.at(lay);

      for(int wi=0;wi<fRCGeometry.fNCells.at(lay); wi++){  

         if( (!placed_one) ) {

            PlaceCells( logicGasDetector, lay, double(wi)*PhiWire, wire_number );

            wire_number++;
         }
      }   
   }

   return gasDetector_phys;
}
//______________________________________________________________________________

G4VPhysicalVolume * RecoilChamberDetectorGeometry::PlacePhysicalVolume(G4LogicalVolume * mother )
{

   using namespace clas12::geo;

   G4ThreeVector gasDetector_pos = G4ThreeVector(gasDetector_posX, gasDetector_posY, gasDetector_posZ);

   G4Tubs* gasDetector = new G4Tubs("GasDetector", 
         innerRadiusOfTheGasDetector,outerRadiusOfTheGasDetector, 
         hightOfTheGasDetector, 
         startAngleOfTheGasDetector, spanningAngleOfTheGasDetector);

   G4LogicalVolume* logicGasDetector =                         
      new G4LogicalVolume(gasDetector, HeiC4H10, "GasDetector");            

   G4VPhysicalVolume * gasDetector_phys = new G4PVPlacement(0, gasDetector_pos, logicGasDetector, 
         "GasDetector", mother, false, 0, true);            

   G4VisAttributes * GasDetectorVisAtt = new G4VisAttributes(G4Colour(0.3,0.1,0.1));
   GasDetectorVisAtt->SetForceWireframe(true);
   logicGasDetector->SetVisAttributes(GasDetectorVisAtt);

   //--Step max
   G4double maxStep = 2.0*mm;
   G4UserLimits * fStepLimit = new G4UserLimits(maxStep);
   fStepLimit->SetMaxAllowedStep(maxStep);
   logicGasDetector->SetUserLimits(fStepLimit);

   // --Placing the wires
   G4Tubs* wiresolid = new G4Tubs("Wire", innerRadiusOfTheWire,outerRadiusOfTheWire, fRCGeometry.WireLength/2.0, startAngleOfTheWire, spanningAngleOfTheWire);
   G4LogicalVolume* wires_even_log = new G4LogicalVolume(wiresolid, Tungsten, "wires_even_log");       
   G4LogicalVolume* wires_odd_log = new G4LogicalVolume(wiresolid, Tungsten, "wires_odd_log");       

   return gasDetector_phys;
}
//______________________________________________________________________________


