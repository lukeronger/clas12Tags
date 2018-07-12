#include "RecoilHodoDetectorGeometry3.hh"

#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"
#include "G4GenericTrap.hh"
#include "G4Box.hh"
#include "G4IntersectionSolid.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"
#include "G4PVPlacement.hh"
#include "CLHEP/Units/PhysicalConstants.h"
#include "G4UserLimits.hh"
#include "G4TwistedTrd.hh"
#include "G4Trd.hh"
#include "G4Trap.hh"
#include "G4Tubs.hh"
#include "G4TwistedTrap.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4Isotope.hh"
#include "G4UnitsTable.hh"
#include "G4Material.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4UserLimits.hh"
#include "G4ExtrudedSolid.hh"

#include <cmath>

RecoilHodoDetectorGeometry3::RecoilHodoDetectorGeometry3()
{
   using namespace CLHEP;

   G4SDManager* SDMan = G4SDManager::GetSDMpointer();

   fWavelengths = {
      373.421917901*nm , 381.454526965*nm , 387.619063761*nm , 391.745703535*nm ,
      394.052066032*nm , 396.802820243*nm , 399.750855903*nm , 401.737907190*nm ,
      402.733466662*nm , 405.542185011*nm , 406.955696417*nm , 408.992576523*nm ,
      411.443340904*nm , 413.117182468*nm , 415.795735738*nm , 419.661028440*nm ,
      431.187756347*nm , 436.657740416*nm , 449.631537914*nm , 456.509948815*nm ,
      467.611606386*nm , 482.131114200*nm , 494.626961798*nm , 505.097115352*nm 
   };

   double hc_over_eV  = 1.239842e-6*m*eV;

   for(int i = 0; i < fWavelengths.size(); i++ ){
      fPhotonEnergies[i] = hc_over_eV/fWavelengths[i];
   }

   fScintLightOutput = {
      0.0000, 0.0423, 0.1954, 0.397,
      0.6152, 0.7504, 0.8936, 0.946,
      0.9682, 0.9882, 0.9800, 0.933,
      0.8646, 0.7392, 0.6178, 0.538,
      0.4390, 0.3719, 0.1973, 0.132,
      0.0769, 0.0315, 0.0064, 0.000
   };

   fRefIndex = {
      1.58, 1.58, 1.58, 1.58,
      1.58, 1.58, 1.58, 1.58,
      1.58, 1.58, 1.58, 1.58,
      1.58, 1.58, 1.58, 1.58,
      1.58, 1.58, 1.58, 1.58,
      1.58, 1.58, 1.58, 1.58
   };

   fScint1_pos  = {0.0, 0.0, 0.0};
   fScint2_pos  = {0.0, 0.0, 0.0};
   fContainer_pos = {0,0,0};
}
//______________________________________________________________________________

RecoilHodoDetectorGeometry3::~RecoilHodoDetectorGeometry3(){}

//______________________________________________________________________________

void RecoilHodoDetectorGeometry3::BuildMaterials()
{
   using std::vector;
   int     natoms           = 0;
   int     ncomponents      = 0;
   double  A                = 0.0;
   double  Z                = 0.0;
   double  density          = 0.0;
   double  pressure         = 0.0;
   double  temperature      = 0.0;
   double  a                = 0.0;

   vector<vector<double>> emissionSpectrum = {
      { // wavelengths
         373.422, 381.455, 387.619, 391.746, 394.052, 396.803, 399.751, 
         401.738, 402.733, 405.542, 406.956, 408.993, 411.443, 413.117, 
         415.796, 419.661, 431.188, 436.658, 449.632, 456.51, 467.612, 
         482.131, 494.627, 505.097
      },
      { // energies 
         3.32064*eV, 3.25072*eV, 3.19902*eV, 3.16532*eV, 
         3.14679*eV, 3.12498*eV, 3.10193*eV, 3.08659*eV, 3.07896*eV, 3.05764*eV, 3.04701*eV, 
         3.03184*eV, 3.01378*eV, 3.00157*eV, 2.98223*eV, 2.95477*eV, 2.87578*eV, 2.83975*eV, 
         2.75781*eV, 2.71626*eV, 2.65177*eV, 2.57191*eV, 2.50694*eV, 2.45497*eV
      },
      { // fast relative light output 
         0.0494898, 4.23782, 19.5498, 39.7241, 61.5267, 75.0436, 89.3671, 94.6063, 
         96.8218, 98.8231, 98.0055, 93.3453, 86.4601, 73.9236, 61.7843, 
         53.8795, 43.9029, 37.1994, 19.7383, 13.2273, 7.69804, 3.15583, 
         0.647436, 0.57693
      },
      { // slow relative light output
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
      },
      { // refractive index 
         1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 
         1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 
         1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58
      }
   };

   // Relative light output for different particles
   // refscint is the reference scintillation yield for the relative measuremnts.
   // It appears to be the scintillation yield from a 160 MeV proton.
   double refscint = fRefScint;//500.0/MeV;
   vector< vector<double>> relativeLightOutput = { 
      { // Particle kinetic energy
         0*MeV, 10*MeV, 20*MeV, 30*MeV, 40*MeV, 50*MeV, 50*MeV, 70*MeV, 90*MeV, 110*MeV, 130*MeV, 150*MeV, 170*MeV, 190*MeV, 190*MeV, 
         240*MeV, 290*MeV, 340*MeV, 390*MeV, 440*MeV, 490*MeV, 540*MeV, 590*MeV, 640*MeV, 690*MeV, 740*MeV, 790*MeV, 840*MeV, 
         890*MeV, 940*MeV, 990*MeV},
      { // 1 electron
         0.0010197*refscint, 0.070038*refscint, 0.139817*refscint, 0.210325*refscint, 0.281509*refscint, 
         0.353295*refscint, 0.353295*refscint, 0.498268*refscint, 0.644225*refscint, 0.789808*refscint, 0.932235*refscint, 1.07667*refscint,
         1.22111*refscint, 1.36555*refscint, 1.36555*refscint, 1.72665*refscint, 2.08774*refscint, 2.44884*refscint, 2.80994*refscint, 
         3.17103*refscint, 3.53213*refscint, 3.89323*refscint, 4.25432*refscint, 4.61542*refscint, 4.97652*refscint, 5.33761*refscint, 
         5.69871*refscint, 6.0598*refscint, 6.4209*refscint, 6.782*refscint, 7.14309*refscint},
      { // 2 proton
         0*refscint, 0.0458257*refscint, 0.0964081*refscint, 
         0.149766*refscint, 0.20579*refscint, 0.264325*refscint, 0.264325*refscint, 0.388084*refscint, 0.518921*refscint, 0.654003*refscint,
         0.786399*refscint, 0.919298*refscint, 1.0522*refscint, 1.1851*refscint, 1.1851*refscint, 1.51734*refscint, 1.84959*refscint, 
         2.18184*refscint, 2.51409*refscint, 2.84633*refscint, 3.17858*refscint, 3.51083*refscint, 3.84308*refscint, 4.17533*refscint, 
         4.50757*refscint, 4.83982*refscint, 5.17207*refscint, 5.50432*refscint, 5.83656*refscint, 6.16881*refscint, 
         6.50106*refscint},
      { // 3 deuteron
         0.000228012*refscint, 0.0357616*refscint, 0.0758948*refscint, 0.120514*refscint, 0.169428*refscint, 
         0.222373*refscint, 0.222373*refscint, 0.338909*refscint, 0.466477*refscint, 0.600215*refscint, 0.725919*refscint, 
         0.853917*refscint, 0.981914*refscint, 1.10991*refscint, 1.10991*refscint, 1.42991*refscint, 1.7499*refscint, 2.06989*refscint, 
         2.38989*refscint, 2.70988*refscint, 3.02988*refscint, 3.34987*refscint, 3.66987*refscint, 3.98986*refscint, 4.30986*refscint, 
         4.62985*refscint, 4.94984*refscint, 5.26984*refscint, 5.58983*refscint, 5.90983*refscint, 6.22982*refscint},
      { // 4 triton
         0.00327668*refscint, 0.0257885*refscint, 0.0518588*refscint, 0.081431*refscint, 0.114411*refscint, 0.150668*refscint, 0.150668*refscint, 
         0.232295*refscint, 0.324506*refscint, 0.424895*refscint, 0.530595*refscint, 0.636537*refscint, 0.742478*refscint, 0.84842*refscint,
         0.84842*refscint, 1.11327*refscint, 1.37813*refscint, 1.64298*refscint, 1.90784*refscint, 2.17269*refscint, 2.43754*refscint, 
         2.7024*refscint, 2.96725*refscint, 3.23211*refscint, 3.49696*refscint, 3.76181*refscint, 4.02667*refscint, 4.29152*refscint, 
         4.55638*refscint, 4.82123*refscint, 5.08609*refscint},
      { // 5 He3
         0.00217241*refscint, 0.0331461*refscint, 0.0682954*refscint, 
         0.107539*refscint, 0.150743*refscint, 0.197719*refscint, 0.197719*refscint, 0.301959*refscint, 0.417672*refscint, 
         0.541409*refscint, 0.667533*refscint, 0.792031*refscint, 0.916529*refscint, 1.04103*refscint, 1.04103*refscint, 1.35227*refscint, 
         1.66352*refscint, 1.97476*refscint, 2.28601*refscint, 2.59725*refscint, 2.9085*refscint, 3.21974*refscint, 3.53099*refscint, 
         3.84223*refscint, 4.15348*refscint, 4.46472*refscint, 4.77597*refscint, 5.08721*refscint, 5.39846*refscint, 5.7097*refscint, 
         6.02094*refscint},
      { // 6 alpha
         0.00438096*refscint, 0.0184309*refscint, 0.0354221*refscint, 0.0553226*refscint, 0.078079*refscint, 
         0.103617*refscint, 0.103617*refscint, 0.16263*refscint, 0.23134*refscint, 0.308381*refscint, 0.393657*refscint, 0.481043*refscint, 
         0.568428*refscint, 0.655813*refscint, 0.655813*refscint, 0.874276*refscint, 1.09274*refscint, 1.3112*refscint, 1.52967*refscint, 
         1.74813*refscint, 1.96659*refscint, 2.18506*refscint, 2.40352*refscint, 2.62198*refscint, 2.84045*refscint, 3.05891*refscint, 
         3.27737*refscint, 3.49584*refscint, 3.7143*refscint, 3.93276*refscint, 4.15123*refscint}
   };

   a = 1.01*g/mole;
   G4Element* elH  = new G4Element("Hydrogen", "H",  1.0, a);
   a = 12.01*g/mole;
   density   = 1.032*g/cm3;
   G4Element* elC  = new G4Element("Carbon"  , "C",  6.0, a);
   fScint1_mat = G4Material::GetMaterial("Scintillator",  /*print warning?*/ false);
   if( !fScint1_mat ) {
      fScint1_mat = new G4Material("Scintillator", density, 2);
      fScint1_mat->AddElement(elC, 9);
      fScint1_mat->AddElement(elH, 10);
   }

   G4MaterialPropertiesTable* Scnt_MPT = new G4MaterialPropertiesTable();
   Scnt_MPT->AddProperty("FASTCOMPONENT", emissionSpectrum[1].data(), emissionSpectrum[2].data(), emissionSpectrum[1].size());
   Scnt_MPT->AddProperty("SLOWCOMPONENT", emissionSpectrum[1].data(), emissionSpectrum[3].data(), emissionSpectrum[1].size());
   Scnt_MPT->AddProperty(     "RINDEX"  , emissionSpectrum[1].data(), emissionSpectrum[4].data(), emissionSpectrum[1].size());
   Scnt_MPT->AddConstProperty("SCINTILLATIONYIELD",  50./MeV  );
   Scnt_MPT->AddConstProperty("RESOLUTIONSCALE"   ,   5.0     );
   Scnt_MPT->AddConstProperty("FASTTIMECONSTANT"  ,   0.5*ns  );
   Scnt_MPT->AddConstProperty("SLOWTIMECONSTANT"  , 100.*ns   );
   Scnt_MPT->AddConstProperty("YIELDRATIO"        ,   1.0    );

   Scnt_MPT->AddProperty("ELECTRONSCINTILLATIONYIELD", relativeLightOutput[0].data(), relativeLightOutput[1].data(), relativeLightOutput[0].size());
   Scnt_MPT->AddProperty("PROTONSCINTILLATIONYIELD"  , relativeLightOutput[0].data(), relativeLightOutput[2].data(), relativeLightOutput[0].size());
   Scnt_MPT->AddProperty("DEUTERONSCINTILLATIONYIELD", relativeLightOutput[0].data(), relativeLightOutput[3].data(), relativeLightOutput[0].size());
   Scnt_MPT->AddProperty("TRITONSCINTILLATIONYIELD"  , relativeLightOutput[0].data(), relativeLightOutput[4].data(), relativeLightOutput[0].size());
   Scnt_MPT->AddProperty("HE3SCINTILLATIONYIELD"     , relativeLightOutput[0].data(), relativeLightOutput[5].data(), relativeLightOutput[0].size());
   Scnt_MPT->AddProperty("ALPHASCINTILLATIONYIELD"   , relativeLightOutput[0].data(), relativeLightOutput[6].data(), relativeLightOutput[0].size());

   // The relative strength of the fast component as a fraction of total scintillation yield is given by the YIELDRATIO. 
   if(fUseOpticalPhotons) {
      fScint1_mat->SetMaterialPropertiesTable(Scnt_MPT);
   }
   fScint2_mat = fScint1_mat;
}
//______________________________________________________________________________

void RecoilHodoDetectorGeometry3::BuildLogicalVolumes()
{
   using namespace CLHEP;

   BuildMaterials();

   bool    checkOverlaps    = false;
   int     natoms           = 0;
   int     ncomponents      = 0;
   double  A                = 0.0;
   double  Z                = 0.0;
   double  thickness        = 0.0;
   double  surface_z        = 0.0;
   double  red              = 0.0;
   double  green            = 0.0;
   double  blue             = 0.0;
   double  alpha            = 0.0;
   double  density          = 0.0;
   double  pressure         = 0.0;
   double  temperature      = 0.0;
   double  a                = 0.0;

   G4NistManager* nist = G4NistManager::Instance();

   // ------------------------------------------------------------------------

   fPhotonDet1_pos   = { fScintLength/2.0 - fPhotonDetThickness/2.0 ,0.0, 0.0 };
   fPhotonDet2_pos   = { -fScintLength/2.0 + fPhotonDetThickness/2.0,0.0, 0.0 };
   fPhotonDet3_pos   = { 0.0,0.0, fScint2Thickness/2.0 - fScintWrapThickness-fPhotonDetThickness/2.0};
   fPhotonDet4_pos   = fPhotonDet2_pos;

   // ------------------------------------------------------------------------
   // Container volume (air)
   // ------------------------------------------------------------------------
   red       = 6.0/256.0;
   green     = 200.0/256.0;
   blue      = 6.0/256.0;
   alpha     = 0.4;

   if(fContainer_log)   delete fContainer_log;
   if(fContainer_solid) delete fContainer_solid;

   fContainer_solid = new G4Tubs("fContainer_solid", fInnerRadius, 
         fInnerRadius + fScint1Thickness + fScint2Thickness + 2.0*cm, fScintLength/2.0,
         0.0, CLHEP::twopi);
   fContainer_mat   = nist->FindOrBuildMaterial("G4_AIR");
   fContainer_log   = new G4LogicalVolume(fContainer_solid, fContainer_mat,"fContainer_log");

   G4Colour            fContainer_color {red, green, blue, alpha };   // Gray 
   G4VisAttributes   * fContainer_vis   = new G4VisAttributes(fContainer_color);
   fContainer_vis->SetForceWireframe(true);
   fContainer_log->SetVisAttributes(fContainer_vis);

   auto user_stepping_action = const_cast<G4UserSteppingAction*>(G4RunManager::GetRunManager()->GetUserSteppingAction()) ;


   // ------------------------------------------------------------------------
   // First layer scintillator bar  
   // ------------------------------------------------------------------------
   /* clear
   red       = 256.0/256.0;
   green     = 1.0/256.0;
   blue      = 1.0/256.0;
   alpha     = 0.4;
   */
   
   red       = 220.0/256.0;
   green     = 23.0/256.0;
   blue      = 246.0/256.0;
   alpha     = 0.94;
  
   if(fScint1_log) delete fScint1_log;
   if(fScint1_solid) delete fScint1_solid;

   double dy1 = fInnerRadius*tan(fScintDeltaTheta/2.0) - fScintWrapThickness;
   double dy2 = dy1 + fScint1Thickness*tan(fScintDeltaTheta/2.0);

   fScint1_solid = new G4Trap("scint_keystone_solid",
         fScintLength/2.0, fScintLength/2.0,
         dy1             , dy2,
         fScint1Thickness/2.0 - fScintWrapThickness
         );

   fScint1_log   = new G4LogicalVolume(fScint1_solid, fScint1_mat,"fScint1_log");

   G4Colour            fScint1_color {red, green, blue, alpha }; 
   G4VisAttributes   * fScint1_vis   = new G4VisAttributes(fScint1_color);
   fScint1_vis->SetForceWireframe(true);
   fScint1_log->SetVisAttributes(fScint1_vis);

   G4UserLimits * fScint1_limits = new G4UserLimits(0.5*mm);
   fScint1_log->SetUserLimits(fScint1_limits);

   // ------------------------------------------------------------------------
   // Second layer scintillator bar  
   // ------------------------------------------------------------------------
   /*
   red       = 0.0/256.0;
   green     = 10.0/256.0;
   blue      = 200.0/256.0;
   alpha     = 0.4;
   */
   red       = 220.0/256.0;
   green     = 23.0/256.0;
   blue      = 246.0/256.0;
   alpha     = 0.54;
  
   if(fScint2_log) delete fScint2_log;
   if(fScint2_solid) delete fScint2_solid;

   double dyy1 = (fInnerRadius + fScint1Thickness )*tan(fScintDeltaTheta/2.0) - fScintWrapThickness;
   double dyy2 = dyy1 + fScint2Thickness*tan(fScintDeltaTheta/2.0);

   fScint2_solid = new G4Trap("scint_keystone_solid",
         fScint2Length/2.0-fScintWrapThickness,fScint2Length/2.0-fScintWrapThickness,// half-length along x at -dz and +dz
         dyy1, dyy2,                                // half-length along y at -dz and +dz
         fScint2Thickness/2.0 - fScintWrapThickness);
   fScint2_log   = new G4LogicalVolume(fScint2_solid, fScint2_mat,"fScint2_log");
   G4Colour            fScint2_color {red, green, blue, alpha };   // Gray 
   G4VisAttributes   * fScint2_vis   = new G4VisAttributes(fScint2_color);
   fScint2_vis->SetForceWireframe(true);
   fScint2_log->SetVisAttributes(fScint2_vis);

   G4UserLimits * fScint2_limits = new G4UserLimits(0.5*mm);
   fScint2_log->SetUserLimits(fScint2_limits);

   // ------------------------------------------------------------------------
   // First scintillator photon detection surface 
   // ------------------------------------------------------------------------
   /*
   red       = 250.0/256.0;
   green     = 0.0/256.0;
   blue      = 1.0/256.0;
   alpha     = 0.4;
   */
   red       = 247.0/256.0;
   green     = 68.0/256.0;
   blue      = 128.0/256.0;
   alpha     = 0.82;
   
   if(fPhotonDet1_phys)  delete fPhotonDet1_phys;
   if(fPhotonDet1_log)   delete fPhotonDet1_log;
   if(fPhotonDet1_solid) delete fPhotonDet1_solid;

   fPhotonDet1_mat   = fScint1_mat;
   fPhotonDet1_solid = new G4Trap("fPhotonDet1_solid",
         fPhotonDetThickness/2.0, fPhotonDetThickness/2.0, // half-length along x at -dz and +dz
         dy1-fScintWrapThickness, dy2-fScintWrapThickness, // half-length along y at -dz and +dz
         fScint1Thickness/2.0-fScintWrapThickness);
   fPhotonDet1_log   = new G4LogicalVolume(fPhotonDet1_solid, fPhotonDet1_mat ,"fPhotonDet1_log");
   fPhotonDet1_phys  = new G4PVPlacement(0,fPhotonDet1_pos, fPhotonDet1_log, "fPhotonDet1_phys",fScint1_log,true,0,checkOverlaps);
   G4Colour            fPhotonDet1_color {red, green, blue, alpha };   // Gray 
   G4VisAttributes   * fPhotonDet1_vis   = new G4VisAttributes(fPhotonDet1_color);
   fPhotonDet1_log->SetVisAttributes(fPhotonDet1_vis);
   // ------------------------------------------------------------------------

   if(fPhotonDet2_phys)  delete fPhotonDet2_phys;
   if(fPhotonDet2_log)   delete fPhotonDet2_log;
   if(fPhotonDet2_solid) delete fPhotonDet2_solid;

   fPhotonDet2_mat   = fScint1_mat;
   fPhotonDet2_solid = fPhotonDet1_solid;
   fPhotonDet2_log   = new G4LogicalVolume(fPhotonDet2_solid, fPhotonDet2_mat,"fPhotonDet2_log");
   fPhotonDet2_phys  = new G4PVPlacement(0,fPhotonDet2_pos, fPhotonDet2_log, "fPhotonDet2_phys",fScint1_log,true,0,checkOverlaps);
   fPhotonDet2_log->SetVisAttributes(fPhotonDet1_vis);


   // ------------------------------------------------------------------------
   // Second scintillator photon detection surface 
   // ------------------------------------------------------------------------
   // Tiled outer scintillator
   /*
   red       = 0.0/256.0;
   green     = 10.0/256.0;
   blue      = 200.0/256.0;
   alpha     = 0.4;
   */
   red       = 247.0/256.0;
   green     = 98.0/256.0;
   blue      = 68.0/256.0;
   alpha     = 0.82;
  
   if(fPhotonDet3_phys)  delete fPhotonDet3_phys;
   if(fPhotonDet3_log)   delete fPhotonDet3_log;
   if(fPhotonDet3_solid) delete fPhotonDet3_solid;

   fPhotonDet3_mat   = fScint2_mat;
   fPhotonDet3_solid = new G4Box("fPhotonDet3_solid",
         1.0*cm/2.0,//fScint2Length/2.0 - fScintWrapThickness/2.0,
         1.0*cm/2.0,//dyy2              - fScintWrapThickness/2.0,
         fPhotonDetThickness/2.0);

   fPhotonDet3_log   = new G4LogicalVolume(fPhotonDet3_solid, fPhotonDet3_mat ,"fPhotonDet3_log");
   fPhotonDet3_phys  = new G4PVPlacement(0,fPhotonDet3_pos, fPhotonDet3_log, "fPhotonDet3_phys",fScint2_log,true,0,checkOverlaps);
   G4Colour            fPhotonDet3_color {red, green, blue, alpha };   // Gray 
   G4VisAttributes   * fPhotonDet3_vis   = new G4VisAttributes(fPhotonDet3_color);
   fPhotonDet3_log->SetVisAttributes(fPhotonDet3_vis);

}
//______________________________________________________________________________

G4VPhysicalVolume * RecoilHodoDetectorGeometry3::PlacePhysicalVolume(
      G4LogicalVolume * mother, 
      G4VPhysicalVolume * mother_phys,
      G4VPhysicalVolume * adjacent_phys  )
{
   bool checkOverlaps = false;

   // scintillator surface
   const    G4int NUM           = 2;
   G4double pp2[NUM]            = {1.6*eV, 8.0*eV};
   G4double specularlobe2[NUM]  = {0.3, 0.3};
   G4double specularspike2[NUM] = {0.2, 0.2};
   G4double backscatter2[NUM]   = {0.1, 0.1};
   G4double rindex2[NUM]        = {1.35, 1.40};
   G4double reflectivity2[NUM]  = {0.9, 0.9};
   G4double efficiency2[NUM]    = {0.0, 0.0};

   G4OpticalSurface * OpSurface = new G4OpticalSurface("scintillator fiber");
   OpSurface->SetType(dielectric_metal);
   OpSurface->SetFinish(ground);
   OpSurface->SetModel(glisur);
   G4double polish = 0.8;

   G4MaterialPropertiesTable *OpSurfaceProperty = new G4MaterialPropertiesTable();
   OpSurfaceProperty->AddProperty("REFLECTIVITY",pp2,reflectivity2,NUM);
   OpSurfaceProperty->AddProperty("EFFICIENCY",  pp2,efficiency2,NUM);
   OpSurface->SetMaterialPropertiesTable(OpSurfaceProperty);

   G4RotationMatrix scint_rot = G4RotationMatrix::IDENTITY;
   scint_rot.rotateY(1.0*CLHEP::pi/2.0);
   G4ThreeVector s1_pos = {fScint1Thickness/2.0 + fInnerRadius,0.0,0.0};
   G4ThreeVector s2_pos = {fScint2Thickness/2.0 + fScint1Thickness + fInnerRadius,0.0,0.0};

   double delta_phi = 360.0*CLHEP::degree/double(fScint1_positions.size());
   G4String name       = "";
   int n_bars = fScint1_positions.size();


   // -----------------------------------------------------------
   name = "Hodoscope3Container";
   G4VPhysicalVolume * phys = new G4PVPlacement( 0, G4ThreeVector(0,0,0), 
            fContainer_log, name, mother, false, 0, checkOverlaps); 

   // -----------------------------------------------------------

   for(int i = 0; i<fScint1_positions.size(); i++) {

      s1_pos.rotateZ(delta_phi);
      s2_pos.rotateZ(delta_phi);
      scint_rot.rotateZ(delta_phi);

      fScint1_positions[i] = fScint1_pos + s1_pos;
      
      name                 = "fScint1_physicals_" + std::to_string(i);
      fScint1_physicals[i] = new G4PVPlacement(
            G4Transform3D(scint_rot,fScint1_positions[i]), 
            fScint1_log, name, fContainer_log, false, i, checkOverlaps); 
      name = "fScint1_borders_" + std::to_string(i);
      fScint1_borders[i]  = new G4LogicalBorderSurface(name, fScint1_physicals[i], phys, OpSurface );

      for(int j = 0; j<10; j++) {
         fScint2_positions[j][i] = fScint2_pos + s2_pos;
         G4ThreeVector Z_offset = {0,0,(30.0*cm/10.0)*(double(j-5)+0.5)};
         fScint2_positions[j][i] += Z_offset;
         name                 = "fScint2_physicals_" + std::to_string(i) + "_" + std::to_string(j);
         fScint2_physicals[j][i] = new G4PVPlacement(
               G4Transform3D(scint_rot,fScint2_positions[j][i]), 
               fScint2_log, name, fContainer_log, false, i*10+j, checkOverlaps); 
         name = "fScint2_borders_" + std::to_string(i);
         fScint2_borders[j][i]  = new G4LogicalBorderSurface(name, fScint2_physicals[j][i], phys, OpSurface );
      }
   }

   return phys;
}
//______________________________________________________________________________

