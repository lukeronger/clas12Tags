
#ifndef clas12_geo_RH_RCGeometry_H
#define clas12_geo_RH_RCGeometry_H

#include "DCWire.hh"

#include "TMath.h"

#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Geometry/Transform3D.h"

#include <array>
#include <map>
#include <tuple>
#include <unordered_map>

namespace clas12 {

  namespace wire {

    using namespace CLHEP;
    using namespace TMath;
    //using namespace HepGeom;

    /*! Recoil Chamber Geometry.
      Numbering conventions:
      - Each "Layer" has NCells number of unit cells
      - Each  unit cell has 1 sense wire and 8 field wires
      - Each  unit cell has 4 sub-cells
      Sub-cell numbering:
      2_______3 
      \__|__/
      \_|_/
      1   0 
      - These sub cells are given the same copy number in geant4
      Field wire number for cell:
      4___5___6 
      \__|__/
      3 \_|_/ 7
      2 1 0 
      - The Field wires surrounding the wire are number going around the wire starting with the inner radius at smallest phi values 
      - Each unit cell has a unique "channel number" starting at 0.
      - Each sub-cell trapeziod has 8 points defining it
      upstream:  2_____3     downstream:  6_____7
      \   /                    \   /
      \_/                      \_/
      1 0                      5 4 
                  

    */
    class RCGeometry {
    public:
      double InnerRadius    = 30.000*mm;
      double OuterRadius    = 79.995*mm;
      double Length         = 30.000*cm;
      double DeltaR         = 2.0*mm;
      int    NLayers        = 8;
      int    NCellLayers    = 3;
      double MaxWireSep     = 2.0*mm; // *desired* wire separation around the circumference. 
      double WireLength     = Length/Cos(10.0*degree);
      double StereoAngle    = 10.0*degree;

      double InnerGap       = 1.0*mm;
      double LayerGap       = 1.4*mm;

      std::array<int   , 8>   fNCells;
      std::array<int   , 8>   fFirstChannelInLayer;
      std::array<double, 8>   fCellInnerRadius;
      std::array<double, 8>   fCellDeltaPhi;
      std::array<double, 8>   fCellStereoAngle;
      std::array<double, 8>   fCellCentralRadius;

      double fNChannels;

    public:
      RCGeometry();
      void Print();

      int        GetWireLayer(int channel) const ;

      double     GetSenseWireRadius(int channel) const ;
      double     GetSenseWirePhi(int channel) const ;
      double     GetSenseWireStereoAngle(int channel) const ;

      Hep3Vector  GetFirstCellPosition(int layer, int subcell=-1) const ;
      HepRotation GetFirstCellRotation(int layer, int subcell=-1) const ;
      HepGeom::Transform3D GetFirstCellTransform(int layer, int subcell=-1) const ;

      Hep3Vector  GetSenseWirePerp(int channel) const ;

      Hep3Vector  GetSenseWirePosition(int channel) const ;
      HepRotation GetSenseWireRotation(int channel) const ;
      HepGeom::Transform3D GetSenseWireTransform(int channel) const ;

      Hep3Vector GetUpstreamSenseWirePosition(int channel) const ;
      Hep3Vector GetDownstreamSenseWirePosition(int channel) const ;

      Hep3Vector GetUnitSubCellPoint(int layer,int subcell, int point) const ;

      std::vector<Hep2Vector> GetSubCellTrapPoints(int layer, int subcell) ;

    };

  }
}

#endif

