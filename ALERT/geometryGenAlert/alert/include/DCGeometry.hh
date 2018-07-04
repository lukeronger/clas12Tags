#ifndef __DCGEO_H
#define __DCGEO_H

#include "TMath.h"

#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Vector/Rotation.h"

#include <array>
#include <map>
#include <tuple>
#include <unordered_map>

namespace clas12 {

  namespace geo {

    using namespace CLHEP;
    using namespace TMath;

    const double  inch = 2.54*cm;

    // rotate a two vector
    Hep2Vector ApplyRotation(Hep2Vector v, double angle);

    namespace DC {

      // The parameters below can be found in:
      // https://clasweb.jlab.org/wiki/images/d/d1/Dc12geom_doc.pdf 
      // all extra parameters come from the drawings found at 
      // https://userweb.jlab.org/~cwiggins/

      const std::array<double,3> FrontGap      = {2.500*cm, 2.5000*cm, 2.5*cm};
      const std::array<double,3> MidGap        = {2.500*cm, 7.1613*cm, 2.5*cm};
      const std::array<double,3> BackGap       = {2.500*cm, 2.5000*cm, 2.5*cm};

      // This is the angle as viewed look down the beampipe
      const std::array<double,3> RegionAngle = {59.0*degree, 60.0*degree, 59.0*degree};

      const std::array<double,3> RegionTilt    = {25.0*degree, 25.0*degree, 25.0*degree};

      // 
      const std::array<double,3> RegionXOffset = {
	8.298*cm, 20.296*cm, 18.054*cm
      };

      // The construction of each chamber is done with the wires streched along the x axis
      // as opposed the SCS which has them along the y axis. There for we rotate the first
      // sector by -90 to bring it back to the SCS definiton.
      const std::array<double,6> SectorRotation = {
	-90.0*degree, -30.0*degree, 30.0*degree,
	90.0*degree, 150.0*degree, -150.0*degree
      };

      // Some (maybe) handy arrays
      const std::array<int,6> SuperLayerRegion      = { 1,1,2,2,3,3 };
      const std::array<int,6> SuperLayerRegionIndex = { 0,0,1,1,2,2 };
      const std::array<std::array<int,2>,3> RegionSuperLayers = { 
	std::array<int,2>{1,2}, 
	std::array<int,2>{3,4}, 
	std::array<int,2>{5,6} };
      const std::array<std::array<int,2>,3> RegionSuperLayerIndex = { 
	std::array<int,2>{0,1}, 
	std::array<int,2>{2,3}, 
	std::array<int,2>{4,5} };

      // The stereo angle for each superlayer
      const std::array<double,6> ThetaStereo = {
	6.0*degree, -6.0*degree,
	6.0*degree, -6.0*degree,
	6.0*degree, -6.0*degree
      };

      // polar angle to the first guard wires's mid-plane point
      const std::array<double,6> ThetaMin = {
	4.694*degree, 4.495*degree,
	4.812*degree, 4.771*degree,
	4.333*degree, 4.333*degree
      };

      // This is the wire separation (for each superlayer)
      // Each side of the hexagon forming the wire volume is of length 2*LayerSep
      const std::array<double,6> LayerSep = {
	0.3861*cm, 0.4042*cm,
	0.6219*cm, 0.6586*cm,
	0.9351*cm, 0.9780*cm
      };

      // This is the spacing between wires that are in the same wire plane. This 
      // measured separation is perpendicular to LayerSep.
      const std::array<double,6> LayerWireSpacing = {
	2.0*TMath::Sqrt(3.0)*LayerSep.at(0), 2.0*TMath::Sqrt(3.0)*LayerSep.at(1),
	2.0*TMath::Sqrt(3.0)*LayerSep.at(2), 2.0*TMath::Sqrt(3.0)*LayerSep.at(3),
	2.0*TMath::Sqrt(3.0)*LayerSep.at(4), 2.0*TMath::Sqrt(3.0)*LayerSep.at(5)
      };

      // 
      const std::array<double,6> LayerWirePlaneLength = {
	113.0*LayerWireSpacing.at(0), 113.0*LayerWireSpacing.at(1),
	113.0*LayerWireSpacing.at(2), 113.0*LayerWireSpacing.at(3),
	113.0*LayerWireSpacing.at(4), 113.0*LayerWireSpacing.at(5)
      };

      // Distance from target to first wire plane in super layer
      // This we call center 
      const std::array<double,6> DistanceToTarget = {
	228.078*cm, 238.687*cm,
	351.544*cm, 371.773*cm,
	489.099*cm, 511.236*cm
      };

      const std::array<double,6> SuperLayerXOffset = {
	8.298*cm,  8.298*cm, 
	14.296*cm, 14.296*cm,
	18.054*cm, 18.054*cm
      };

      // Width of superlayer from front to back gaurd wires  
      const std::array<double,6> SuperLayerWidth = {
	21.0*LayerSep.at(0), 21.0*LayerSep.at(1), 
	21.0*LayerSep.at(2), 21.0*LayerSep.at(3),
	21.0*LayerSep.at(4), 21.0*LayerSep.at(5)
      };

      // offset of first superlayer from front window for given region.
      const std::array<double,6> SuperLayerZOffset = {
	FrontGap.at(0) + 1.5*LayerSep.at(0), 
	FrontGap.at(0) + SuperLayerWidth.at(0) + MidGap.at(0)+1.5*LayerSep.at(1),
	FrontGap.at(1) + 3.5*LayerSep.at(2), 
	FrontGap.at(1) + SuperLayerWidth.at(2) + MidGap.at(1)+1.5*LayerSep.at(3),
	FrontGap.at(2) + 1.5*LayerSep.at(4), 
	FrontGap.at(2) + SuperLayerWidth.at(4) + MidGap.at(2)+1.5*LayerSep.at(5)
      };


      // Distance from target to first guard wire (at ThetaMin)
      const std::array<double,6> CenterToRefWire = {
	DistanceToTarget.at(0)*Tan( RegionTilt.at(0) - ThetaMin.at(0) ),
	DistanceToTarget.at(1)*Tan( RegionTilt.at(0) - ThetaMin.at(1) ),
	DistanceToTarget.at(2)*Tan( RegionTilt.at(1) - ThetaMin.at(2) ),
	DistanceToTarget.at(3)*Tan( RegionTilt.at(1) - ThetaMin.at(3) ),
	DistanceToTarget.at(4)*Tan( RegionTilt.at(2) - ThetaMin.at(4) ),
	DistanceToTarget.at(5)*Tan( RegionTilt.at(2) - ThetaMin.at(5) )
      };

      // Distance from first guard wire (at ThetaMin) to target 
      const std::array<double,6> DistanceToRefWire = {
	Sqrt( Power(DistanceToTarget.at(0),2.0) +  Power(CenterToRefWire.at(0),2.0) ),
	Sqrt( Power(DistanceToTarget.at(1),2.0) +  Power(CenterToRefWire.at(1),2.0) ),
	Sqrt( Power(DistanceToTarget.at(2),2.0) +  Power(CenterToRefWire.at(2),2.0) ),
	Sqrt( Power(DistanceToTarget.at(3),2.0) +  Power(CenterToRefWire.at(3),2.0) ),
	Sqrt( Power(DistanceToTarget.at(4),2.0) +  Power(CenterToRefWire.at(4),2.0) ),
	Sqrt( Power(DistanceToTarget.at(5),2.0) +  Power(CenterToRefWire.at(5),2.0) )
      };

      // Offset in the end-plate coordinates (tilted)
      const std::array<Hep2Vector,6> RefWireOffset = {
	Hep2Vector( Cos(RegionTilt.at(0))*(DistanceToRefWire.at(0)*Sin( ThetaMin.at(0) ) - SuperLayerXOffset.at(0)),
		    Sin(RegionTilt.at(0))*(DistanceToRefWire.at(0)*Sin( ThetaMin.at(0) ) - SuperLayerXOffset.at(0)) ),
	Hep2Vector( Cos(RegionTilt.at(0))*(DistanceToRefWire.at(1)*Sin( ThetaMin.at(1) ) - SuperLayerXOffset.at(1)),
		    Sin(RegionTilt.at(0))*(DistanceToRefWire.at(1)*Sin( ThetaMin.at(1) ) - SuperLayerXOffset.at(1)) ),
	Hep2Vector( Cos(RegionTilt.at(1))*(DistanceToRefWire.at(2)*Sin( ThetaMin.at(2) ) - SuperLayerXOffset.at(2)),
		    Sin(RegionTilt.at(1))*(DistanceToRefWire.at(2)*Sin( ThetaMin.at(2) ) - SuperLayerXOffset.at(2)) ),
	Hep2Vector( Cos(RegionTilt.at(1))*(DistanceToRefWire.at(3)*Sin( ThetaMin.at(3) ) - SuperLayerXOffset.at(3)),
		    Sin(RegionTilt.at(1))*(DistanceToRefWire.at(3)*Sin( ThetaMin.at(3) ) - SuperLayerXOffset.at(3)) ),
	Hep2Vector( Cos(RegionTilt.at(2))*(DistanceToRefWire.at(4)*Sin( ThetaMin.at(4) ) - SuperLayerXOffset.at(4)),
		    Sin(RegionTilt.at(2))*(DistanceToRefWire.at(4)*Sin( ThetaMin.at(4) ) - SuperLayerXOffset.at(4)) ),
	Hep2Vector( Cos(RegionTilt.at(2))*(DistanceToRefWire.at(5)*Sin( ThetaMin.at(5) ) - SuperLayerXOffset.at(5)),
		    Sin(RegionTilt.at(2))*(DistanceToRefWire.at(5)*Sin( ThetaMin.at(5) ) - SuperLayerXOffset.at(5)) )
      };

      const std::array<double,6> LeftFeedThruOffset = {
	0.4745*cm, 0.5254*cm,
	0.0282*cm, 0.0789*cm,
	0.4745*cm, 0.5254*cm
      };

      const std::array<double,6> RightFeedThruOffset = {
	0.5254*cm, 0.4745*cm,
	0.0789*cm, 0.0282*cm,
	0.5254*cm, 0.4745*cm
      };
      //_____________________________________________________________________
      // Nose piece dimensions (width inside gas volue)
      // and width outside the gas volume (nearest to beam)
      const std::array<double,3> NoseWidthInside  = { 
	4.25*inch, 
	5.47*inch, // at it largest point. Need to cut out SL areas
	20.7467*cm
      };
      const std::array<double,3> NoseWidthOutside = {
	1.987*inch,
	1.772*inch,
	12.1260*cm
      };

      // the depth along z in the untilted frame
      const std::array<double,3> NoseDepth = {
	(12.405 - 2.0*Tan(RegionTilt.at(0)) + (1.547*Tan(RegionTilt.at(0))-0.311) )*2.54*cm,
	21.152*inch,
	23.419*inch  // from drawing 01-10-0201
      };
      const std::array<double,3> NoseHeight = {
	2.0*2.54*cm,
	3.156*inch,
	3.0*2.54*cm
      };
      const std::array<double,3> NoseTiltShift = {
	NoseHeight.at(0)*Tan(RegionTilt.at(0)), 
	NoseHeight.at(1)*Tan(RegionTilt.at(1)), 
	NoseHeight.at(2)*Tan(RegionTilt.at(2))
      };
      // Opening shift is the x length looking down beam pipe
      const std::array<double,3> NoseOpeningShift = {
	NoseHeight.at(0)*Tan(RegionAngle.at(0)/2.0), 
	NoseHeight.at(1)*Tan(RegionAngle.at(1)/2.0), 
	NoseHeight.at(2)*Tan(RegionAngle.at(2)/2.0)
      };
      const std::array<double,3> NoseSideLength = {
	Sqrt(NoseHeight.at(0)*NoseHeight.at(0)+ 
	     NoseTiltShift.at(0)*NoseTiltShift.at(0)+
	     NoseOpeningShift.at(0)*NoseOpeningShift.at(0) ),
	Sqrt(NoseHeight.at(1)*NoseHeight.at(1)+ 
	     NoseTiltShift.at(1)*NoseTiltShift.at(1)+
	     NoseOpeningShift.at(1)*NoseOpeningShift.at(1) ),
	Sqrt(NoseHeight.at(2)*NoseHeight.at(2)+
	     NoseTiltShift.at(2)*NoseTiltShift.at(2)+
	     NoseOpeningShift.at(2)*NoseOpeningShift.at(2) )
      };
      const std::array<double,3> NoseEndPlateAngle = {
	90*degree-(Hep3Vector(NoseTiltShift.at(0), NoseOpeningShift.at(0), NoseHeight.at(0)).angle(
												   Hep3Vector(0.0, NoseOpeningShift.at(0), NoseHeight.at(0))) ),
	90*degree-(Hep3Vector(NoseTiltShift.at(1), NoseOpeningShift.at(1), NoseHeight.at(1)).angle(
												   Hep3Vector(0.0, NoseOpeningShift.at(1), NoseHeight.at(1))) ),
	90*degree-(Hep3Vector(NoseTiltShift.at(2), NoseOpeningShift.at(2), NoseHeight.at(2)).angle(
												   Hep3Vector(0.0, NoseOpeningShift.at(2), NoseHeight.at(2))) )
      };
      const std::array<double,3> NoseToEndPlateRotationAngle = {
	-(Hep3Vector(NoseTiltShift.at(0), NoseOpeningShift.at(0), NoseHeight.at(0)).angle(
											  Hep3Vector(0.0, NoseOpeningShift.at(0), NoseHeight.at(0))) ),
	-(Hep3Vector(NoseTiltShift.at(1), NoseOpeningShift.at(1), NoseHeight.at(1)).angle(
											  Hep3Vector(0.0, NoseOpeningShift.at(1), NoseHeight.at(1))) ),
	-(Hep3Vector(NoseTiltShift.at(2), NoseOpeningShift.at(2), NoseHeight.at(2)).angle(
											  Hep3Vector(0.0, NoseOpeningShift.at(2), NoseHeight.at(2))) )
      };

      // Ref Locations (below) in the non-tilted nose-plate plane (as in 01-10-0201)
      const std::array<double,3> NoseRefLocation_Y = {
	0.375*2.54*cm,
	NoseSideLength.at(1)*Sin(NoseEndPlateAngle.at(1)) - (2.265*2.54*cm),
	NoseSideLength.at(2)*Sin(NoseEndPlateAngle.at(2)) - (2.265*2.54*cm)
      };
      const std::array<double,3> NoseRefLocation_X = {
	NoseDepth.at(0) - (22.586*2.54*cm),
	NoseDepth.at(1) - (22.586*2.54*cm),
	NoseDepth.at(2) - (22.586*2.54*cm)
      };
      const std::array<Hep2Vector,3> NoseRefLocation = {
	Hep2Vector(NoseRefLocation_X.at(0), NoseRefLocation_Y.at(0)),
	Hep2Vector(NoseRefLocation_X.at(1), NoseRefLocation_Y.at(1)),
	Hep2Vector(NoseRefLocation_X.at(2), NoseRefLocation_Y.at(2))
      };
      //_____________________________________________________________________

      // Now working in x-y plane in tilted endplate coords
      // y is along the wire plane and x is normal to the front window
      //
      // The reference location for the nose to end plate is the closest upstream and to the beam
      // of the three screw locations for the endplate to mount to the nose plate.
      // Note that  the hole placement on the plate side of the 3 screw holes 
      // is asymmetric in drawing 01-10-0201 and 01-08-2001.
      // The end plate coordinate system for the reference point is in the tilted system
      // with the orgin at the nose/plate corner upstream and closest to the beamline 
      //const std::array<double,3> NoseEndPlateRefLocation = {
      //   NoseSideLength 2.265*2.54*
      //};
      const std::array<Hep2Vector,3> NoseEndPlateRefLocation = {
	ApplyRotation(Hep2Vector(NoseRefLocation_X.at(0), NoseRefLocation_Y.at(0)),NoseToEndPlateRotationAngle.at(0)),
	ApplyRotation(Hep2Vector(NoseRefLocation_X.at(1), NoseRefLocation_Y.at(1)),NoseToEndPlateRotationAngle.at(1)),
	ApplyRotation(Hep2Vector(NoseRefLocation_X.at(2), NoseRefLocation_Y.at(2)),NoseToEndPlateRotationAngle.at(2))
      };


      // 
      const std::array<Hep2Vector,3> EndPlateRefLocation_D = {
	NoseEndPlateRefLocation.at(0) + Hep2Vector( -0.280*2.54*cm, 1.655*2.54*cm),
	NoseEndPlateRefLocation.at(1) + Hep2Vector( -0.280*2.54*cm, 1.655*2.54*cm),
	NoseEndPlateRefLocation.at(2) + Hep2Vector( -0.280*2.54*cm, 1.655*2.54*cm)
      };

      const std::array<Hep2Vector,6> EndPlateSLRefWire_Left = {
	EndPlateRefLocation_D.at(0) + Hep2Vector(  1.447*2.54*cm,  0.453*2.54*cm),
	EndPlateRefLocation_D.at(0) + Hep2Vector(  5.710*2.54*cm, -0.588*2.54*cm),
	EndPlateRefLocation_D.at(1) + Hep2Vector(  1.447*2.54*cm,  0.453*2.54*cm),
	EndPlateRefLocation_D.at(1) + Hep2Vector(  5.710*2.54*cm, -0.588*2.54*cm),
	EndPlateRefLocation_D.at(2) + Hep2Vector(  1.447*2.54*cm,  0.453*2.54*cm),
	EndPlateRefLocation_D.at(2) + Hep2Vector(  5.710*2.54*cm, -0.588*2.54*cm)
      };

      const std::array<Hep2Vector,6> EndPlateSLRefWire_Right = {
	EndPlateRefLocation_D.at(0) + Hep2Vector(  1.440*2.54*cm,  1.079*2.54*cm),
	EndPlateRefLocation_D.at(0) + Hep2Vector(  5.717*2.54*cm, -1.221*2.54*cm),
	EndPlateRefLocation_D.at(1) + Hep2Vector(  1.440*2.54*cm,  1.079*2.54*cm),
	EndPlateRefLocation_D.at(1) + Hep2Vector(  5.717*2.54*cm, -1.221*2.54*cm),
	EndPlateRefLocation_D.at(2) + Hep2Vector(  1.440*2.54*cm,  1.079*2.54*cm),
	EndPlateRefLocation_D.at(2) + Hep2Vector(  5.717*2.54*cm, -1.221*2.54*cm)
      };

      // Wire spacing in a layer
      const std::array<double,6> EndPlateWireSpacing_Left = {
	0.565*2.54*cm, 0.658*2.54*cm,
	0.600*2.54*cm, 0.600*2.54*cm,
	0.600*2.54*cm, 0.600*2.54*cm
      };

      const std::array<double,6> EndPlateWireSpacing_Right = {
	0.629*2.54*cm, 0.591*2.54*cm,
	0.700*2.54*cm, 0.700*2.54*cm,
	0.700*2.54*cm, 0.700*2.54*cm
      };

      // Spacing 
      const std::array<double,6> EndPlateLayerSpacing_Left = {
	0.155*2.54*cm, 0.163*2.54*cm,
	0.100*2.54*cm, 0.100*2.54*cm,
	0.100*2.54*cm, 0.100*2.54*cm
      }; 
      const std::array<double,6> EndPlateLayerSpacing_Right = {
	0.155*2.54*cm, 0.163*2.54*cm,
	0.100*2.54*cm, 0.100*2.54*cm,
	0.100*2.54*cm, 0.100*2.54*cm
      };

      // vector to a wire in the next wire layer with the wire being above the current(in +y)
      const std::array<Hep2Vector,6> EndPlateLayerTranslate1_Right = {
	Hep2Vector(0.155*2.54*cm, 0.336*2.54*cm),
	Hep2Vector(0.163*2.54*cm, 0.309*2.54*cm),
	Hep2Vector(0.155*2.54*cm, 0.336*2.54*cm),
	Hep2Vector(0.163*2.54*cm, 0.309*2.54*cm),
	Hep2Vector(0.155*2.54*cm, 0.336*2.54*cm),
	Hep2Vector(0.163*2.54*cm, 0.309*2.54*cm)
      };
      const std::array<Hep2Vector,6> EndPlateLayerTranslate2_Right = {
	Hep2Vector(2.0*0.155*2.54*cm, 0.042*2.54*cm),
	Hep2Vector(2.0*0.163*2.54*cm, 0.026*2.54*cm),
	Hep2Vector(2.0*0.155*2.54*cm, 0.042*2.54*cm),
	Hep2Vector(2.0*0.163*2.54*cm, 0.042*2.54*cm),
	Hep2Vector(2.0*0.155*2.54*cm, 0.042*2.54*cm),
	Hep2Vector(2.0*0.163*2.54*cm, 0.042*2.54*cm)
      };

      const std::array<double,3> EndPlateLongSideLength  = {
	81.305*inch,
	(125.468+8.00)*inch,
	191.275*inch
      };
      const std::array<double,3> EndPlateShortSideLength = {
	(80.130-4.468)*inch,
	(123.432+0.551)*inch,
	(187.59+2.308)*inch // Here we don't cut the corner at 20 degrees
      };
      const std::array<double,3> EndPlateWidth  = {
	11.01*inch,
	(18.102+0.335)*inch,
	21.78*inch
      };
      const std::array<double,3> EndPlateTiltAngle  = {
	ATan(4.468/11.01),
	ATan((8.0-0.551)/(18.102+0.335)),
	ATan(8.839/21.78)
      };
      const std::array<double,3> EndPlateThickness  = {
	0.375*inch,
	0.432*inch,
	0.375*inch
      };

      // Extra gap between endplate and container volume
      const std::array<double,3> EndPlateExtraGap  = {
	0.3*inch,
	1.0*inch,
	1.0*inch
      };
      // top/Back plate extra gap between the gas and container volumes
      const std::array<double,3> BackPlateExtraGap  = {
	2.0*inch,
	8.0*inch,
	10.0*inch
      };
      // 
      const std::array<double,3> FrontBackExtraGap  = {
	3.0*inch,
	5.0*inch,
	8.0*inch
      };

    }
    //________________________________________________________________________

    std::vector<Hep2Vector> ContainerTrapPoints(int region);
    double                  ContainerTrapWidth(int region);

    // This is the shift in the y direciont added to the continer so that 
    // there is a gap (EndPlateExtraGap) on the side.
    double                  ContainerExtraShift(int region);

    std::vector<Hep2Vector> RegionTrapPoints(int region);
    double                  RegionTrapWidth(int region);

    double                  RegionTrapOffset(int region);
    double                  RegionTrapCorner_z(int region);
    double                  RegionTrapCorner_y(int region);
    HepRotation             RegionRotation(int sec, int region);
    Hep3Vector              RegionTranslation(int sec, int region);

    double                  SuperLayerRefWire_z(int sl);
    double                  SuperLayerRefWire_y(int sl);
    HepRotation             LayerWireRotation(int sl);

    Hep3Vector              ToWireMidPlane(int sl, int layer, int wire);
    double                  WireLength(int sl,int layer,int wire);
    double                  WireShift(int sl,int layer,int wire);
    Hep3Vector              WireStereoShift(int sl, int layer, int wire);

    HepRotation             GetLeftEndplateRotation(int region);
    Hep3Vector              GetLeftEndplatePosition( int sec, int region);
    HepRotation             GetRightEndplateRotation(int region);
    Hep3Vector              GetRightEndplatePosition(int sec, int region);
    std::vector<Hep2Vector> EndplateTrapPoints(int region);

    HepRotation             SectorZRotation(int sec);

    // th is the sector angle (eg 59 deg)
    double  WindowAngle(double th, double th_tilt=25*CLHEP::degree);

    // Total depth of region from front to back windows
    // (entering normal to window, ie in the tilted coords) 
    double  RegionDepth(int region);

  }
}

#endif

