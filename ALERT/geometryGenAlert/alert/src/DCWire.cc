#include "DCGeometry.hh"
#include "TMath.h"
#include <iomanip>


namespace clas12 {

  namespace  wire {

    using namespace CLHEP;
    using namespace TMath;
    using namespace clas12::geo;
    
    Hep2Vector ApplyRotation(Hep2Vector v, double angle)
    {
      v.rotate(angle);
      return v;
    }
    //______________________________________________________________________________

    std::vector<Hep2Vector> ContainerTrapPoints(int region)
    {
      // returns the 8 points defining a trapezoid that contains all parts of a 
      // given region's drift chamber. This includes the 
      // endplates, the front and back windows, and noseplate..
      // The first point is "nose" point were the endplates meet 
      // closest to the z axis and furthest upstream towards the target.
      // The origin of the solid's  coordinate system is this first gaurd wire
      // in the midplane (used to define theta_min).
      using namespace clas12::geo::DC;

      int index = region-1;
      std::vector<Hep2Vector> points;
      if( index < 0 || index >= 3 )
	{
	  std::cout << "clas12::geo::RegionTrapPoints(region=" << region << ") : Region out of range\n";
	  return  points;
	}
      double RegionLength = RegionDepth(region) + FrontBackExtraGap.at(region-1);

      // The mid plane cross section is a rectangle plus a right triangle. 
      // extra_back is the length of the triangle on the back window plane
      //    ______
      //   /| 
      //  /_|_____ 
      //   e        
      // The line segment above e is extra_back
      double extra_back = RegionLength*TMath::Tan(RegionTilt[index]);

      // This returns the nose angle (between two end plates) for the window when viewed normal to it
      double actual_plate_angle = WindowAngle( RegionAngle[index]/2.0 , RegionTilt[index] ) ;

      // this is the shift along with the desired gap in the top end plate
      double db_plus_G = ContainerExtraShift(region) + BackPlateExtraGap.at(index);
      double delta_l   = db_plus_G/Cos(actual_plate_angle);

      double l1 = EndPlateShortSideLength[index] + delta_l;
      double l2 = EndPlateLongSideLength[index] + delta_l;

      // nose height cut from triangular tip.
      double nose_width  =  NoseWidthOutside[index]/2.0;
      double nose_height =  nose_width/TMath::Tan(actual_plate_angle);
      double nose_side   = TMath::Sqrt(nose_width*nose_width + nose_height*nose_height);

      double inside_nose_width  =  NoseWidthInside[index]/2.0;
      double inside_nose_height =  inside_nose_width/TMath::Tan(actual_plate_angle);
      double inside_nose_side   = TMath::Sqrt(inside_nose_width*inside_nose_width + inside_nose_height*inside_nose_height);

      double h1 = (l1+nose_side)*TMath::Cos(actual_plate_angle);
      double x1 = (l1+nose_side)*TMath::Sin(actual_plate_angle);
      double h2 = (l2+nose_side)*TMath::Cos(actual_plate_angle);
      double x2 = (l2+nose_side)*TMath::Sin(actual_plate_angle);

      Hep2Vector p0(  nose_width, nose_height );
      Hep2Vector p1( -nose_width, nose_height );
      Hep2Vector p2( -x1,  h1);
      Hep2Vector p3(  x1,  h1);

      double dx = extra_back;//RegionLength*TMath::Tan( RegionTilt[index]  );

      Hep2Vector p4(  nose_width, nose_height-dx );
      Hep2Vector p5( -nose_width, nose_height-dx );
      Hep2Vector p6( -x2, h2-dx);
      Hep2Vector p7(  x2, h2-dx);

      points = { p1, p2, p3, p0, p5, p6, p7, p4 };

      return points;
    }
    //______________________________________________________________________________

    double  ContainerTrapWidth(int region)
    {
      int index = region-1;
      if( index < 0 || index >= 3 )
	{
	  std::cout << "clas12::geo::RegionTrapWidth(region=" << region << ") : Region out of range\n";
	  return  0.0;
	}
      using namespace clas12::geo::DC;
      double width = RegionDepth(region);
      width += FrontBackExtraGap.at(region-1);
      return width/2.0;
    }
    //______________________________________________________________________________

    double  ContainerExtraShift(int region)
    {
      // This is the shift in the y direciont added to the continer so that 
      // there is a gap (EndPlateExtraGap) on the side.
      using namespace clas12::geo::DC;
      int index = region-1;
      double actual_plate_angle = WindowAngle( RegionAngle[index]/2.0 , RegionTilt[index] ) ;
      return( EndPlateExtraGap.at(index)/Sin(actual_plate_angle) );
    }
    //______________________________________________________________________________

    std::vector<Hep2Vector> RegionTrapPoints(int region)
    {
      // returns the 8 points defining a trapezoid volume defined by the two 
      // endplates and the front and back windows.
      // The first point is "nose" point were the endplates meet 
      // closest to the z axis and furthest upstream towards the target.
      // The origin of the solid's  coordinate system is this first gaurd wire
      // in the midplane (used to define theta_min).
      using namespace clas12::geo::DC;

      int index = region-1;
      std::vector<Hep2Vector> points;
      if( index < 0 || index >= 3 )
	{
	  std::cout << "clas12::geo::RegionTrapPoints(region=" << region << ") : Region out of range\n";
	  return  points;
	}
      double RegionLength = RegionDepth(region);

      // The mid plane cross section is a rectangle plus a right triangle. 
      // extra_back is the length of the triangle on the back window plane
      //    ______
      //   /| 
      //  /_|_____ 
      //   e        
      // The line segment above e is extra_back
      double extra_back = RegionLength*TMath::Tan(RegionTilt[index]);

      double l1 = EndPlateShortSideLength[index];
      double l2 = EndPlateLongSideLength[index];

      // This returns the nose angle (between two end plates) for the window when viewed normal to it
      double actual_plate_angle = WindowAngle( RegionAngle[index]/2.0 , RegionTilt[index] );

      //std::cout << actual_plate_angle/degree << " vs " << RegionAngle[index]/degree << std::endl;

      // nose height cut from triangular tip.
      double nose_width  = NoseWidthOutside[index]/2.0;
      double nose_height = nose_width/TMath::Tan(actual_plate_angle);
      double nose_side   = TMath::Sqrt(nose_width*nose_width + nose_height*nose_height);

      double inside_nose_width  =  NoseWidthInside[index]/2.0;
      double inside_nose_height =  inside_nose_width/TMath::Tan(actual_plate_angle);
      double inside_nose_side   = TMath::Sqrt(inside_nose_width*inside_nose_width + inside_nose_height*inside_nose_height);

      double h1 = (l1+inside_nose_side)*TMath::Cos(actual_plate_angle);
      double x1 = (l1+inside_nose_side)*TMath::Sin(actual_plate_angle);
      double h2 = (l2+inside_nose_side)*TMath::Cos(actual_plate_angle);
      double x2 = (l2+inside_nose_side)*TMath::Sin(actual_plate_angle);

      //std::cout << "Rough check : " << extra_back << " vs " << l2-l1 << std::endl;

      Hep2Vector p0(  inside_nose_width, inside_nose_height );
      Hep2Vector p1( -inside_nose_width, inside_nose_height );
      Hep2Vector p2( -x1,  h1);
      Hep2Vector p3(  x1,  h1);

      double dx = extra_back;//RegionLength*TMath::Tan( RegionTilt[index]  );

      Hep2Vector p4(  inside_nose_width, inside_nose_height-dx );
      Hep2Vector p5( -inside_nose_width, inside_nose_height-dx );
      Hep2Vector p6( -x2, h2-dx);
      Hep2Vector p7(  x2, h2-dx);

      points = { p1, p2, p3, p0, p5, p6, p7, p4 };

      return points;
    }
    //______________________________________________________________________________

    double  RegionTrapWidth(int region)
    {
      // returns the half length of the region's trapezoid
      int index = region-1;
      if( index < 0 || index >= 3 )
	{
	  std::cout << "clas12::geo::RegionTrapWidth(region=" << region << ") : Region out of range\n";
	  return  0.0;
	}
      using namespace clas12::geo::DC;
      double RegionLength = RegionDepth(region);
      return RegionLength/2.0;
    }
    //______________________________________________________________________________

    double  RegionTrapOffset(int region)
    {
      int index = region-1;
      if( index < 0 || index >= 3 )
	{
	  std::cout << "clas12::geo::RegionTrapOffset(region=" << region << ") : Region out of range\n";
	  return  0.0;
	}
      double RegionLength = RegionDepth(region);
      return RegionLength/2.0;
    }
    //______________________________________________________________________________
    HepRotation RegionRotation(int sec, int region)
    {
      using namespace clas12::geo::DC;
      HepRotation  rot;
      rot.rotateX( - RegionTilt[region-1] );
      rot.rotateZ( SectorRotation[sec-1] );
      return rot;
    }
    //______________________________________________________________________________

    Hep3Vector RegionTranslation(int sec, int region)
    {
      using namespace clas12::geo::DC;
      CLHEP::Hep3Vector  vec(
			     0.0,
			     RegionTrapCorner_y(region) - ContainerExtraShift(region) ,
			     RegionTrapCorner_z(region)
			     );
      vec.rotateZ( SectorRotation[sec-1] );
      return vec;
    }
    //______________________________________________________________________________

    double RegionTrapCorner_z(int region)
    {
      // coordinate to place the geant4 trap 
      // note this is the front corner of the edge where the endplates would meet
      using namespace clas12::geo::DC;
      int index = region-1;
      double z_sl =  SuperLayerRefWire_z(RegionSuperLayers[index][0]);
      double c = Cos( RegionTilt[index] );
      return( z_sl - FrontGap[index]*c );
    }
    //______________________________________________________________________________

    double RegionTrapCorner_y(int region)
    {
      // coordinate to place the geant4 trap
      // note this is the front corner of the edge where the endplates would meet
      using namespace clas12::geo::DC;
      int index = region-1;
      return( RegionXOffset[index] );
    }
    //______________________________________________________________________________

    double SuperLayerRefWire_z(int sl)
    {
      using namespace clas12::geo::DC;
      int index = sl-1;
      return DistanceToRefWire.at(index)*Cos(ThetaMin.at(index));
    }
    //______________________________________________________________________________

    double SuperLayerRefWire_y(int sl)
    {
      using namespace clas12::geo::DC;
      using namespace TMath; 
      int index = sl-1;
      return DistanceToRefWire.at(index)*Sin(ThetaMin.at(index));
    }
    //______________________________________________________________________________

    HepRotation LayerWireRotation(int sl)
    {
      using namespace clas12::geo::DC;
      HepRotation  rot;
      rot.rotateZ( ThetaStereo[sl-1] );
      return rot;
    }
    //______________________________________________________________________________

    Hep3Vector ToWireMidPlane(int sl, int layer, int wire)
    {
      // Returns vector that translates the midplane location 
      // relative to the refernece which is the upstream vertex of the front window
      // plane. 
      using namespace clas12::geo::DC;
      using namespace CLHEP;

      int index = SuperLayerRegionIndex.at(sl-1);
      double RegionLength = RegionDepth(index+1);//FrontGap[index] + MidGap[index] + BackGap[index];

      double arb_offset = 4*cm;
      int    i = (layer%2);
      double y = arb_offset + (0.5*double(i) + double(wire-1))*LayerWireSpacing[sl-1];
      double z = -RegionLength/2.0;
      z += SuperLayerZOffset.at(sl-1);
      z += double(layer)*3.0*LayerSep[sl-1]; Hep3Vector  vec(   0.0, y, z);
      Hep3Vector  offset(0.0, RefWireOffset.at(sl-1).y(), RefWireOffset.at(sl-1).x());
      vec += offset;
      return vec;
    }
    //______________________________________________________________________________

    double  WireLength(int sl,int layer,int wire)
    {
      // Returns the wire's length
      using namespace clas12::geo::DC;
      using namespace CLHEP;
      using namespace TMath;

      int region = (sl-1)/2 + 1;
      int index  = (sl-1)/2;
      if( index < 0 || index >= 3 )
	{
	  std::cout << "clas12::geo::WireLength(region=" << region << ") : Region out of range\n";
	  return  0.0;
	}

      Hep3Vector wire_position = ToWireMidPlane(sl,layer,wire);

      //double x       = wire_position.x() // zero
      double y1      = wire_position.y();
      double z_depth = wire_position.z();
      double RegionLength = RegionDepth(region);

      // The mid plane cross section is a rectangle plus a right triangle. 
      // extra_back is the length of the triangle on the back window plane
      //    ______
      //   /| 
      //  /_|_____ 
      //   e        
      // The line segment above e is extra_back
      // coordinate system is normal to the window faces
      //double extra_back = RegionLength*Tan(RegionTilt[index]);
      double y0 = z_depth*Tan(RegionTilt[index]);

      // This returns the nose angle (between two end plates) for the window when viewed normal to it
      double actual_plate_angle = WindowAngle( RegionAngle[index]/2.0 , RegionTilt[index] );
      //std::cout << " actual_plate_angle " << actual_plate_angle/CLHEP::degree << std::endl;

      double l_total  = (y1 + y0)*Tan(actual_plate_angle);

      // Now we rotate the wire by the stereo angle.
      double th_stereo = ThetaStereo[sl-1];
      double phi_1     = pi/2.0 - actual_plate_angle - th_stereo;
      double phi_2     = pi/2.0 + actual_plate_angle - th_stereo;

      double w1 = l_total*( Cos(th_stereo) + Sin(th_stereo)/Tan(phi_1) );
      double w2 = l_total*( Cos(th_stereo) + Sin(th_stereo)/Tan(phi_2) );

      return(w1 + w2);
    }
    //______________________________________________________________________________

    double  WireShift(int sl,int layer,int wire)
    {
      // Shift of wire's center, along the wire, due to stereo angle
      using namespace clas12::geo::DC;
      using namespace CLHEP;
      using namespace TMath;

      int region = (sl-1)/2 + 1;
      int index  = (sl-1)/2;
      if( index < 0 || index >= 3 )
	{
	  std::cout << "clas12::geo::WireLength(region=" << region << ") : Region out of range\n";
	  return  0.0;
	}

      Hep3Vector wire_position = ToWireMidPlane(sl,layer,wire);

      //double x       = wire_position.x() // zero
      double y1      = wire_position.y();
      double z_depth = wire_position.z();
      double RegionLength = RegionDepth(region);

      // The mid plane cross section is a rectangle plus a right triangle. 
      // extra_back is the length of the triangle on the back window plane
      //    ______
      //   /| 
      //  /_|_____ 
      //   e        
      // The line segment above e is extra_back
      // coordinate system is normal to the window faces
      //double extra_back = RegionLength*Tan(RegionTilt[index]);
      double y0 = z_depth*Tan(RegionTilt[index]);

      // This returns the nose angle (between two end plates) for the window when viewed normal to it
      double actual_plate_angle = WindowAngle( RegionAngle[index]/2.0 , RegionTilt[index] );
      //std::cout << " actual_plate_angle " << actual_plate_angle/CLHEP::degree << std::endl;

      double l_total  = (y1 + y0)*Tan(actual_plate_angle);

      // Now we rotate the wire by the stereo angle.
      double th_stereo = ThetaStereo[sl-1];
      double phi_1     = pi/2.0 - actual_plate_angle - th_stereo;
      double phi_2     = pi/2.0 + actual_plate_angle - th_stereo;

      double w1 = l_total*( Cos(th_stereo) + Sin(th_stereo)/Tan(phi_1) );
      double w2 = l_total*( Cos(th_stereo) + Sin(th_stereo)/Tan(phi_2) );

      return(w1 - w2);
    }
    //______________________________________________________________________________

    Hep3Vector  WireStereoShift(int sl,int layer,int wire)
    {
      // Shift of wire center along wire due to stereo angle
      using namespace clas12::geo::DC;
      using namespace CLHEP;
      using namespace TMath;

      int region = (sl-1)/2 + 1;
      int index  = (sl-1)/2;
      if( index < 0 || index >= 3 )
	{
	  std::cout << "clas12::geo::WireLength(region=" << region << ") : Region out of range\n";
	  return  {0.0,0.0,0.0};
	}

      Hep3Vector wire_position = ToWireMidPlane(sl,layer,wire);

      //double x       = wire_position.x() // zero
      double y1      = wire_position.y();
      double z_depth = wire_position.z();
      double RegionLength = RegionDepth(region);

      // The mid plane cross section is a rectangle plus a right triangle. 
      // extra_back is the length of the triangle on the back window plane
      //    ______
      //   /| 
      //  /_|_____ 
      //   e        
      // The line segment above e is extra_back
      // coordinate system is normal to the window faces
      //double extra_back = RegionLength*Tan(RegionTilt[index]);
      double y0 = z_depth*Tan(RegionTilt[index]);

      // This returns the nose angle (between two end plates) for the window when viewed normal to it
      double actual_plate_angle = WindowAngle( RegionAngle[index]/2.0 , RegionTilt[index] );
      //std::cout << " actual_plate_angle " << actual_plate_angle/CLHEP::degree << std::endl;

      double l_total  = (y1 + y0)*Tan(actual_plate_angle);

      // Now we rotate the wire by the stereo angle.
      double th_stereo = ThetaStereo[sl-1];
      double phi_1     = pi/2.0 - actual_plate_angle - th_stereo;
      double phi_2     = pi/2.0 + actual_plate_angle - th_stereo;

      double w1 = l_total*( Cos(th_stereo) + Sin(th_stereo)/Tan(phi_1) );
      double w2 = l_total*( Cos(th_stereo) + Sin(th_stereo)/Tan(phi_2) );
      double shift = (w1 - w2)/2.0;
      Hep3Vector res = {shift*Cos(th_stereo),shift*Sin(th_stereo), 0};
      return(res);
    }
    //______________________________________________________________________________

    HepRotation  GetLeftEndplateRotation(int region)
    {
      using namespace clas12::geo::DC;
      int index = region-1;
      HepRotation rot           = HepRotation::IDENTITY;
      double RegionLength       = RegionDepth(region);
      double actual_plate_angle = WindowAngle( RegionAngle[index]/2.0 , RegionTilt[index] );
      //    ______
      //   /| 
      //  /_|_____ 
      //   e        
      // The line segment above e is extra_back
      double e_extra            = RegionLength*TMath::Tan(RegionTilt[index]);
      double eta = ATan( (e_extra/RegionLength)*Tan(actual_plate_angle)*Sin(90.0*degree-actual_plate_angle));
      double nx = (NoseWidthInside[index]/2.0);
      double ny = (NoseWidthInside[index]/2.0)/TMath::Tan(actual_plate_angle);
      Hep3Vector n_rotation = {nx, ny, 0.0};
      // First orient the trap z direction 
      rot.rotateY(-90*degree);
      rot.rotateZ(-actual_plate_angle);
      rot.rotate(eta,n_rotation);
      return( rot );
    }
    //______________________________________________________________________________

    Hep3Vector  GetLeftEndplatePosition(int sec, int region)
    {
      using namespace clas12::geo::DC;
      int index = region-1;
      double actual_plate_angle = WindowAngle( RegionAngle[index]/2.0 , RegionTilt[index] );
      double x_pos = (NoseWidthOutside[index]/2.0) + EndPlateThickness.at(index)/2.0;
      double y_pos = (NoseWidthOutside[index]/2.0)/TMath::Tan(actual_plate_angle) + ContainerExtraShift(region);
      double z_pos = (ContainerTrapWidth(region)-RegionTrapWidth(region))/2.0;
      // This is the position of the upstream corner
      Hep3Vector  endplate_pos = {
	x_pos,
	y_pos,
	z_pos,
      };
      // rotate the 
      return( endplate_pos );
    }
    //______________________________________________________________________________

    HepRotation  GetRightEndplateRotation(int region)
    {
      using namespace clas12::geo::DC;
      int index = region-1;
      HepRotation rot           = HepRotation::IDENTITY;
      double RegionLength       = RegionDepth(region);
      double actual_plate_angle = WindowAngle( RegionAngle[index]/2.0 , RegionTilt[index] );
      //    ______
      //   /| 
      //  /_|_____ 
      //   e        
      // The line segment above e is extra_back
      double e_extra            = RegionLength*TMath::Tan(RegionTilt[index]);
      double eta = ATan( (e_extra/RegionLength)*Tan(actual_plate_angle)*Sin(90.0*degree-actual_plate_angle));
      double nx = (NoseWidthInside[index]/2.0);
      double ny = (NoseWidthInside[index]/2.0)/TMath::Tan(actual_plate_angle);
      Hep3Vector n_rotation = {-nx, ny, 0.0};
      // First orient the trap z direction 
      rot.rotateY(-90*degree);
      rot.rotateZ(actual_plate_angle);
      rot.rotate(-eta,n_rotation);
      return( rot );
    }
    //______________________________________________________________________________
    Hep3Vector  GetRightEndplatePosition(int sec, int region)
    {
      using namespace clas12::geo::DC;
      int index = region-1;
      double actual_plate_angle = WindowAngle( RegionAngle[index]/2.0 , RegionTilt[index] );
      double x_pos = (NoseWidthOutside[index]/2.0) + EndPlateThickness.at(index)/2.0;
      double y_pos = (NoseWidthOutside[index]/2.0)/TMath::Tan(actual_plate_angle) + ContainerExtraShift(region);
      double z_pos = (ContainerTrapWidth(region)-RegionTrapWidth(region))/2.0;
      // This is the position of the upstream corner
      Hep3Vector  endplate_pos = {
	-x_pos,
	y_pos,
	z_pos
      };
      // rotate the 
      return( endplate_pos );
    }
    //______________________________________________________________________________

    std::vector<Hep2Vector> EndplateTrapPoints(int region)
    {
      // origin will be the upstream corner
      // coordinates are y perp to wires, x goes with along increaseing layer number
      // z is normal to plate
      using namespace TMath;
      using namespace DC;
      int index = region-1;
      double bottom_extra = EndPlateWidth.at(index)*Tan(EndPlateTiltAngle.at(index));
      double top_extra = EndPlateLongSideLength.at(index) - (EndPlateShortSideLength.at(index)+bottom_extra);
      if(top_extra<0) std::cout << top_extra << " ERROR TOP EXTRA SHOULD BE POSITIVE\n";

      double w = EndPlateWidth.at(index);
      std::vector<Hep2Vector> points = {
	{0.0,0.0},
	{0.0, EndPlateShortSideLength[index]}, 
	{w  , EndPlateShortSideLength.at(index)+top_extra},
	{w  , -bottom_extra},
	{0.0,0.0},
	{0.0, EndPlateShortSideLength[index]}, 
	{w  , EndPlateShortSideLength.at(index)+top_extra},
	{w  , -bottom_extra}
      };
      return points;
    }
    //______________________________________________________________________________

    HepRotation SectorZRotation(int sec)
    {
      using namespace clas12::geo::DC;
      HepRotation  rot;
      rot.rotateZ( SectorRotation[sec-1] );
      return rot;
    }
    //______________________________________________________________________________

    double  WindowAngle(double th, double th_tilt)
    {
      // The angle of the window where the endplate meets.
      using namespace clas12::geo::DC;
      Hep3Vector Y1(0,1.0,0);
      Hep3Vector X1(1,0,0);
      X1.setRThetaPhi(1.0,90*degree,90*degree - th);
      Hep3Vector dz(0,0,TMath::Tan(th_tilt)) ;
      //std::cout << " Angle between Y and vector is " << X1.angle(Y1)/degree << std::endl;
      Hep3Vector X2  = X1 + dz;
      Hep3Vector Y2  = Y1 + dz;
      //std::cout << " Now angle is " << X2.angle(Y2)/degree << std::endl;
      return X2.angle(Y2);
    }
    //______________________________________________________________________________

    double RegionDepth(int region){
      // Distance from region's front window to its back window in the
      // tilted frame (normal to window)
      using namespace clas12::geo::DC;
      int index = region-1;
      double RegionLength = FrontGap.at(index) + MidGap.at(index) + BackGap.at(index);
      RegionLength += 21.0*LayerSep.at( RegionSuperLayerIndex.at(index).at(0) ) ;
      RegionLength += 21.0*LayerSep.at( RegionSuperLayerIndex.at(index).at(1) ) ;
      return RegionLength;
    }
    //______________________________________________________________________________

  }
}


