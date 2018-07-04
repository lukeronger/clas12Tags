#include "RCGeometry.hh"

#include "CLHEP/Geometry/Transform3D.h"

#include "TMath.h"

#include <iomanip>


namespace clas12 {

  namespace  wire {

    using namespace CLHEP;
    //using namespace HepGeom;

    RCGeometry::RCGeometry()
    {
      // Moving radial the structure goes like this:
      // window | gap | field wire | field gap | sense wire | field gap | field wire | gap | field wire ...
      // g fsf g ... g fsf g 
      fNChannels = 0;
      for(int i = 0; i< NLayers; i++){
	// The inner radius for a given layer of wires
	fCellInnerRadius[i]     = (InnerRadius + InnerGap + double(i)*(2.0*DeltaR + LayerGap ));
	fNCells[i]              = int(CLHEP::pi*fCellInnerRadius[i]/MaxWireSep);
	fCellDeltaPhi[i]        = (2.*CLHEP::pi/double(fNCells[i]));
	fFirstChannelInLayer[i] = fNChannels;
	fNChannels += fNCells[i];
	double sign = 1.0;
	if(i%2==0) sign=1.;
	else sign=-1.;
	fCellStereoAngle[i] = sign*StereoAngle;
	fCellCentralRadius[i] = fCellInnerRadius[i] + DeltaR;
      }

    }
    //______________________________________________________________________________

    void RCGeometry::Print()
    {
      std::cout << "layer   Radius   cells  start_channel delta_phi   stereo"  << std::endl;
      for(int i = 0; i< NLayers; i++){
	std::cout << i << " "
		  << fCellInnerRadius[i]/CLHEP::mm << " " 
		  << fNCells[i]       << " " 
		  << fFirstChannelInLayer[i]       << " " 
		  << fCellDeltaPhi[i]/CLHEP::degree << " " 
		  << fCellStereoAngle[i]/CLHEP::degree << std::endl;
      }
    }
    //______________________________________________________________________________

    int  RCGeometry::GetWireLayer(int channel) const {
      int res  = 0;
      int chan = channel;
      for(auto n: fNCells) {
	//std::cout << n  << " cells" << std::endl;
	if( chan < n) break;
	chan -= n;
	res++;
      }  
      return res;
    }
    //______________________________________________________________________________

    double RCGeometry::GetSenseWireRadius(int channel) const {
      //std::cout << " ch " << channel << std::endl;
      int lay = GetWireLayer(channel);
      //std::cout << " lay " << lay << std::endl;
      return fCellCentralRadius.at(lay);
    }
    //______________________________________________________________________________

    double RCGeometry::GetSenseWirePhi(int channel) const {
      //std::cout << " ch " << channel << std::endl;
      int     lay = GetWireLayer(channel);
      //std::cout << " lay " << lay << std::endl;
      int     wire_number = channel - fFirstChannelInLayer.at(lay);
      return wire_number*fCellDeltaPhi.at(lay);
    }
    //______________________________________________________________________________

    double RCGeometry::GetSenseWireStereoAngle(int channel) const {
      int     lay = GetWireLayer(channel);
      return fCellStereoAngle.at(lay);
    }
    //______________________________________________________________________________

    Hep3Vector RCGeometry::GetFirstCellPosition(int layer, int subcell) const {
      double r    = fCellCentralRadius.at(layer);
      double dphi = 0.0;
      //std::cout << " subcell " << subcell << std::endl;
      if( subcell != -1 ){
	if( (subcell == 1) || (subcell == 2) ) {
	  dphi = fCellDeltaPhi.at(layer)/4.0;
	} else {
	  dphi = -fCellDeltaPhi.at(layer)/4.0;
	}
	//if( subcell > 1) {
	//   r += DeltaR;
	//}
      }
      //std::cout << " dph " << dphi/degree << std::endl;
      Hep3Vector res = {r,0,0};
      res.rotateZ(-dphi);
      return res;
    }
    //______________________________________________________________________________

    HepRotation RCGeometry::GetFirstCellRotation(int layer, int subcell) const {
      HepRotation rot = HepRotation::IDENTITY;
      double dphi = 0.0;
      if( subcell != -1 ){
	if( (subcell == 1) || (subcell == 2) ) {
	  dphi = fCellDeltaPhi.at(layer)/4.0;
	} else {
	  dphi = -fCellDeltaPhi.at(layer)/4.0;
	}
      }
      rot.rotateZ(dphi);
      rot.rotateX(-1.0*fCellStereoAngle.at(layer));
      return rot;
    }
    //______________________________________________________________________________

    HepGeom::Transform3D RCGeometry::GetFirstCellTransform(int layer, int subcell) const {
      HepRotation rot = GetFirstCellRotation(layer,subcell);
      return( HepGeom::Transform3D(rot, Hep3Vector(0,0,0))*
	      HepGeom::Transform3D( HepRotation::IDENTITY,
				    -1.0*GetFirstCellPosition(layer,subcell) ) ) ;
    }
    //______________________________________________________________________________

    Hep3Vector RCGeometry::GetSenseWirePerp(int channel) const {
      //returns a vector perpendicular to \vec{r}_wire and \vec{v}_wire , where the
      // later is the vector along the 
      Hep3Vector res = {0,1,0};
      res.rotateX( GetSenseWireStereoAngle(channel) );
      res.rotateZ( GetSenseWirePhi(channel) );
      return res;
    }
    //______________________________________________________________________________

    Hep3Vector RCGeometry::GetSenseWirePosition(int channel) const {
      double  r   = GetSenseWireRadius(channel);
      Hep3Vector res = {r,0,0};
      res.rotateZ(GetSenseWirePhi(channel));
      return res;
    }
    //______________________________________________________________________________

    HepRotation RCGeometry::GetSenseWireRotation(int channel) const {
      int     lay = GetWireLayer(channel);
      double  phi = GetSenseWirePhi(channel);
      Hep3Vector res = {1,0,0};
      res.rotateZ(phi);
      HepRotation rot = HepRotation::IDENTITY;
      rot.rotateZ(phi);
      rot.rotate(fCellStereoAngle.at(lay),res);
      return rot;
    }
    //______________________________________________________________________________

    HepGeom::Transform3D RCGeometry::GetSenseWireTransform(int channel) const {
      return HepGeom::Transform3D(
				  GetSenseWireRotation(channel),
				  GetSenseWirePosition(channel) );
    }
    //______________________________________________________________________________

    Hep3Vector RCGeometry::GetUpstreamSenseWirePosition(int channel) const {
      int     lay = GetWireLayer(channel);
      double  r   = GetSenseWireRadius(channel);
      //Hep3Vector res2 = {0,0,-1.0*WireLength/2.0};
      //Hep3Vector res = {r,0,0};
      //res.rotateZ( GetSenseWirePhi(channel) );
      //res2.rotateZ( GetSenseWirePhi(channel) );
      //res2.rotateX( fCellStereoAngle.at(lay) );
      //res += res2;
      Hep3Vector res2 = { 0, 0, -1.0*WireLength/2.0 };
      res2.rotateX( fCellStereoAngle.at(lay)   );
      Hep3Vector res = {r,0,0};
      res += res2;
      res.rotateZ( GetSenseWirePhi(channel) );
      return res;
    }
    //______________________________________________________________________________

    Hep3Vector RCGeometry::GetDownstreamSenseWirePosition(int channel) const {
      int     lay = GetWireLayer(channel);
      double  r   = GetSenseWireRadius(channel);
      Hep3Vector res2 = { 0, 0, WireLength/2.0 };
      res2.rotateX( fCellStereoAngle.at(lay)   );
      Hep3Vector res = {r,0,0};
      res += res2;
      res.rotateZ( GetSenseWirePhi(channel) );
      return res;
    }
    //______________________________________________________________________________

    Hep3Vector RCGeometry::GetUnitSubCellPoint(int layer,int subcell, int point) const {
      double  r   = fCellInnerRadius.at(layer);
      double phi_sign = -1.0;
      double Lsign = -1.0;
      int    p    = point%4;
      int    up   = point/4;
      if( (p == 1) || (p == 2) ) {
	phi_sign  = 1.0;
      }
      if( p > 1 ) {
	r += DeltaR;
      }
      if( subcell > 1) {
	r += DeltaR;
      }
      if( up > 0 ) {
	Lsign = 1.0;
      }
      Hep3Vector res2 = { 0, 0, Lsign*WireLength/2.0 };
      res2.rotateX( fCellStereoAngle.at(layer)   );
      Hep3Vector res = {r,0,0};
      res += res2;
      res.rotateZ( phi_sign*fCellDeltaPhi.at(layer)/4.0 );
      return( res );
    }
    //______________________________________________________________________________

    std::vector<Hep2Vector> RCGeometry::GetSubCellTrapPoints(int layer, int subcell) {
      std::vector<Hep2Vector>  res;
      for(int i = 0; i <8; i++){
	Hep3Vector p = GetUnitSubCellPoint(layer, subcell, i);
	res.push_back( {p.x(),p.y()} );
      }
      return res;
    }
    //______________________________________________________________________________

  }
}
