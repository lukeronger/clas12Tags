#ifndef __DCWIRE_HH
#define __DCWIRE_HH

#include "DCGeometry.hh"

#include "TNamed.h"
#include "TBrowser.h"
#include "TMath.h"

#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Vector/Rotation.h"

#include <array>
#include <map>
#include <tuple>
#include <unordered_map>

//______________________________________________________________________________

namespace clas12 {

  namespace wire {

    class DCSuperLayer {
    public:
      int fSector     = 0;
      int fRegion     = 0;
      int fSuperLayer = 0;

    public:
      DCSuperLayer() { } 
      DCSuperLayer(int s, int r, int sl) :
	fSector(s), fRegion(r), fSuperLayer(sl) { }
      DCSuperLayer(const DCSuperLayer&) = default;
      DCSuperLayer(DCSuperLayer&&)      = default;
      DCSuperLayer& operator=(const DCSuperLayer&) = default;
      DCSuperLayer& operator=(DCSuperLayer&&)      = default;
      virtual ~DCSuperLayer() = default;

      int GlobalChannel() const { return( 6*(fSector-1) + (fSuperLayer-1) ); }

      bool operator< (const DCSuperLayer& rhs) const {return( GlobalChannel() < rhs.GlobalChannel() ); }
      bool operator> (const DCSuperLayer& rhs) const {return rhs < (*this);}
      bool operator<=(const DCSuperLayer& rhs) const {return !((*this) > rhs);}
      bool operator>=(const DCSuperLayer& rhs) const {return !((*this) < rhs);}
      bool operator==(const DCSuperLayer& rhs) const {return( GlobalChannel() == rhs.GlobalChannel() ); }
      bool operator!=(const DCSuperLayer& rhs) const {return !((*this) == rhs);}

      ClassDef(DCSuperLayer,1)
    };
    //________________________________________________________________________

    // Wire channel numbering convention
    // the smallest channel number in a layer starts 
    // nearest the beam pipe.
    // Numbering starts first by completing a sector.
    // Each sector has 6 superlayers, and each superlayer
    // has 6 layers, each layer has 
    class DCWire {

    public :
      int  fSector;
      int  fRegion;
      int  fSuperLayer;   // redundant information Region I: (1,2), Region II: (3,4) R III: (5,6)
      int  fLayer;
      int  fWire;       
      int  fChannel; // Unique channel number

      DCWire(int sec=0, int reg=0, int sl=0, int l=0, int w=0);
      virtual ~DCWire();

      static int GetChannel(int sector, int sl, int layer, int wire) {
	int sec_ind  = sector-1;
	int sl_ind   = sl-1;
	int lay_ind  = layer-1;
	int wire_ind = wire-1;

	const int WiresPerLayer   = 112;
	const int WiresPerSL      = 6*112;
	const int WiresPerSector  = 6*6*112;
	const int TotalWires      = 6*6*6*112;

	int res = WiresPerSector*sec_ind + WiresPerSL*sl_ind + WiresPerLayer*lay_ind + WiresPerLayer*wire_ind;
	return res;
      }

      static int GetSector    (int channel) {
	//const int WiresPerLayer   = 112;
	//const int WiresPerSL      = 6*112;
	const int WiresPerSector  = 6*6*112;
	int sec_ind  = channel/WiresPerSector;
	return( sec_ind+1 );
      }
      static int GetSuperLayer(int channel) {
	//const int WiresPerLayer   = 112;
	const int WiresPerSL      = 6*112;
	const int WiresPerSector  = 6*6*112;
	int sec_ind   = channel/WiresPerSector;
	int sec_chan  = channel - sec_ind*WiresPerSector;
	int sl_ind    = sec_chan/WiresPerSL;
	return( sl_ind+1 );
      }
      static int GetLayer     (int channel) {
	const int WiresPerLayer   = 112;
	const int WiresPerSL      = 6*112;
	const int WiresPerSector  = 6*6*112;
	int sec_ind   = channel/WiresPerSector;
	int sec_chan  = channel - sec_ind*WiresPerSector;
	int sl_ind    = sec_chan/WiresPerSL;
	int sl_chan   = sec_chan - sl_ind*WiresPerSL;
	int lay_ind   = sl_chan/WiresPerLayer;
	return( lay_ind+1 );
      }
      static int GetWire      (int channel) {
	const int WiresPerLayer   = 112;
	const int WiresPerSL      = 6*112;
	const int WiresPerSector  = 6*6*112;
	int sec_ind   = channel/WiresPerSector;
	int sec_chan  = channel - sec_ind*WiresPerSector;
	int sl_ind    = sec_chan/WiresPerSL;
	int sl_chan   = sec_chan - sl_ind*WiresPerSL;
	int lay_ind   = sl_chan/WiresPerLayer;
	int wire_ind  = sl_chan - lay_ind*WiresPerLayer ;
	return( wire_ind+1 );
      }
      static int GetRegion    (int channel) {
	return( GetSuperLayer(channel)/2 + 1);
      }

      void Print(Option_t * opt = "") const;

      ClassDef(DCWire,3)
    };

  }


}

#endif

