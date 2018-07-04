#include "G4GDMLParser.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VisAttributes.hh"

#include <iostream>

using namespace std;

int main (){

  G4GDMLParser parser;

  try {
    parser.Read("test.gdml");
    G4VPhysicalVolume* W = parser.GetWorldVolume(); //world volume
    W->GetLogicalVolume()->SetVisAttributes(G4VisAttributes::Invisible);
  }catch (const char* eMsg){
    cout << eMsg << endl;
  }
  return 0;
}
