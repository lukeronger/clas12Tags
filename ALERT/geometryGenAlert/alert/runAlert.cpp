#include "ALERTConstruction.hh"

#include "G4RunManager.hh"
#include "QBBC.hh"
#include "G4TransportationManager.hh"
#include "G4GDMLParser.hh"

int main()
{
  G4RunManager * runManager = new G4RunManager;
  
  runManager->SetUserInitialization(new ALERT());
  
  G4VModularPhysicsList* physicsList = new QBBC;
  
  runManager->SetUserInitialization(physicsList);

  runManager->Initialize();

  G4VPhysicalVolume* pWorld = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->GetWorldVolume();

  G4GDMLParser Parser;
  
  Parser.Write("alert.gdml",pWorld); 
 
  delete runManager;  
  return 0;
}
