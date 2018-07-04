void viewALERT() {
  TGeoManager::Import("alert.gdml");
  gGeoManager->GetTopVolume()->Draw("ogl");
  gGeoManager->Export("output.gdml","","gf");
}
