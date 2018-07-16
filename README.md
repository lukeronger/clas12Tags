# clas12Tags with ALERT

The ALERT geometry for gemc is included with this fork of clas12Tags. It is still very much in testing and should not be taken as a final implementation. 

The geometry was built in the ```ALERT/geometryGenAlert/alert-build``` directory, directly from Geant4 defined geometry from a modified version of Whitney Armstrong's (Argonne National Labratory) [c12sim](https://gitlab.com/whit2333/c12sim). The outputted GDML was imported to gemc in ```ALERT/alert_clas12.gcard```. The GDML has to be modified so the world volume is named "World" to properly work with gemc.
