// -----------------------------------------------------------------------------
//
//  Gmsh C++ tutorial 1
//
//  Geometry basics, elementary entities, physical groups
//
// -----------------------------------------------------------------------------

#include <set>

// The Gmsh C++ API is entirely defined in the `gmsh.h' header (which contains
// the full documentation of all the functions in the API):
#include <gmsh.h>
#include <iostream>
#include "src_G4/GmshLYSO.hh"


int main(int argc, char **argv)
{ 

  int Znode=2;
  double ptsY[Znode+1];
  for(int i =0;i<3;i++){
	  ptsY[i]=3.;
  }
  ptsY[1]=2.;

  double Ztot=28.5;
  double Xtot=3;

  GmshLYSO *LYSOmesh = new GmshLYSO( Znode, Xtot, Ztot, ptsY, "testname");

// GUI Plotting
  std::set<std::string> args(argv, argv + argc);
  if(!args.count("-nopopup")) gmsh::fltk::run();
  gmsh::finalize();



  return 0;
}
