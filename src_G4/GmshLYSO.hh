#ifndef GMSHLYSO_HH
#define GMSHLYSO_HH

#include <set>
#include <gmsh.h>
#include <iostream>


#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UserLimits.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4Tet.hh"
#include <G4UIcommand.hh>


class GmshLYSO 
{

    void CreateG4LYSO(G4Material *material, G4LogicalVolume *logicWorld);

    double* GetNodeTags() const {return mesh_nodetags;}
    double* GetNodeCoords() const {return mesh_nodecoords;}
    int GetNumberOfNodes() const {return mesh_totalnodes;}


public:
    GmshLYSO(int Znode, double Xtot, double Ztot,double* ptsY, char* modelname);
    ~GmshLYSO();
    
private:

int mesh_totalnodes;
double *mesh_nodetags, *mesh_nodecoords;
G4Tet *LYSOTet;

std::vector<std::size_t> tags;
std::vector<double> coord, param;
std::vector<std::vector<std::size_t> > elemTags, elemNodeTags;

G4ThreeVector *p1,*p2,*p3,*p4;


};    

#endif
