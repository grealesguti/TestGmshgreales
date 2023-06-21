#include "GmshLYSO.hh"

GmshLYSO :: GmshLYSO(int Znode, double Xtot, double Ztot, double* ptsY, char* modelname)
{
	gmsh::initialize();
	gmsh::model::add(modelname);

// Geometry
  double dZ=Ztot/Znode;
  std::vector<int> pp;
  std::vector<int> pm;
  for(int i = 0; i < Znode+1; i++) {
    gmsh::model::geo::addPoint(-Xtot/2, +1*ptsY[i],-1*Ztot+dZ*i ,0,
                               1000 + i);
    gmsh::model::geo::addPoint(-Xtot/2, -1*ptsY[i],-1*Ztot+dZ*i ,0,
                               2000 + i);
    std::cout<<-Xtot/2<< " " <<+1*ptsY[i]<< " " <<+-1*Ztot+dZ*i<< std::endl;
    pp.push_back(1000 + i);
    pm.push_back(2000 + i);
  }
  int splp, splm, l0, lm, cl, sf, ext;
splp=gmsh::model::geo::addSpline(pp);
splm=gmsh::model::geo::addSpline(pm);
l0 = gmsh::model::geo::addLine(pp[Znode], pm[Znode]);
lm = gmsh::model::geo::addLine(pp[0], pm[0]);
cl = gmsh::model::geo::addCurveLoop({-lm, splp, l0, -splm});
sf = gmsh::model::geo::addPlaneSurface({cl});
std::vector<std::pair<int, int> > ov2;
gmsh::model::geo::extrude({{2, sf}}, Xtot, 0, 0,ov2,{1},{});
gmsh::model::geo::mesh::setTransfiniteCurve(l0,2);
gmsh::model::geo::mesh::setTransfiniteCurve(lm,2);

// Meshing
int nsecmesh=1;
gmsh::model::geo::mesh::setTransfiniteCurve(splp,nsecmesh+1);
gmsh::model::geo::mesh::setTransfiniteCurve(splm,nsecmesh+1);
gmsh::model::geo::synchronize();
gmsh::model::mesh::generate(3);

// Getting nodal information

  // Print the model name and dimension:
  std::string name;
  gmsh::model::getCurrent(name);
  std::cout << "Model " << name << " (" << gmsh::model::getDimension()
            << "D)\n";

  std::vector<std::pair<int, int> > entities;
  gmsh::model::getEntities(entities);

  for(auto e : entities) {
    // Dimension and tag of the entity:
    int dim = e.first, tag = e.second;

    // Mesh data is made of `elements' (points, lines, triangles, ...), defined
    // by an ordered list of their `nodes'. Elements and nodes are identified by
    // `tags' as well (strictly positive identification numbers), and are stored
    // ("classified") in the model entity they discretize. Tags for elements and
    // nodes are globally unique (and not only per dimension, like entities).

    // A model entity of dimension 0 (a geometrical point) will contain a mesh
    // element of type point, as well as a mesh node. A model curve will contain
    // line elements as well as its interior nodes, while its boundary nodes
    // will be stored in the bounding model points. A model surface will contain
    // triangular and/or quadrangular elements and all the nodes not classified
    // on its boundary or on its embedded entities. A model volume will contain
    // tetrahedra, hexahedra, etc. and all the nodes not classified on its
    // boundary or on its embedded entities.

    // Get the mesh nodes for the entity (dim, tag):
    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoords, nodeParams;
    gmsh::model::mesh::getNodes(nodeTags, nodeCoords, nodeParams, dim, tag);

    // Get the mesh elements for the entity (dim, tag):
    std::vector<int> elemTypes;
    gmsh::model::mesh::getElements(elemTypes, elemTags, elemNodeTags, dim, tag);

    std::string type;
    gmsh::model::getType(dim, tag, type);
    std::string name;
    gmsh::model::getEntityName(dim, tag, name);
    if(name.size()) name += " ";
    if(type == "Volume"){
		std::cout << "Entity " << name << "(" << dim << "," << tag << ") of type "
				  << type << "\n";

		// * Number of mesh nodes and elements:
		int numElem = 0;
		for(auto &tags : elemTags) numElem += tags.size();
		std::cout << " - Mesh has " << nodeTags.size() << " nodes and " << numElem
				  << " elements\n";

		// * List all types of elements making up the mesh of the entity:
		for(auto elemType : elemTypes) {
		  std::string name;
		  int d, order, numv, numpv;
		  std::vector<double> param;
		  gmsh::model::mesh::getElementProperties(elemType, name, d, order, numv,
												  param, numpv);
		  std::cout << " - Element type: " << name << ", order " << order << "\n";
		  std::cout << "   with " << numv << " nodes in param coord: (";
		  for(auto p : param) std::cout << p << " ";
		  std::cout << ")\n";
		}
		//gmsh::model::mesh::getElements(elemTypes, elemTags, elemNodeTags, dim, tag);

		for(auto eNodeTags : elemNodeTags) {
			std::cout<<" - Nodal Tags:: ";
			for (int i : eNodeTags) std::cout<< i<< " ";
			std::cout<<std::endl;
		}
		//gmsh::model::mesh::getNodes(nodeTags, nodeCoords, nodeParams, dim, tag);
		gmsh::model::mesh::getNodes(tags, coord, param);
		std::cout<<tags.size()<<std::endl;
				std::cout<<coord.size()<<std::endl;
    for(std::size_t i = 0; i < coord.size(); i += 3) {
	std::cout<<"Coords::"<<coord[i]<<" "<<coord[i+1]<<" "<<coord[i+2]<<std::endl;
    }
    for(std::size_t i = 0; i < tags.size(); i += 1) {
	std::cout<<"Tags:: "<<i<<" ";
    }
	std::cout<<std::endl;
		
	}
  }
	
	}
	
GmshLYSO :: ~GmshLYSO()
{
	}

// Given a material assignment creates the logical and solid volumes 
// from the gmsh mesh
void GmshLYSO ::CreateG4LYSO(G4Material *material, G4LogicalVolume *logicWorld){
		int etag[4], gidx;
	double x[4],y[4],z[4];
	gidx = 1000;
	G4String tetname="G4Tet_";
	for(itd::size_t i = 0; i < elemNodeTags.size(); i += 4) {
		etag[0] = elemNodeTags[i];etag[1] = elemNodeTags[i+1];etag[2] = elemNodeTags[i+2];etag[3] = elemNodeTags[i+3];
        for(int j=0;j<4;j++){
            x[j]=coords[etag[j]*3+j];
            y[j]=coords[etag[j]*3+j+1];
            z[j]=coords[etag[j]*3+j+2];
			}
			p1 = new G4ThreeVector(x[0],y[0],z[0]);
			p2 = new G4ThreeVector(x[1],y[1],z[1]);
			p3 = new G4ThreeVector(x[2],y[2],z[2]);
			p4 = new G4ThreeVector(x[3],y[3],z[3]);
			LYSOTet_Solid = new G4Tet(tetname+ G4String("solid_")+G4UIcommand::ConvertToString(gidx),p1,p2,p3,p4,0);
			LYSOTet_Logic = new G4LogicalVolume(LYSOTet_Solid, material, tetname+ G4String("logical_")+G4UIcommand::ConvertToString(gidx))
			LYSOTet_Phys = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),LYSOTet_Logic,tetname+ G4String("phys_")+G4UIcommand::ConvertToString(gidx),logicWorld,false,0,true);       
			gidx+=1;
		}
	
	
	}
