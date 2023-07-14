#include "GmshLYSO.hh"

GmshLYSO::GmshLYSO(int Znode, double Xtot, double Ztot, double* ptsY, char* modelname)
{
    gmsh::initialize();
    gmsh::model::add(modelname);

    // Geometry
    double dZ = Ztot / Znode;
    std::vector<int> pp;
    std::vector<int> pm;
    for (int i = 0; i < Znode + 1; i++) {
        gmsh::model::geo::addPoint(-Xtot / 2, +1 * ptsY[i], -1 * Ztot + dZ * i, 0, 1000 + i);
        gmsh::model::geo::addPoint(-Xtot / 2, -1 * ptsY[i], -1 * Ztot + dZ * i, 0, 2000 + i);
        std::cout << -Xtot / 2 << " " << +1 * ptsY[i] << " " << +-1 * Ztot + dZ * i << std::endl;
std::cout << "Point 1: " << -Xtot / 2 << " " << +1 * ptsY[i] << " " << -1 * Ztot + dZ * i << std::endl;
std::cout << "Point 2: " << -Xtot / 2 << " " << -1 * ptsY[i] << " " << -1 * Ztot + dZ * i << std::endl;

        pp.push_back(1000 + i);
        pm.push_back(2000 + i);
    }
    int splp, splm,splp1, splm1, l0, lm, cl, sf, ext, extr;
    splp = gmsh::model::geo::addSpline(pp);
    splm = gmsh::model::geo::addSpline(pm);
    splp1 = gmsh::model::geo::addSpline(pp);
    splm1 = gmsh::model::geo::addSpline(pm);
    l0 = gmsh::model::geo::addLine(pp[Znode], pm[Znode]);
    lm = gmsh::model::geo::addLine(pp[0], pm[0]);
    cl = gmsh::model::geo::addCurveLoop({-lm, splp, l0, -splm});
    sf = gmsh::model::geo::addPlaneSurface({cl});
    gmsh::model::geo::symmetrize({{1,splp},{1,splm},{1,l0}}, 0,0,1,0);  // Symmetrize the extruded volume

    // Extrude the surface to create a 3D volume
    //std::vector<std::pair<int, int>> ov2;
    //gmsh::model::geo::extrude({{2, sf}}, Xtot, 0, 0, ov2, {1}, {});
/*
    // Set transfinite curves for symmetry
    gmsh::model::geo::mesh::setTransfiniteCurve(l0, 10);
    gmsh::model::geo::mesh::setTransfiniteCurve(lm, 10);

    // Meshing
    int nsecmesh = 1;
    gmsh::model::geo::mesh::setTransfiniteCurve(splp, nsecmesh + 1);
    gmsh::model::geo::mesh::setTransfiniteCurve(splm, nsecmesh + 1);
*/
    gmsh::model::geo::synchronize();
//    gmsh::model::mesh::generate(3);
}

GmshLYSO::~GmshLYSO()
{
}
