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
#include <iostream>
#include <vector>
#include <algorithm>
#include <G4UIcommand.hh>

std::vector<std::vector<double>> appendReversedInitialVector(const std::vector<std::vector<double>>& inputVector)
{
    std::vector<std::vector<double>> result = inputVector;  // Make a copy of the input vector

    if (!result.empty())
    {
        std::vector<double>& initialVector = result.front();  // Get a reference to the initial vector

        if (initialVector.size() > 1)
        {
            std::reverse(initialVector.begin(), initialVector.end() - 1);  // Reverse the initial vector excluding its last value
        }
    }

    return result;  // Return the modified vector
}

std::vector<std::vector<double>> reverseSubvectors(const std::vector<std::vector<double>>& inputVector) {
    std::vector<std::vector<double>> resultVector = inputVector; // Create a copy of the input vector
    
    std::reverse(resultVector.begin(), resultVector.end()); // Reverse the order of the subvectors
    
    return resultVector; // Return the resulting vector of vectors
}

std::vector<std::vector<double>> combineAllIndices(const std::vector<double>& vec1, const std::vector<double>& vec2) {
    std::vector<std::vector<double>> combinedVector;
    
    size_t size = vec1.size();
    
    for (size_t n = 0; n < size; ++n) {
        std::vector<double> innerVec;
        innerVec.push_back(vec2[n]);
        innerVec.push_back(vec2[(n + 1) % size]);
        innerVec.push_back(vec1[(n + 1) % size]);
        innerVec.push_back(vec1[n]);
        
        combinedVector.push_back(innerVec);
    }
    
    return combinedVector;
}

std::vector<std::vector<int>> combineAllIndicesint(const std::vector<int>& vec1, const std::vector<int>& vec2) {
    std::vector<std::vector<int>> combinedVector;
    
    size_t size = vec1.size();
    
    for (size_t n = 0; n < size; ++n) {
        std::vector<int> innerVec;
        innerVec.push_back(vec1[n]);
        innerVec.push_back(vec1[(n + 1) % size]);
        innerVec.push_back(vec2[(n + 1) % size]);
        innerVec.push_back(vec2[n]);
        
        combinedVector.push_back(innerVec);
    }
    
    return combinedVector;
}
std::vector<std::vector<double>> calculatePoints(const std::vector<double>& Y, double Xin, double Zmin, double Zmax, double Yzero)
{
    std::vector<std::vector<double>> points;

    // Calculate the step size for Z values
    double step = (Zmax - Zmin) / (Y.size() - 1);

    // Create symmetrical points with respect to Zmin
    for (int i = Y.size() - 1; i >= 1; --i)
    {
        double Z = Zmin - (i ) * step;  // Symmetrical Z values
        points.push_back({Xin, Y[i], Z});
        std::cout<< Xin<< " "<< Y[i]<< " "<< Z<< std::endl;
    }

    for (size_t i = 0; i < Y.size(); ++i)
    {
        double Z = Zmin + i * step;
        points.push_back({Xin, Y[i], Z});
        std::cout<< Xin<< " "<< Y[i]<< " "<< Z<< std::endl;
    }


    std::vector<std::vector<double>> additionalPoints = {{Xin, Yzero, Zmax}, {Xin, Yzero, -Zmax}};
    points.insert(points.end(), additionalPoints.begin(), additionalPoints.end());

    return points;
}

std::vector<int> createGmshPoints(const std::vector<std::vector<double>>& points, int taginit)
{
    std::vector<int> pointTags;

    // Create points
    int count = 1 ;
    for (const auto& point : points)
    {
        double x = point[0];
        double y = point[1];
        double z = point[2];

        int tag;
        tag = taginit+count;
        gmsh::model::geo::addPoint(x, y, z, 0.0, tag);
        pointTags.push_back(tag);
        count=count+1;
    }

    return pointTags;
}

std::vector<int> createGmshLines(const std::vector<int>& pointTags, int taginit, const std::set<int>& excludedIndices = {})
{
    std::vector<int> lineTags;
    int lineTag;

    int count = 1,l0;
    int pointCount = pointTags.size();

    // Create lines between successive points, excluding the lines with excluded indices
    for (int i = 0; i < pointCount - 1; i++)
    {
        if (excludedIndices.count(i) > 0)
            continue;

        lineTag = taginit + count;
        std::cout<<lineTag<<" ";
        l0=gmsh::model::geo::addLine(pointTags[i], pointTags[i + 1], lineTag);
        gmsh::model::geo::mesh::setTransfiniteCurve(l0,1);

        lineTags.push_back(lineTag);
        count++;
    }

    // Create line between the last and first points, excluding the line with the excluded index
    if (pointCount > 1 && excludedIndices.count(pointCount - 1) == 0)
    {
        lineTag = taginit + count;
        std::cout<<lineTag<<" ";
        gmsh::model::geo::addLine(pointTags[pointCount - 1], pointTags[0], lineTag);
        lineTags.push_back(lineTag);
    }
	std::cout<<std::endl;
    return lineTags;
}


int createGmshSurface(const std::vector<int>& lineTags, int taginit)
{
    int surfaceTag = taginit;

    // Create line loop
    int lineLoopTag = surfaceTag;
    gmsh::model::geo::addCurveLoop(lineTags, lineLoopTag);

    // Create surface
    gmsh::model::geo::addPlaneSurface({lineLoopTag}, surfaceTag);

    return surfaceTag;
}
std::vector<double> extractRange(const std::vector<double>& Y_all, int sec, int nodesec)
{
    std::vector<double> Y(Y_all.begin() + sec*nodesec, Y_all.begin() + sec*nodesec + nodesec);
    return Y;
}
void printVectorOfVectors(const std::vector<std::vector<int>>& vec)
{
	int i=0;
    for (const auto& innerVec : vec)
    {
		std::cout<< i<< " ";
        for (const auto& value : innerVec)
        {
            std::cout << value << " ";
        }
        std::cout << std::endl;
        i++;
    }
}

void printVector(const std::vector<int>& vec)
{
    for (const auto& value : vec)
    {
        std::cout << value << " ";
    }
    std::cout << std::endl;
}
void printVectordouble(const std::vector<double>& vec)
{
    for (const auto& value : vec)
    {
        std::cout << value << " ";
    }
    std::cout << std::endl;
}
std::vector<int> getElementsAtIndex(const std::vector<std::vector<int>>& vec, int index, int a, int b)
{
    std::vector<int> result;

    for (int i = a; i <= b; ++i)
    {
        if (i >= 0 && i < vec.size())
        {
            const auto& subvector = vec[i];
            if (index >= 0 && index < subvector.size())
            {
                result.push_back(subvector[index]);
            }
        }
    }

    return result;
}

std::vector<int> getOddElementsAtIndex(const std::vector<std::vector<int>>& vec, int index, int a, int b)
{    
	std::vector<int> result;

    for (int i = a; i <= b; ++i)
    {
        if (i >= 0 && i < vec.size())
        {
            const auto& subvector = vec[i];
            if (index >= 0 && index < subvector.size())
            {
                if (i % 2 != 0) // Check if the value is odd
                {
                    result.push_back(subvector[index]);
                }
            }
        }
    }

    return result;
}

std::vector<int> getEvenElementsAtIndex(const std::vector<std::vector<int>>& vec, int index, int a, int b)
{    
	std::vector<int> result;

    for (int i = a; i <= b; ++i)
    {
        if (i >= 0 && i < vec.size())
        {
            const auto& subvector = vec[i];
            if (index >= 0 && index < subvector.size())
            {
                if (i % 2 == 0) // Check if the value is odd
                {
                    result.push_back(subvector[index]);
                }
            }
        }
    }

    return result;
}

std::vector<double> generateEquispacedSegments(double Xmin, double Xmax, int numElements)
{
    std::vector<double> segments;

    if (numElements <= 1)
    {
        segments.push_back(Xmin);
        return segments;
    }

    double step = (Xmax - Xmin) / (numElements - 1);
    for (int i = 0; i < numElements; ++i)
    {
        double X = Xmin + i * step;
        segments.push_back(X);
    }

    return segments;
}
std::vector<int> getValuesAtIndices(const std::vector<int>& vec, const std::vector<int>& indices) {
    std::vector<int> result;
    for (int index : indices) {
        if (index >= 0 && index < vec.size()) {
            result.push_back(vec[index]);
        }
    }
    return result;
}
void generateTile(double Xmin, double Xmax, int nX, double Zmin, double Zmax, double Yzero, std::vector<double> Y_all, int nodesec  ){

    std::vector<double> Xins = generateEquispacedSegments(Xmin, Xmax, nX);
    std::vector<double> XSym = generateEquispacedSegments(Xmax+(Xins[1]-Xins[0]), Xmax+(Xmax-Xmin), nX-1);
	std::reverse(XSym.begin(), XSym.end());
	
	printVectordouble(Xins);
	printVectordouble(XSym);
	
    int surfaceTagi,idx;
	std::vector<std::vector<int>> lineTagSecs,lineTag4Surf,PointTagSecs, combinedpts, lineTagsAll;
	std::vector<int> SurfTagSecs, lineTags1,lineTags2, pointTags, line4, ptsv;
	int ptc=1000,ltc=2000,stc=3000;
	
	
	for (int sec = 0; sec < Xins.size(); sec++) {
		std::vector<double> Y = extractRange(Y_all, sec, nodesec);
		double Xin = Xins[sec];

		std::vector<std::vector<double>> points = calculatePoints(Y, Xin, Zmin, Zmax, Yzero);
		pointTags = createGmshPoints(points, ptc);
		PointTagSecs.push_back(pointTags);
		ptc=pointTags.back();
	}
	
	for (int sec = Xins.size()-2; sec > -1 ; sec--) {
		std::vector<double> Y = extractRange(Y_all, sec, nodesec);
		double Xin = XSym[sec];

		std::vector<std::vector<double>> points = calculatePoints(Y, Xin, Zmin, Zmax, Yzero);
		pointTags = createGmshPoints(points, ptc);
		PointTagSecs.push_back(pointTags);
		ptc=pointTags.back();
	}    
	
	std::vector<int> tr1,tr2;
	//////////////////
    // Creating Lines
    //////////////////
	std::cout<<"nodesec iter: "<<nodesec<<std::endl;
    int nodesec1=nodesec*2-1;
	std::cout<<"nodesec iter: "<<nodesec1<<std::endl;

	for (int sec = 0; sec < Xins.size()*2-2; sec++) {
		combinedpts =combineAllIndicesint(PointTagSecs[sec], PointTagSecs[sec+1]);
		for (int surf = 0; surf < nodesec1+2; surf++) {
			if (sec==0 && surf==0){
				ptsv=combinedpts[surf];
				tr1 = getValuesAtIndices(ptsv,{0,1,2});
				tr2 = getValuesAtIndices(ptsv,{0,2,3});

				lineTags1 = createGmshLines(tr1, ltc,{});
				lineTagsAll.push_back(lineTags1);
				line4=lineTags1;
				lineTag4Surf.push_back(line4);
				ltc=lineTags1.back();

				lineTags2 = createGmshLines(tr2, ltc,{0});
				lineTagsAll.push_back(lineTags2);
				line4={-lineTags1[2],lineTags2[0],lineTags2[1]};
				lineTag4Surf.push_back(line4);
				ltc=lineTags2.back();
			}else if (sec==0 && surf<nodesec1+1){
				ptsv=combinedpts[surf];

				if(surf<nodesec-1){
					tr1 = getValuesAtIndices(ptsv,{0,1,2});
					tr2 = getValuesAtIndices(ptsv,{0,2,3});
					
					lineTags1 = createGmshLines(tr1, ltc,{});
					lineTagsAll.push_back(lineTags1);
					line4=lineTags1;
					lineTag4Surf.push_back(line4);
					ltc=lineTags1.back();

					lineTags2 = createGmshLines(tr2, ltc,{0,2});
					lineTagsAll.push_back(lineTags2);
					line4={-lineTags1[2],lineTags2[0],-lineTag4Surf[(surf-1)*2+0][1]};
					lineTag4Surf.push_back(line4);
					ltc=lineTags2.back();
				}else{
					tr1 = getValuesAtIndices(ptsv,{0,1,3});
					tr2 = getValuesAtIndices(ptsv,{1,2,3});

					lineTags1 = createGmshLines(tr1, ltc,{2});
					lineTagsAll.push_back(lineTags1);
					ltc=lineTags1.back();
					lineTags2 = createGmshLines(tr2, ltc,{2});
					lineTagsAll.push_back(lineTags2);
					ltc=lineTags2.back();
					
					line4={lineTags1[0],lineTags1[1],-lineTag4Surf[(surf-1)*2+0][1]};
					lineTag4Surf.push_back(line4);
					
					line4={-lineTags1[1],lineTags2[0],lineTags2[1]};
					lineTag4Surf.push_back(line4);
					}
			}else if(sec==0 && surf==nodesec1+1){
				ptsv=combinedpts[surf];
				tr1 = getValuesAtIndices(ptsv,{0,1,2});
				tr2 = getValuesAtIndices(ptsv,{0,2,3});

				lineTags1 = createGmshLines(tr1, ltc,{1});
				lineTagsAll.push_back(lineTags1);
				line4={lineTags1[0],-lineTag4Surf[1][2],lineTags1[1]};
				lineTag4Surf.push_back(line4);
				ltc=lineTags1.back();

				lineTags2 = createGmshLines(tr2, ltc,{0,2});
				lineTagsAll.push_back(lineTags2);
				line4={-lineTags1[1],lineTags2[0],-lineTag4Surf[(surf-1)*2+0][1]};
				lineTag4Surf.push_back(line4);
				ltc=lineTags2.back();
				
			}else if(sec > 0  && surf==0){
				ptsv=combinedpts[surf];
				tr1 = getValuesAtIndices(ptsv,{0,1,2});
				tr2 = getValuesAtIndices(ptsv,{0,2,3});

				lineTags1 = createGmshLines(tr1, ltc,{0});
				lineTagsAll.push_back(lineTags1);
				line4={-lineTag4Surf[(sec-1)*(nodesec1+2)*2+surf*2+1][1],lineTags1[0],lineTags1[1]};;
				lineTag4Surf.push_back(line4);
				ltc=lineTags1.back();

				lineTags2 = createGmshLines(tr2, ltc,{0});
				lineTagsAll.push_back(lineTags2);
				line4={-lineTags1[1],lineTags2[0],lineTags2[1]};
				lineTag4Surf.push_back(line4);
				ltc=lineTags2.back();
				 
			}else if(sec > 0 && surf<nodesec1+1){
				
				ptsv=combinedpts[surf];
				tr1 = getValuesAtIndices(ptsv,{0,1,2});
				tr2 = getValuesAtIndices(ptsv,{0,2,3});

				lineTags1 = createGmshLines(tr1, ltc,{0});
				lineTagsAll.push_back(lineTags1);
				line4={-lineTag4Surf[(sec-1)*(nodesec1+2)*2+surf*2+1][1],lineTags1[0],lineTags1[1]};;
				lineTag4Surf.push_back(line4);
				ltc=lineTags1.back();

				lineTags2 = createGmshLines(tr2, ltc,{0,2});
				lineTagsAll.push_back(lineTags2);
				line4={-lineTags1[1],lineTags2[0],-lineTag4Surf[(sec)*(nodesec1+2)*2+(surf-1)*2+0][1]};
				lineTag4Surf.push_back(line4);
				ltc=lineTags2.back();
				
			}else if(sec > 0  && surf==nodesec1+1){
				ptsv=combinedpts[surf];
				tr1 = getValuesAtIndices(ptsv,{0,1,2});
				tr2 = getValuesAtIndices(ptsv,{0,2,3});

				lineTags1 = createGmshLines(tr1, ltc,{0,1});
				lineTagsAll.push_back(lineTags1);
				line4={-lineTag4Surf[(sec-1)*(nodesec1+2)*2+surf*2+1][1],-lineTag4Surf[(sec)*(nodesec1+2)*2+surf*0+1][2],lineTags1[0]};;
				lineTag4Surf.push_back(line4);
				ltc=lineTags1.back();

				lineTags2 = createGmshLines(tr2, ltc,{0,2});
				lineTagsAll.push_back(lineTags2);
				line4={-lineTags1[0],lineTags2[0],-lineTag4Surf[(sec)*(nodesec1+2)*2+(surf-1)*2+0][1]};
				lineTag4Surf.push_back(line4);
				ltc=lineTags2.back();
			}
 		}
	}
	std::cout<<ltc<<std::endl;
	std::cout<<"### Surface Tags:"<<std::endl;
	printVectorOfVectors(lineTag4Surf);
	
    //////////////////
    // Creating Surfaces 
    //////////////////
		/*
	for (const auto& group : lineTag4Surf)
	{
		// Create surface using the four point tags in the group
		createGmshSurface(group,stc);SurfTagSecs.push_back(stc);stc++;
	}


	std::vector<int> s1 = getEvenElementsAtIndex(lineTag4Surf,0, 0, nodesec1*2+2);
	std::cout<<"s1: "<<std::endl;
	printVector(s1);
    std::reverse(s1.begin(), s1.end()); // Reverse the copy
	createGmshSurface(s1,stc);SurfTagSecs.push_back(stc);stc++;

	
	std::vector<int> s2 = getOddElementsAtIndex(lineTag4Surf,1, lineTag4Surf.size()-(nodesec1*2+3), lineTag4Surf.size()+1);
	std::cout<<"s2: "<<std::endl;
	printVector(s2);	
	createGmshSurface(s2,stc);SurfTagSecs.push_back(stc);stc++;
	*/
    gmsh::model::geo::synchronize();
	

	
    //////////////////
    // Creating Solid
    //////////////////
    
    // Create a surface loop
    //int surfaceLoop = gmsh::model::geo::addSurfaceLoop(SurfTagSecs);

    // Create a volume from the surface loop
    //int volume = gmsh::model::geo::addVolume({surfaceLoop});	

}


int main(int argc, char **argv)
{ 
    
    //////////////////
    // Creating points
    //////////////////
    
    std::vector<double> Y_all = {	1.0, 1.0,1.0, 
									1.0, 1.0,1.0};
    int nodesec = 3, nX=2;
    double Xmin = 0;
    double Xmax = 1.5;	
    double Zmin = 0.0;
    double Zmax = 28.5;
    double Yzero = -0.1;
    

	// Initialize gmsh
    gmsh::initialize();

    // Create a new model
    gmsh::model::add("myModel");

    // Disable terminal outputs (optional)
    gmsh::option::setNumber("General.Terminal", 1);
    
	generateTile(Xmin, Xmax, nX, Zmin, Zmax, Yzero, Y_all, nodesec );
    //////////////////
    // Meshing
    //////////////////
    // Generate the mesh
    gmsh::model::geo::synchronize();
    
    //gmsh::model::mesh::generate(3);
    
	//generateMesh();
    //////////////////
    // Running the GUI
    //////////////////
    gmsh::model::geo::synchronize();
	std::set<std::string> args(argv, argv + argc);
	if(!args.count("-nopopup")) gmsh::fltk::run();
	gmsh::finalize();

    // Generate the mesh (optional)
    //gmsh::model::mesh::generate();



    // Finalize gmsh
    gmsh::finalize();

  /*
  int Znode=2;
  double ptsY[Znode+1];
  for(int i =0;i<3;i++){
	  ptsY[i]=3.;
  }
  double Ztot=28.5;
  double Xtot=3;
  */

    
    
    
 // GmshLYSO *LYSOmesh = new GmshLYSO( Znode, Xtot, Ztot, ptsY, "testname");

// GUI Plotting
  /*std::set<std::string> args(argv, argv + argc);
  if(!args.count("-nopopup")) gmsh::fltk::run();
  gmsh::finalize();*/



  return 0;
}
