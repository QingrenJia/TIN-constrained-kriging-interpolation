#pragma once
//STL
#include <vector>
using std::vector;
#include <fstream>
//Eigen
#include <Eigen/Dense>
using Eigen::Vector3d;
//
typedef vector<Vector3d> Triangle;
//
class CFileOperator
{
public:
	CFileOperator(void);
	~CFileOperator(void);
	//solid name
	//facet normal ni nj nk
	//	outer loop
	//	vertex v1x v1y v1z
	//	vertex v2x v2y v2z
	//	vertex v3x v3y v3z
	//	endloop
	//	endfacet
	//endsolid name
	void WriteSTLFile(const vector<Triangle> &triangles)
	{
		std::ofstream fileOut;
		std::string file_name = "../../exe/surfacemodel.stl";
		//write into file
		fileOut.open(file_name.c_str());
		fileOut << "solid geologicalSurface" << std::endl;
		vector<Triangle>::const_iterator it = triangles.begin();
		for (; it != triangles.end() ; it++)
		{
			fileOut << "facet" << std::endl;
			fileOut << "  outer loop" << std::endl;
			vector<Vector3d>::const_iterator itN = it->begin();
			for (; itN != it->end() ; itN++)
			{
					fileOut << "    vertex " << (*itN)[0] << ' ' << (*itN)[1] << ' ' << (*itN)[2] << std::endl;;
			}
			fileOut << "  endloop" << std::endl;
			fileOut << "endfacet" << std::endl;
		}
		fileOut << "endsolid geologicalSurface" << std::endl;
		fileOut.close();
	}
};
