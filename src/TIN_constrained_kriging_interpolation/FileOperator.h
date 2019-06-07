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
	void WriteSTLFile(unsigned short method,const vector<Triangle> &triangles)
	{
		std::ofstream fileOut;
		std::string file_name;
		if (method == 0)
			file_name = "../../exe/surfaceModel/surfacemodel_ok.stl";
		else if (method == 1)
			file_name = "../../exe/surfaceModel/surfacemodel_ck.stl";
		else
			return;
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
