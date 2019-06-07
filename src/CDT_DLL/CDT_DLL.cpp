// CDT_DLL.cpp : 定义 DLL 应用程序的导出函数。
//

#include "stdafx.h"
#include "CDTAlgorithm.h"
#include <vector>
#include <Eigen\Dense>
using std::vector;
using Eigen::Vector3d;


vector<TINTriangle*> ConstructConstrainedDTIN(vector<Vector3d>* points,
						   vector<vector<Vector3d>>* headholes,
						   vector<vector<Vector3d>>* chestholes)
{
	if (!points)
		return vector<TINTriangle*>();
	//formate data
	vector<TINPoint> boundingBox;
	vector<TINPoint> innerPoints;
	vector<vector<TINPoint>> head_holes;
	vector<vector<TINPoint>> chest_holes;
	vector<Vector3d>::iterator it = points->begin();
	double minx((*points)[0][0]),miny((*points)[0][1]),maxx((*points)[0][0]),maxy((*points)[0][1]);
	for (it;it!=points->end();++it)
	{
		innerPoints.push_back(TINPoint((*it)[0],(*it)[1],(*it)[2]));
		minx = min(minx,(*it)[0]);
		miny = min(miny,(*it)[1]);
		maxx = max(maxx,(*it)[0]);
		maxy = max(maxy,(*it)[1]);
	}
	boundingBox.push_back(TINPoint(minx-10,miny-10,0));
	boundingBox.push_back(TINPoint(maxx+10,miny-10,0));
	boundingBox.push_back(TINPoint(maxx+10,maxy+10,0));
	boundingBox.push_back(TINPoint(minx-10,maxy+10,0));
	//
	vector<vector<Vector3d>>::iterator itH;
	if (headholes)
	{
		itH = headholes->begin();
		for (itH;itH!=headholes->end();++itH)
		{
			vector<TINPoint> hole;
			it = itH->begin();
			for (it;it!=itH->end();++it)
			{
				hole.push_back(TINPoint((*it)[0],(*it)[1],(*it)[2]));
			}
			head_holes.push_back(hole);
		}
	}
	if (chestholes)
	{
		itH = chestholes->begin();
		for (itH;itH!=chestholes->end();++itH)
		{
			vector<TINPoint> hole;
			it = itH->begin();
			for (it;it!=itH->end();++it)
			{
				hole.push_back(TINPoint((*it)[0],(*it)[1],(*it)[2]));
			}
			chest_holes.push_back(hole);
		}
	}
	//construct TIN and insert holes
	CCDTAlgorithm Algorithm;
	CDT* cdt = Algorithm.Entrance(boundingBox, innerPoints, head_holes, chest_holes);
	return cdt->GetTriangles();
}
