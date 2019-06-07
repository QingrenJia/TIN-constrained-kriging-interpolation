#pragma once

#include <cstdlib>
#include <time.h>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <iostream>
using namespace std;
#include <Eigen\Dense>
using Eigen::Vector3d;
//
#include "poly2tri.h"
using namespace p2t;
class CCDTAlgorithm
{
public:
	CCDTAlgorithm(void);
	~CCDTAlgorithm(void);

	CDT* Entrance(vector<TINPoint> &boundingPoly
		,vector<TINPoint> &steinerPoints
		,vector<vector<TINPoint>> &head_holes
		,vector<vector<TINPoint>> &chest_holes);

	/// Dude hole examples
	vector<TINPoint*> CreateHeadHole();
	vector<TINPoint*> CreateChestHole();

	/// Constrained triangles
	vector<TINTriangle*> triangles;
	/// TINTriangle map
	list<TINTriangle*> map;
	/// Polylines
	vector< vector<TINPoint*> > polylines;

	template <class C> void FreeClear( C & cntr ) {
		for ( typename C::iterator it = cntr.begin(); 
				  it != cntr.end(); ++it ) {
			delete * it;
		}
		cntr.clear();
	}

	double CCDTAlgorithm::Fun(double x)
	{
		return 2.5 + sin(10 * x) / x;
	}

};
