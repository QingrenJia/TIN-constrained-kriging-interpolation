#include "StdAfx.h"
#include "CDTAlgorithm.h"

CCDTAlgorithm::CCDTAlgorithm(void)
{
}

CCDTAlgorithm::~CCDTAlgorithm(void)
{
}
bool PointInArray(TINPoint* p,vector<TINPoint*> &parray)
{
	vector<TINPoint*>::iterator it = parray.begin();
	for (it;it!=parray.end();++it)
		if (fabs((*it)->xyz[0]-p->xyz[0])<0.001 && fabs((*it)->xyz[1]-p->xyz[1])<0.001)
			return true;
	return false;

}
CDT* CCDTAlgorithm::Entrance(vector<TINPoint> &boundingPoly
							 ,vector<TINPoint> &steinerPoints
							 ,vector<vector<TINPoint>> &head_holes
							 ,vector<vector<TINPoint>> &chest_holes)
{
	vector<TINPoint*> tinPoints;
	//bounding box
	vector<TINPoint*> polyline;
	vector<TINPoint>::iterator it = boundingPoly.begin();
	for (it;it!=boundingPoly.end();++it)
		polyline.push_back(new TINPoint(it->xyz[0],it->xyz[1],it->xyz[2]));
	polylines.push_back(polyline);
	/* Perform triangulation!*/
	/*
	* STEP 1: Create CDT and add primary polyline
	* NOTE: polyline must be a simple polygon. The polyline's points
	* constitute constrained edges. No repeat points!!!
	*/
	CDT* cdt = new CDT(polyline);

	/* STEP 2: Add holes or Steiner points if necessary */
	// Add head hole
	vector<vector<TINPoint>>::iterator itH = head_holes.begin();
	for (itH;itH!=head_holes.end();++itH)
	{
		vector<TINPoint*> head_hole;
		for (it=itH->begin();it!=itH->end();++it)
			head_hole.push_back(new TINPoint(it->xyz[0],it->xyz[1],it->xyz[2]));
		polylines.push_back(head_hole);
		cdt->AddHole(head_hole);
		tinPoints.insert(tinPoints.end(),head_hole.begin(),head_hole.end());
	}
	// Add chest hole
	itH = chest_holes.begin();
	for (itH;itH!=chest_holes.end();++itH)
	{
		vector<TINPoint*> chest_hole;
		for (it=itH->begin();it!=itH->end();++it)
			chest_hole.push_back(new TINPoint(it->xyz[0],it->xyz[1],it->xyz[2]));
		polylines.push_back(chest_hole);
		cdt->AddHole(chest_hole);
		tinPoints.insert(tinPoints.end(),chest_hole.begin(),chest_hole.end());
	}
	polyline.clear();
	it = steinerPoints.begin();
	for (it;it!=steinerPoints.end();++it)
	{
		polyline.push_back(new TINPoint(it->xyz[0],it->xyz[1],it->xyz[2]));
		cdt->AddPoint(polyline.back());
		tinPoints.push_back(polyline.back());
	}
	polylines.push_back(polyline);
	/* STEP 3: Triangulate! */
	cdt->Triangulate();

	//triangles = cdt->GetTriangles();
	map = cdt->GetMap();
	
	//// Cleanup
	//delete cdt;
	//// Free points
	//for(int i = 0; i < polylines.size(); i++) {
	//	vector<TINPoint*> poly = polylines[i];
	//	FreeClear(poly);
	//}
	return cdt;
}
