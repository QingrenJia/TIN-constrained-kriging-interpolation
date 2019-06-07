#pragma once
#include "Interpolator.h"
using std::string;
using namespace Interpolator;

// to simulate and storage a fault geometry
class FaultGeometry
{
private:
	//get a medium plan
	void GetMediumPlan();
	//get distance to medium plan
	double GetSignedDistanceToMediumPlan(Vector3d point);
	//ax+by+cz+d=0
	double m_a,m_b,m_c,m_d;
	//
	Matrix3d m_vm,m_vm_inverse;
public:
	//set fault parameters, including center point O {ox, oy, oz}.
	//and its three vectors {vs, vf, vd} (aligned on the fault strike, dip and normal direction, respectively)
	//and fault observation points
	void SetPara(Vector3d centerPoint, Vector3d vs, Vector3d vf, Vector3d vd,vector<Vector3d>& observationPoints);
	//construct fault surface
	void InterpolateSurface(double step);
	//Translate coordinate
	Vector3d LocalToGlobal(Vector3d point);
	Vector3d GlobalToLocal(Vector3d point);
	//get gd value of fault cutoff line in the local coordinate system
	double GetGdOnCutofflines(double gs, double drop, bool upOrDown);
	//get hanging wall and footwall of fault (at one geological layer or geological layer)
	void GetCutofflines(Vector3d centerAtCoal);
	//get gd value of fault deformation in the local coordinate system
	double GetGdInLocalSystem(double gs, double dropGs);
	//Deform points in the fault local system
	void Deformate(vector<Vector3d>& points);
	//set a fault
	void Set(Vector3d centerPoint, Vector3d vs, Vector3d vf, Vector3d vd
		,vector<Vector3d>& observationPoints);
	//get normals in the fault local system
	vector<Vector3d> GetNormals(vector<Vector3d>& points);
	//center point O {ox, oy, oz}
	Vector3d m_centerPoint;
	//fault strike, dip and normal direction
	Vector3d m_vs,m_vf,m_vd;
	//half of extend of fault
	double m_halfextend;
	//fault drop
	double m_drop;
	//bool true = normal fault
	bool m_isNormal;
	//
	vector<Vector3d> m_observePoints;
	vector<Vector3d> m_observePointsLocal;
	//fault hangingwall and fFootwall lines on a geological layer
	vector<Vector3d> pointsHanging;
	vector<Vector3d> pointsFoot;
	vector<Vector3d> pointsCutoffPolygon;
	//fault surface triangles
	vector<Triangle> faultTris;
	//an ordinary kriging interpolator
	OKInterpolator okinterpolator;
};

class GeologicalLayerGeometry
{
public:
	GeologicalLayerGeometry(void):m_visualizedGridSize(4.0f){};
	~GeologicalLayerGeometry(void){};
	//data to be interpolated
	vector<Vector3d> m_samplePoints;
	vector<vector<Vector3d>> m_faultsCutofflines;
	//interpolator
	OKInterpolator okinterpolator;
	CKInterpolator ckinterpolator;
	//grid size for estimation and visualization
	float m_visualizedGridSize;
	//simulated surface model
	vector<Triangle> m_surfaceTriangles;
	//surface model built from interpolating sample points and cutoff lines
	vector<Triangle> m_interpolatedSurfaceTriangles;
	//set collected sample points of geological layer surface 
	void SetSamplePoints(vector<Vector3d>& samplePoints)
	{
		m_samplePoints = samplePoints;
	}
	//set computed intersection lines (fault cutoff lines) of geological layer surface and fault surface 
	void SetFaultCutofflines(vector<vector<Vector3d>>& faultsCutofflines)
	{
		m_faultsCutofflines = faultsCutofflines;
	}
	//construct a interpolator
	void ConstructInterpolator(InterpolatorType interpolator);
	//interpolate geological layer sample points (and fault cutoff line points)
	void Interpolate(InterpolatorType interpolator);
};

class FaultedSurfaceModel
{
public:
	FaultGeometry fault;
	GeologicalLayerGeometry geologicalLayer;
	void VisualizeLayer(GeologicalLayerGeometry& geologicalLayer);
	// operating a cross validation
	void CrossValidation(Interpolator::InterpolatorType interpolator);
	void KrigingValidation(Interpolator::InterpolatorType interpolator
		,vector<Vector3d> samplePoints
		,vector<vector<int>> &sub_indexs_all
		,vector<pair<int,double>> &predZ_vec2);
	void WriteCVResults(Interpolator::InterpolatorType interpolator
		,unsigned int i
		,vector<Vector3d> samplePoints
		,std::vector<pair<int,double>> predZ);
};
