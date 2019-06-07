// TIN_constrained_kriging_interpolation.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "TINConstrainedKriging.h"
#include "SimulateData.h"
#include "FileOperator.h"

typedef vector<Vector3d> Triangle;
typedef vector<Triangle> Surface;

int _tmain(int argc, _TCHAR* argv[])
{
	unsigned int sampleSizeTest = 50; //100;500;1000;3000;5000;7000;9000;10000;20000;30000;40000;50000;
	//simulate data
	CSimulateData sd;
	CSimulateData::LayerSamplePoints layerData = sd.GetLayerSamples(sampleSizeTest);
	CSimulateData::FaultData faultData = sd.DefineFaultData();

	//set data
	FaultedSurfaceModel Modelor;
	Modelor.fault.Set(faultData.centerPoint,faultData.vs,faultData.vf,faultData.vd,faultData.faultObservations);
	Modelor.fault.Deformate(layerData.samplePoints);
	Modelor.geologicalLayer.SetSamplePoints(layerData.samplePoints);
	Modelor.geologicalLayer.SetFaultCutofflines(vector<vector<Vector3d>>(1,Modelor.fault.pointsCutoffPolygon));
	//Modelor.VisualizeLayer(Simulator.geologicalLayer);
	
	//interpolation
	//1:Constrained kriging
	Modelor.geologicalLayer.ConstructInterpolator(Interpolator::ConstrainedKriging);
	Modelor.geologicalLayer.Interpolate(Interpolator::ConstrainedKriging);
	
	//2:Ordinary kriging
	/*Modelor.geologicalLayer.ConstructInterpolator(Interpolator::OrdinaryKriging);
	Modelor.geologicalLayer.Interpolate(Interpolator::OrdinaryKriging);
	//Modelor.geologicalLayer.Intersecting(Simulator.fault);
	//Modelor.geologicalLayer.Deformate(Simulator.fault);
	*/
	//write faulted geological surface to .stl file
	CFileOperator fo;
	fo.WriteSTLFile(Modelor.geologicalLayer.m_interpolatedSurfaceTriangles);

	//cross validation
	Modelor.CrossValidation();
	return 0;
}

