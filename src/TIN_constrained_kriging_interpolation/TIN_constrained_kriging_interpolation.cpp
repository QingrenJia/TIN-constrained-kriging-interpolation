// TIN_constrained_kriging_interpolation.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "TINConstrainedKriging.h"
#include "SimulateData.h"
#include "FileOperator.h"
#include <iostream>

typedef vector<Vector3d> Triangle;
typedef vector<Triangle> Surface;

int _tmain(int argc, _TCHAR* argv[])
{
	unsigned int sampleSizeTest = 50; //100;500;1000;3000;5000;7000;9000;10000;20000;30000;40000;50000;
	std::cout << "------------------- Start -------------------" << std::endl;
	std::cout << std::endl << "Simulate simple points of geological layer" << std::endl;
	std::cout << "--Please input point size:";
	std::cin >> sampleSizeTest;
	if (0 == sampleSizeTest)
		return 0;
	//simulate data
	CSimulateData sd;
	CSimulateData::LayerSamplePoints layerData = sd.GetLayerSamples(sampleSizeTest);
	CSimulateData::FaultData faultData = sd.DefineFaultData();
	
	std::cout << "--Simpled "<< sampleSizeTest << " points in range [minx,miny,maxx,maxy] = [-50, -50, 50, 50]" << std::endl;
	std::cout << std::endl << "Simpled fault data:" << std::endl;
	std::cout << "--CenterPoint:" << faultData.centerPoint[0]<< ',' << faultData.centerPoint[1]<< ',' << faultData.centerPoint[2]<< ',' << std::endl;
	std::cout << "--Direction vector aligned on the fault strike:" << faultData.vs[0]<< ',' << faultData.vs[1]<< ',' << faultData.vs[2]<< ',' << std::endl;
	std::cout << "--Direction vector of dip:" << faultData.vf[0]<< ',' << faultData.vf[1]<< ',' << faultData.vf[2]<< ',' << std::endl;
	std::cout << "--Direction vector of normal:" << faultData.vd[0]<< ',' << faultData.vd[1]<< ',' << faultData.vd[2]<< ',' << std::endl;
	std::cout << "--Fault observations:"<< std::endl;
	for (vector<Vector3d>::iterator it=faultData.faultObservations.begin(); it!=faultData.faultObservations.end();++it)
		std::cout << "  " << (*it)[0]<< ',' << (*it)[1]<< ',' << (*it)[2]<< ',' << std::endl;
	//set data
	FaultedSurfaceModel Modelor;
	Modelor.fault.Set(faultData.centerPoint,faultData.vs,faultData.vf,faultData.vd,faultData.faultObservations);
	Modelor.fault.Deformate(layerData.samplePoints);
	Modelor.geologicalLayer.SetSamplePoints(layerData.samplePoints);
	Modelor.geologicalLayer.SetFaultCutofflines(vector<vector<Vector3d>>(1,Modelor.fault.pointsCutoffPolygon));
	//Modelor.VisualizeLayer(Simulator.geologicalLayer);
	
	//interpolation
	char op[CHAR_MAX];
	std::cout << std::endl << "Please select your operation:" << std::endl;
	std::cout << "--IN: Apply ordinary-kriging / TIN-constrained-kriging method to dataset" << std::endl;
	std::cout << "--CV: Apply cross validation" << std::endl;
	std::cin >> op;
	if ( (op[0] == 'i'|| op[0] == 'I')
		&&(op[1] == 'n'|| op[1] == 'N'))
	{
		//1:Constrained kriging
		Modelor.geologicalLayer.ConstructInterpolator(Interpolator::ConstrainedKriging);
		Modelor.geologicalLayer.Interpolate(Interpolator::ConstrainedKriging);
		//write faulted geological surface to .stl file
		CFileOperator fo;
		fo.WriteSTLFile(Interpolator::ConstrainedKriging,Modelor.geologicalLayer.m_interpolatedSurfaceTriangles);
		
		//2:Ordinary kriging
		Modelor.geologicalLayer.ConstructInterpolator(Interpolator::OrdinaryKriging);
		Modelor.geologicalLayer.Interpolate(Interpolator::OrdinaryKriging);
		//write faulted geological surface to .stl file
		fo.WriteSTLFile(Interpolator::OrdinaryKriging,Modelor.geologicalLayer.m_interpolatedSurfaceTriangles);
		
	}else if ( (op[0] == 'c'|| op[0] == 'C')
		&&(op[1] == 'v'|| op[1] == 'V'))
	{
		//cross validation
		Modelor.CrossValidation(Interpolator::OrdinaryKriging);
		Modelor.CrossValidation(Interpolator::ConstrainedKriging);
	}
	else
	{
		std::cout << "You can only choose either interpolation or cross validation operation" << std::endl;
	}
	std::cout << "Files outputted to exe folder." << std::endl;
	std::cout << "Press enter to quit." << std::endl;
	system("pause");
	return 0;
}

