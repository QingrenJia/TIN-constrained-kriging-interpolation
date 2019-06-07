#pragma once

//STL
#include <vector>
using std::vector;
//Eigen
#include <Eigen/Dense>
using Eigen::Vector3d;

#define Range_min -50.0
#define Range_max 50.0

class CSimulateData
{
public:
	CSimulateData(void);
	~CSimulateData(void);
	struct LayerSamplePoints{
		vector<Vector3d> samplePoints;
	};
	struct FaultData{
		Vector3d centerPoint;
		Vector3d vs;
		Vector3d vf;
		Vector3d vd;
		vector<Vector3d> faultObservations;
	};
	//simulate samples for geological layer
	LayerSamplePoints GetLayerSamples(unsigned int simpleCount)
	{
		LayerSamplePoints lp;
		srand(0);
		lp.samplePoints.reserve(simpleCount);
		for (unsigned int i(0);i<simpleCount;++i)
		{
			Vector3d sample(0.0f,0.0f,0.0f);
			sample[0] = (double)rand() / (RAND_MAX + 1) *(Range_max - Range_min) + Range_min;
			sample[1] = (double)rand() / (RAND_MAX + 1) *(Range_max - Range_min) + Range_min;
			lp.samplePoints.push_back(sample);
		}
		return lp;
	}
	//define simulate data for fault
	FaultData DefineFaultData()
	{
		FaultData fd;
		fd.centerPoint = Vector3d (0.0,0.0,0.0);
		fd.vs = Vector3d ( 0.0  ,1.0, 0.0  );
		fd.vf = Vector3d (-0.866,0.0, 0.5  );
		fd.vd = Vector3d (-0.5  ,0.0,-0.866);
		fd.faultObservations.push_back(Vector3d(-6.044765081,9.723434391,-10.46953312));
		fd.faultObservations.push_back(Vector3d(-5.375518375,-39.61537091,-9.310397825));
		fd.faultObservations.push_back(Vector3d(-14.72559852,22.50759897,-25.50473663));
		fd.faultObservations.push_back(Vector3d(11.71057711,-27.21068041,20.28271956));
		fd.faultObservations.push_back(Vector3d(3.28475248,15.47642223,5.689191295));
		fd.faultObservations.push_back(Vector3d(-19.29997759,18.84648705,-33.42756119));
		fd.faultObservations.push_back(Vector3d(-8.664725274,-49.11408414,-15.00730418));
		fd.faultObservations.push_back(Vector3d(-0.07806791,47.63417789,-0.135213619));
		fd.faultObservations.push_back(Vector3d(16.82668148,-12.35379372,29.14381233));
		fd.faultObservations.push_back(Vector3d(2.781830667,-30.88659123,4.818130715));
		return fd;
	}
};
