#pragma once
#include<windows.h>
//STL
#include <set>
#include <vector>
#include <list>
using std::set;
using std::pair;
using std::vector;
using std::list;
//Eigen
#include <Eigen/Dense>
using Eigen::Vector3d;
using Eigen::Matrix3d;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::LDLT;
//CD-TIN
#include "CDT/include/shapes.h"
#pragma  comment(lib, "CDT_LIB.lib")
//
typedef vector<Vector3d> Triangle;
typedef vector<Triangle> Surface;
//
namespace Interpolator
{
	enum InterpolatorType{ OrdinaryKriging = 0, ConstrainedKriging };
	enum VariogramModel{Sph = 0, Exp, Gau};
	//create a TIN with the constraint of holes
	static vector<TINTriangle*> CreateCDTIN(vector<Vector3d> &samplePoint,vector<vector<Vector3d>> &constraintLines)
	{
		typedef vector<TINTriangle*> (*CallDllFun1)(vector<Vector3d>* ,
			vector<Triangle>* ,//head holes
			vector<Triangle>* );//chest holes
		HINSTANCE htin = LoadLibrary(_T("CDT_DLL.dll"));
		if(!htin) 
			return vector<TINTriangle*>();
		CallDllFun1 CreatTin=(CallDllFun1)GetProcAddress(htin,"ConstructConstrainedDTIN");
		vector<TINTriangle*> interTinTriangles = CreatTin(&samplePoint,NULL,&constraintLines);
		FreeLibrary(htin);
		return interTinTriangles;
	}
	// return all triangles of TIN
	static vector<Triangle> CovertTIN2Tris(vector<TINTriangle*>& TIN)
	{
		vector<Triangle> tris;
		vector<TINTriangle*>::iterator itPoly = TIN.begin();
		for(itPoly;itPoly!=TIN.end();++itPoly) {
			Triangle tri;
			tri.push_back((*itPoly)->GetPoint(0)->xyz);
			tri.push_back((*itPoly)->GetPoint(1)->xyz);
			tri.push_back((*itPoly)->GetPoint(2)->xyz);
			tris.push_back(tri);
		}
		return tris;
	}
	struct Variogram
	{
		Variogram():e(2.718281828){}
		Variogram(VariogramModel type,double sill,double range, double nugget):e(2.718281828)
		{
			this->m_sill = sill;
			this->m_range = range;
			this->m_nugget = nugget;
			this->m_modelType = type;
		}
		double m_sill;
		double m_range;
		double m_nugget;
		VariogramModel m_modelType;
		//
		const double e;
		Variogram& operator = (Variogram &model)
		{
			this->m_sill = model.m_sill;
			this->m_range = model.m_range;
			this->m_nugget = model.m_nugget;
			this->m_modelType = model.m_modelType;
			return *this;
		}
		double GetGamma(const double d)
		{
			double dist = fabs(d);
			double gamma;
			if (dist < 0.001)
				return 0.0;
			else if( dist >= m_range )
				return m_sill;
			else 
			{
				switch(m_modelType) 
				{
				case Sph: 
					{
						double d1 = dist/m_range*3/2;
						double d2 =	pow(dist/m_range,3)/2;
						gamma = m_nugget+((m_sill-m_nugget)*(d1-d2));
					}
					break;
				case Exp:
					gamma =  m_sill *(1-pow(e,(-3*dist/m_range)));
					break;
				case Gau: 
					{
						//(h/a)^2 ; r = sqrt(3)*a
						double d = pow(3*dist/m_range,2);
						gamma = m_nugget + m_sill * (1-pow(e,-d));
					}
					break;
				default:
					dist = 0;
					break;
				}
			}
			return gamma;
		}
	};
	//Ordinary kriging
	class OKInterpolator
	{
	public:
		OKInterpolator(void){};
		//
		Variogram m_variogram;
		void SetVariogramModel(Variogram model)
		{
			m_variogram = model;
		}
		// Data to be interpolated
		vector<Vector3d>	m_samplePoints;
		void SetSamplePoints(vector<Vector3d>& points)
		{
			m_samplePoints = points;
		}
		//give an estimation of value at unsampled location
		double EstimateVauleOfLocation(Vector3d& loc)
		{
			set<Vector3d*> input;
			SelNeighborsByDist(input,loc);
			return EstimateValue(input,loc);
		}
		//estimation
		double EstimateValue(set<Vector3d*>& input,Vector3d& loc)
		{
			vector<Vector3d> nearPoints;
			set<Vector3d*>::iterator it;
			for (it = input.begin();it!=input.end();++it)
				nearPoints.push_back(**it);
			//
			MatrixXd matK;
			MatrixXd matM;
			MatrixXd matWeight;
			int nSize = nearPoints.size() + 1;
			//Constructing a covariance matrix
			matK.resize(nSize, nSize);
			for(int i=0; i<nSize-1; i++)
			{
				matK(i, nSize-1) = 1.0;
				matK(nSize-1, i) = 1.0;
			}
			matK(nSize-1, nSize-1) = 0.0;
			for(int i=0; i < nSize-1; i++) 
			{
				for(int j=0; j < nSize-1; j++) 
				{
					if (i == j)
					{
						matK(i, j) = m_variogram.m_sill;
						continue;
					}
					Vector3d& point_i = nearPoints.at(i);
					Vector3d& point_j = nearPoints.at(j);
					double dis = GetPointsDistance_2D(point_i[0],point_i[1],point_j[0],point_j[1]);
					double g = m_variogram.GetGamma(dis);
					matK(i, j) = m_variogram.m_sill - g;
				}
			}
			//Constructing the right matrix of the ok interpolation equation
			matM.resize(nSize,1);
			matM(nSize-1,0)= 1.0;
			for (int i = 0;i<nSize-1;i++)
			{
				Vector3d& temp = nearPoints.at(i);
				double dis = GetPointsDistance_2D(loc[0],loc[1],temp[0],temp[1]);
				double tempGamma = m_variogram.GetGamma(dis);
				matM(i,0) = m_variogram.m_sill - tempGamma;
			}
			//solve
			loc[2] = 0.0;
			matWeight.resize(nSize,1);
			LDLT<MatrixXd> llt;
			llt.compute(matK);
			matWeight = llt.solve(matM);
			//prediction values
			for (int i=0;i<nSize-1;i++)
			{
				double w = matWeight(i,0);
				loc[2] += nearPoints[i][2] * w;
			}
			//prediction variances
			MatrixXd M = matWeight.transpose()*matM;
			return M(0,0) - m_variogram.m_nugget;
		}
	private:
		static bool Compare(pair<double,Vector3d*>& lhs, pair<double,Vector3d*>& rhs)
		{
			return lhs.first < rhs.first;
		}
		double GetPointsDistance_2D(double xpos1, double ypos1, double xpos2, double ypos2)
		{
			return sqrt(pow((xpos1-xpos2), 2) + pow((ypos1-ypos2), 2));
		}
		//select 8 neighbors according to their distance to the location of point to be estimated
		void SelNeighborsByDist(set<Vector3d*>& nearPoint,Vector3d &loc)
		{
			vector<Vector3d>::iterator it = m_samplePoints.begin();
			nearPoint.clear();
			vector<pair<double,Vector3d*>> dist_pt;
			for (it;it!=m_samplePoints.end();it++)
			{
				Vector3d& temp = *it;
				double dis = sqrt(pow(temp[0] - loc[0],2)+pow(temp[1] - loc[1],2));
				dist_pt.push_back(pair<double,Vector3d*>(dis,&temp));
			}
			sort(dist_pt.begin(),dist_pt.end(),Compare);
			for (int i=0 ; i < 8; ++i )
			{
				Vector3d* temp = dist_pt[i].second;
				nearPoint.insert(temp);
			}
		}
	};
	//TIN constrained kriging
	//TIN is constructed by CDT algorithm using points and constraint lines, such as fault cutoff lines
	//, and is used to constraint the neighborhood selection process of OK method
	class CKInterpolator : public OKInterpolator
	{
	public:
		// constraint lines, such as fault cutoff lines, all these lines forming a closed polygon
		vector<vector<Vector3d>> m_constraintLines;
		//
		vector<TINTriangle*> m_interTinPolys;
		//
		void SetConstraintLines(vector<vector<Vector3d>> &constraintlines)
		{
			m_constraintLines = constraintlines;
		}
		// construct searching network, i.e., a TIN serving the neighborhood selection process
		void ConstructSearchingNetwork()
		{
			m_interTinPolys = CreateCDTIN(m_samplePoints,m_constraintLines);
		}
		//give an estimation of value at unsampled location
		double EstimateVauleOfLocation(Vector3d& newPoint)
		{
			set<Vector3d*> input;
			SelNeighborsInTIN(input,newPoint);
			return EstimateValue(input,newPoint);
		}
	private:
		// judge if triangle contains point to be estimated
		bool PointInTriangle(const Vector3d &loc, TINTriangle *tri)
		{
			Vector3d pointA(tri->GetPoint(0)->xyz);
			Vector3d pointB(tri->GetPoint(1)->xyz);
			Vector3d pointC(tri->GetPoint(2)->xyz);
			Vector3d PA = pointA - loc;
			Vector3d PB = pointB - loc;
			Vector3d PC = pointC - loc;
			double t1 = PA[0]*PB[1] - PA[1]*PB[0];
			double t2 = PB[0]*PC[1] - PB[1]*PC[0];
			double t3 = PC[0]*PA[1] - PC[1]*PA[0];
			return t1*t2 >= 0 && t1*t3 >= 0;
		}
		void SelNeighborsInTIN(set<Vector3d*>& nearPoint,Vector3d& point)
		{
			TINTriangle* pTri(NULL);
			set<TINPoint*> nearTinPoint;
			vector<TINTriangle*>::iterator itPoly = m_interTinPolys.begin();
			for(itPoly;itPoly!=m_interTinPolys.end();++itPoly) {
				pTri = *itPoly;
				//
				double minx = min(pTri->GetPoint(0)->xyz[0],pTri->GetPoint(1)->xyz[0]);
				minx = min(minx,pTri->GetPoint(2)->xyz[0]);
				double maxx = max(pTri->GetPoint(0)->xyz[0],pTri->GetPoint(1)->xyz[0]);
				maxx = max(maxx,pTri->GetPoint(2)->xyz[0]);
				double miny = min(pTri->GetPoint(0)->xyz[1],pTri->GetPoint(1)->xyz[1]);
				miny = min(miny,pTri->GetPoint(2)->xyz[1]);
				double maxy = max(pTri->GetPoint(0)->xyz[1],pTri->GetPoint(1)->xyz[1]);
				maxy = max(maxy,pTri->GetPoint(2)->xyz[1]);
				if (point[0]<minx || point[0]>maxx
					|| point[1]<miny || point[1]>maxy)
					continue;
				//
				if (PointInTriangle(point,pTri))
					break;
				pTri = NULL;
			}
			if (!pTri){
				return;
			}
			// sample points from pTri's nodes
			nearTinPoint.insert(pTri->GetPoint(0));
			nearTinPoint.insert(pTri->GetPoint(1));
			nearTinPoint.insert(pTri->GetPoint(2));
			// sample points from pTri's neighbors
			if (TINTriangle* neighbor0 = pTri->GetNeighbor(0))
			{
				nearTinPoint.insert(neighbor0->GetPoint(0));
				nearTinPoint.insert(neighbor0->GetPoint(1));
				nearTinPoint.insert(neighbor0->GetPoint(2));
			}
			if (TINTriangle* neighbor1 = pTri->GetNeighbor(1))
			{
				nearTinPoint.insert(neighbor1->GetPoint(0));
				nearTinPoint.insert(neighbor1->GetPoint(1));
				nearTinPoint.insert(neighbor1->GetPoint(2));
			}
			if (TINTriangle* neighbor2 = pTri->GetNeighbor(1))
			{
				nearTinPoint.insert(neighbor2->GetPoint(0));
				nearTinPoint.insert(neighbor2->GetPoint(1));
				nearTinPoint.insert(neighbor2->GetPoint(2));
			}
			//
			set<TINPoint*>::iterator itP = nearTinPoint.begin();
			for(itP;itP!=nearTinPoint.end();++itP) {
				nearPoint.insert(&((*itP)->xyz));
			}
		}
	};
}