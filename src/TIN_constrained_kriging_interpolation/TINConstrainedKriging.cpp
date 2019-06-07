#include "StdAfx.h"
#include "TINConstrainedKriging.h"
#include <fstream>
#include <sstream>

//set fault parameters, including center point O {ox, oy, oz}.
//and its three vectors {vs, vf, vd} (aligned on the fault strike, dip and normal direction, respectively)
void FaultGeometry::SetPara(Vector3d centerPoint, Vector3d vs, Vector3d vf, Vector3d vd
							,vector<Vector3d>& observationPoints)
{
	m_centerPoint = centerPoint;
	m_vs = vs;
	m_vf = vf;
	m_vd = vd;
	m_observePoints = observationPoints;
	GetMediumPlan();
}
//get a medium plan
void FaultGeometry::GetMediumPlan()
{
	//m_vd[0]*(x-x0)+m_vd[1]*(x-y0)+m_vd[2]*(x-z0)=0
	m_a = m_vd[0];
	m_b = m_vd[1];
	m_c = m_vd[2];
	m_d = -1 * m_vd.transpose()*m_centerPoint;
	m_vm.row(0).template segment<3>(0) = m_vs;
	m_vm.row(1).template segment<3>(0) = m_vf;
	m_vm.row(2).template segment<3>(0) = m_vd;
	m_vm_inverse = m_vm.inverse();
}
double FaultGeometry::GetSignedDistanceToMediumPlan(Vector3d point)
{
	return double(point.transpose()*m_vf) / double(sqrt(m_vf.transpose()*m_vf));
}
//construct fault surface
void FaultGeometry::InterpolateSurface(double step)
{
	//construct interpolator
	vector<Vector3d>::iterator it = m_observePoints.begin();
	for (it;it!=m_observePoints.end();++it)
	{
		Vector3d p = m_vm * (*it-m_centerPoint);
		Vector3d pChange(p[0], p[2] ,p[1]);
		m_observePointsLocal.push_back(pChange);
	}
	Variogram vgm(Sph,22,1000,2);
	okinterpolator.SetVariogramModel(vgm);
	okinterpolator.SetSamplePoints(m_observePointsLocal);
	//
	std::ofstream file("F:/points.txt");
	for (it=m_observePointsLocal.begin();it!=m_observePointsLocal.end();++it)
	{
		file << (*it)[0] << "," << (*it)[1] << ","<< (*it)[2] << std::endl;
	}
	file.close();
	//
	//compute grid in plane
	double min = -1 * m_halfextend;
	int iSize = static_cast<int>(m_halfextend * 2 / step) + 1;
	for (int i(0); i<iSize ; ++i)
	{
		double gs = min + step * i;
		for (int j(0); j<iSize ; ++j)
		{
			double gd = min + step * j;
			Triangle tri(3);
			tri[0] = Vector3d(gs,0.0,gd);
			tri[1] = Vector3d(gs+step,0.0,gd);
			tri[2] = Vector3d(gs+step,0.0,gd+step);
			faultTris.push_back(tri);
			//
			tri[0] = Vector3d(gs,0.0,gd);
			tri[1] = Vector3d(gs+step,0.0,gd+step);
			tri[2] = Vector3d(gs,0.0,gd+step);
			faultTris.push_back(tri);
		}
	}
	//interpolate and change coordinate
	vector<Triangle>::iterator itT = faultTris.begin();
	for (itT;itT!=faultTris.end();++itT)
		for (int i(0);i<3;++i)
			(*itT)[i] = LocalToGlobal((*itT)[i] );
}
//Translate coordinate
Vector3d FaultGeometry::LocalToGlobal(Vector3d point)
{
	//interpolate grid using the signed distance
	Vector3d pChange(point[0], point[2] ,point[1]);
	okinterpolator.EstimateVauleOfLocation(pChange);
	//
	Vector3d gCoord;
	gCoord = m_vm_inverse*(point-Vector3d(0,pChange[2],0))
			+ m_centerPoint;
	return gCoord;
}
Vector3d FaultGeometry::GlobalToLocal(Vector3d point)
{
	Vector3d lCoord = m_vm * (point-m_centerPoint);
	//interpolate grid using the signed distance
	Vector3d pChange(lCoord[0], lCoord[2] ,lCoord[1]);
	okinterpolator.EstimateVauleOfLocation(pChange);
	//
	lCoord[1] += pChange[2];
	return lCoord;
}
//get gd value of fault cutoff line in the local coordinate system
double FaultGeometry::GetGdOnCutofflines(double gs, double drop, bool upOrDown)
{
	double gd(0.0);
	if (fabs(gs) < m_halfextend)
	{
		if (upOrDown)//x^2 + (y+d)^2 = R^2, gd > 0
			gd = sqrt(pow(25.0,2)-pow(gs,2))-(25.0-drop);
		else//x^2 + (y-d)^2 = R^2, gd < 0
			gd = -1*sqrt(pow(25.0,2)-pow(gs,2))+(25.0-drop);
	}
	return gd;
}
//get hanging wall and footwall of fault (at one geological layer or geological layer)
void FaultGeometry::GetCutofflines(Vector3d centerAtCoal)
{
	pointsHanging.clear();
	pointsFoot.clear();
	//compute at local system
	//
	Vector3d point(-1*m_halfextend,0,0);
	for (int i(0);i<11;++i)
	{
		point[2] = GetGdOnCutofflines(point[0], m_drop, true);
		pointsHanging.push_back(point);
		point[2] = GetGdOnCutofflines(point[0], m_drop, false);
		pointsFoot.push_back(point);
		//
		point[0] += (m_halfextend/5);
	}
	//covert to global system
	vector<Vector3d>::iterator it = pointsHanging.begin();
	for (it;it!=pointsHanging.end();it++)
	{
		*it = LocalToGlobal(*it);
		pointsCutoffPolygon.push_back(*it);
	}
	vector<Vector3d>::reverse_iterator itr = pointsFoot.rbegin();
	for (itr; itr != pointsFoot.rend(); itr++)
	{
		*itr = LocalToGlobal(*itr);
		if (itr!=pointsFoot.rbegin() && itr!= (pointsFoot.rend()-1))
			pointsCutoffPolygon.push_back(*itr);
	}
}
//get gd value of fault deformation in the local coordinate system
double FaultGeometry::GetGdInLocalSystem(double gs, double dropGs)
{
	double gd(0);
	if (fabs(gs) < 2*fabs(dropGs))
	{
		gd = dropGs*(1-0.5/dropGs*gs);
	}
	return gd;
}
//Deform points in the fault local system
void FaultGeometry::Deformate(vector<Vector3d>& points)
{
	vector<Vector3d>::iterator it = points.begin();
	for (it;it!=points.end();++it)
	{
		Vector3d pointLocal = GlobalToLocal(*it);
		double dropGs(-0.0);
		if (pointLocal[2] > 0)
			dropGs = GetGdOnCutofflines(pointLocal[0],m_drop,true);
		else
			dropGs = GetGdOnCutofflines(pointLocal[0],m_drop,false);
		if (abs(dropGs) > 0)
		{
			double move = GetGdInLocalSystem(pointLocal[2],dropGs);
			if (fabs(move) > 0)
			{
				pointLocal[2] += move;
				*it = LocalToGlobal(pointLocal);
			}
		}
	}
}
//Get normals in the fault local system
vector<Vector3d> FaultGeometry::GetNormals(vector<Vector3d>& points)
{
	vector<Vector3d> normals;
	vector<Vector3d>::iterator it = points.begin();
	for (it;it!=points.end();++it)
	{
		Vector3d pointLocal = GlobalToLocal(*it);
		double dropGs(-0.0);
		if (pointLocal[2] > 0)
			dropGs = GetGdOnCutofflines(pointLocal[0],m_drop,true);
		else
			dropGs = GetGdOnCutofflines(pointLocal[0],m_drop,false);
		if (abs(dropGs) > 0)
		{
			double move = GetGdInLocalSystem(pointLocal[2],dropGs);
			if (fabs(move) > 0)
			{
				Vector3d deviation(pointLocal[0], pointLocal[2], pointLocal[1]);
				if (pointLocal[2] > 0)
					deviation = Vector3d(-2*deviation[0],-2*(deviation[1]+(25.0-m_drop)),0.0);
				else
					deviation = Vector3d(2*deviation[0],2*(deviation[1]-(25.0-m_drop)),0.0);
				//
				Matrix3d roate;
				roate <<  1.0, 0.0, 0.0
					, 0.0, 2.0/sqrt(5.0),-1.0/sqrt(5.0)
					, 0.0, 1.0/sqrt(5.0), 2.0/sqrt(5.0);
				deviation = roate*deviation;
				//
				deviation = LocalToGlobal(Vector3d(deviation[0],deviation[2],deviation[1]));
				deviation.normalize();
				normals.push_back(deviation);
				continue;
			}
		}
		normals.push_back(Vector3d(0.0,0.0,1.0));
	}
	return normals;
}
void FaultGeometry::Set(Vector3d centerPoint, Vector3d vs, Vector3d vf, Vector3d vd
						,vector<Vector3d>& observationPoints)
{
	m_isNormal = true;
	SetPara(centerPoint,vs,vf,vd,observationPoints);
	//
	m_drop = 10.0;
	m_halfextend = sqrt(pow(25.0,2)-pow(25.0-m_drop,2));
	InterpolateSurface(5.0);
	GetCutofflines(m_centerPoint);
}
//construct a interpolator
void GeologicalLayerGeometry::ConstructInterpolator(InterpolatorType interpolator)
{
	Variogram vgm(Sph,8,50,0.2);
	switch (interpolator)
	{
	case OrdinaryKriging:
		{
			okinterpolator.SetVariogramModel(vgm);
			okinterpolator.SetSamplePoints(m_samplePoints);
		}
		break;
	case ConstrainedKriging:
		{
			ckinterpolator.SetVariogramModel(vgm);
			ckinterpolator.SetSamplePoints(m_samplePoints);
			ckinterpolator.SetConstraintLines(m_faultsCutofflines);
			ckinterpolator.ConstructSearchingNetwork();
		}
		break;
	default:
		break;
	}
}
//scale positions of two points based on there center position 
void ScalePoints(Vector3d& p1,Vector3d& p2, float scale)
{
	Vector3d cen = p1 + ( p2 - p1 ) / 2;
	p1 = cen + ( p1 - cen ) * scale ;
	p2 = cen + ( p2 - cen ) * scale ;
}
//interpolate geological layer with sample points (and fault cutoff line points)
void GeologicalLayerGeometry::Interpolate(InterpolatorType interpolator)
{
	//create a 20*20 grid with in x[-50,50],y[-50,50] 
	//and the step of both of x and y axis are 5
	vector<Vector3d> gridPoints;
	double step = m_visualizedGridSize;
	double min = -50.0;
	int num = 25;
	for (int i(0);i<num;i++)
	{
		double gridx = min+i*step;
		for (int j(0);j<num;j++)
		{
			Vector3d v(gridx, min+j*step, 0.0);
			if (interpolator == Interpolator::OrdinaryKriging)
				okinterpolator.EstimateVauleOfLocation(v);
			else if(interpolator == Interpolator::ConstrainedKriging)
				ckinterpolator.EstimateVauleOfLocation(v);
			gridPoints.push_back(v);
		}
	}
	//insert faults into surface 
	m_interpolatedSurfaceTriangles = CovertTIN2Tris(CreateCDTIN(gridPoints,m_faultsCutofflines));
}

//visualize the geological layer surface (whole process)
void FaultedSurfaceModel::VisualizeLayer(GeologicalLayerGeometry& geologicalLayer)
{
	//create a 20*20 grid with in x[-50,50],y[-50,50] 
	//and the step of both of them are 5 
	//around x = 0, we have a fault there, so add one more line of data
	vector<Vector3d> gridPoints;
	double step = 5.0;
	double min = -50.0;
	int num = 20;
	for (int i(0);i<num;i++)
	{
		double gridx = min+i*step;
		if (i == num/2)
		{
			for (int j(0);j<num;j++)
				gridPoints.push_back(Vector3d(gridx-0.01, min+j*step, 0.0));
			for (int j(0);j<num;j++)
				gridPoints.push_back(Vector3d(gridx+0.01, min+j*step, 0.0));
		}else
		{
			for (int j(0);j<num;j++)
				gridPoints.push_back(Vector3d(gridx, min+j*step, 0.0));
		}
	}
	//deformate all points by fault
	fault.Deformate(gridPoints);
	//trianglization
	for (int i(0);i<num ;i++)
	{
		for (int j(0);j<num-1;j++)
		{
			//ignore triangles with x coordinate including both -0.01 and 0.01 
			//and located on fault
			Vector3d vp = fault.GlobalToLocal(gridPoints[i*num+j+1]);
			if (vp[1] > 0.0 && vp[1] < 0.1 && abs(vp[2]) > 0.01)
				continue;
			vp = fault.GlobalToLocal(gridPoints[(i+1)*num+j]);
			if (vp[1] < 0.0 && vp[1] > -0.1 && abs(vp[2]) > 0.01)
				continue;
			//
			Triangle tri(3);
			tri[2] = gridPoints[i*num+j];
			tri[1] = gridPoints[i*num+j+1];
			tri[0] = gridPoints[(i+1)*num+j];
			geologicalLayer.m_surfaceTriangles.push_back(tri);
			//
			tri[2] = gridPoints[i*num+j+1];
			tri[1] = gridPoints[(i+1)*num+j+1];
			tri[0] = gridPoints[(i+1)*num+j];
			geologicalLayer.m_surfaceTriangles.push_back(tri);
		}
	}
}

//cross validation
void FaultedSurfaceModel::CrossValidation(Interpolator::InterpolatorType interpolator)
{
	vector<float> stat_aver;
	vector<float> stat_stan;
	//
	vector<Vector3d> samplePoints = geologicalLayer.m_samplePoints;
	int numAll = (int)samplePoints.size();
	int num = int(ceil(numAll/10.0));
	std::vector<int> indexs(numAll);
	std::vector<int>::iterator it = indexs.begin();
	for (it;it!=indexs.end();it++)
		*it = it-indexs.begin();
	for (unsigned int i0(0);i0<10;++i0) {
		std::vector<pair<int,double>> predZ;
		//separate the data set into 10 subset
		srand( 0 );
		std::random_shuffle(indexs.begin(),indexs.end());
		std::vector<std::vector<int>> sub_indexs_all;
		for (unsigned int i(0);i<10;++i)
		{
			std::vector<int> sub_indexs;
			int iStart,iEnd;
			iStart = i*num;
			iEnd = iStart + num;
			sub_indexs.insert(sub_indexs.end(), indexs.begin() + iStart, indexs.begin() + iEnd);
			sort(sub_indexs.begin(),sub_indexs.end());
			sub_indexs_all.push_back(sub_indexs);
		}
		KrigingValidation(interpolator,samplePoints,sub_indexs_all,predZ);
		WriteCVResults(interpolator,i0,samplePoints,predZ);
		predZ.clear();
	}
	geologicalLayer.SetSamplePoints(samplePoints);
}
void FaultedSurfaceModel::KrigingValidation(Interpolator::InterpolatorType interpolator
											,vector<Vector3d> samplePoints
											,vector<vector<int>> &sub_indexs_all
											,vector<pair<int,double>> &predZ_vec2)
{
	vector<vector<int>>::iterator it;
	for (int i(0); i<10 ; ++i)
	{
		vector<Vector3d> coalDataPoint;
		vector<Vector3d> coalValidPoint;
		it = sub_indexs_all.begin();
		for (;it!=sub_indexs_all.end();++it)
		{
			vector<int>::iterator it2 = (*it).begin();
			if (it-sub_indexs_all.begin() == i)
			{
				for (;it2!=(*it).end();it2++)
					coalValidPoint.push_back(samplePoints[*it2]);
			}
			else
			{
				for (;it2!=(*it).end();it2++)
					coalDataPoint.push_back(samplePoints[*it2]);
			}
		}
		//
		geologicalLayer.SetSamplePoints(coalDataPoint);
		if (interpolator == Interpolator::ConstrainedKriging)
		{
			geologicalLayer.ConstructInterpolator(Interpolator::ConstrainedKriging);
		}
		else
		{
			geologicalLayer.ConstructInterpolator(Interpolator::OrdinaryKriging);
		}
		//
		vector<Vector3d>::iterator it3 = coalValidPoint.begin();
		for (;it3!=coalValidPoint.end();++it3)
		{
			Vector3d& point = *it3;
			if (interpolator == Interpolator::ConstrainedKriging)
			{
				geologicalLayer.ckinterpolator.EstimateVauleOfLocation(point);
			}
			else
			{
				geologicalLayer.okinterpolator.EstimateVauleOfLocation(point);
			}
			predZ_vec2.push_back(std::make_pair(sub_indexs_all[i][it3-coalValidPoint.begin()],point[2]));
		}
	}
}
void FaultedSurfaceModel::WriteCVResults(Interpolator::InterpolatorType interpolator
										 ,unsigned int i
										 ,vector<Vector3d> samplePoints
										 ,std::vector<pair<int,double>> predZ)
{
	std::stringstream ss;
	std::string file_name;
	std::ofstream fileOut;
	if (interpolator == Interpolator::OrdinaryKriging)
	{
		ss << "../../exe/crossValidation/CV_OK_" << i << ".csv";
	}else if (interpolator == Interpolator::ConstrainedKriging)
	{
		ss << "../../exe/crossValidation/CV_CK_" << i << ".csv";
	}
	else
		return;
	//
	file_name = ss.str();
	//write into file
	fileOut.open(file_name.c_str());
	fileOut << "x,y,z,predicted_z,difference" << std::endl;
	int numAll = (int)samplePoints.size();
	for (int j=0; j<numAll; j++){
		Vector3d &pt = samplePoints[predZ[j].first];
		fileOut << pt[0] << ',' << pt[1] << ',' << pt[2] << ',' << predZ[j].second << ',' << fabs(pt[2]-predZ[j].second);
		fileOut << std::endl;
	}
	fileOut.close();
}
