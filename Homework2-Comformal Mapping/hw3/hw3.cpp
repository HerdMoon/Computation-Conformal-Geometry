// hw3.cpp : 定义控制台应用程序的入口点。
//
#include "stdafx.h"
#include <fstream>
#include <vector>

#include "OBJFileReader.h"
#include "Solid.h"
#include "iterators.h"

#include "SolidDelegate.h"
#include <math.h>
#include <map>
#include <iostream>

using namespace MeshLib;
using namespace std;
double t = 1e-5;
double steplength = 1e-3;

int r()
{

	Point center;

	E0 = E;
	SolidVertexIterator iterhSD(&newMesh);
	SolidVertexIterator iterhL(&mesh);
	SolidEdgeIterator e1_hiter(&newMesh);
	Edge * temp;

	//HARMONIC ENERGY CALCULATIONS
	while (1)
	{
		newMesh.UpdateNormals();


		iterhL.reset();
		for (; !iterhL.end(); ++iterhL)
		{
			Vertex *v = *iterhL;
			Vertex *newv = newMesh.idVertex(v->id());
			L[0] = 0; L[1] = 0; L[2] = 0;
			for (VertexVertexIterator vv_iter(v); !vv_iter.end(); ++vv_iter)
			{
				Vertex *vv = *vv_iter;
				Vertex *newvv = newMesh.idVertex(vv->id());
				L[0] += newvv->point()[0] - newv->point()[0];
				L[1] += newvv->point()[1] - newv->point()[1];
				L[2] += newvv->point()[2] - newv->point()[2];
			}
			//  if( L.norm() > 1e-5 )
			// {
			L /= L.norm();
			// }
			newv->Laplacian() = L;
		}
		iterhSD.reset();
		for (; !iterhSD.end(); ++iterhSD)
		{
			Vertex *v = *iterhSD;
			Point D = -(v->Laplacian() - (v->point() * (v->Laplacian() * v->point()))) * steplength;
			//D /= D.norm();
			v->point() -= D;
			v->point() /= v->point().norm();
		}

		//Harmonic mass center

		center[0] = 0; center[1] = 0; center[2] = 0;
		//STEP NUMBER : 4'-1
		for (SolidVertexIterator iter_centre(&newMesh); !iter_centre.end(); ++iter_centre)
		{
			Vertex *v = *iter_centre;
			center = center + v->point()*v->v_farea();
		}
		center = center / newMesh.numVertices();
		//STEP NUMBER : 4'- 2 & 4' - 3
		for (SolidVertexIterator pv(&newMesh); !pv.end(); ++pv)
		{
			Vertex *pV = *pv;
			pV->point() -= center;
			pV->point() /= pV->point().norm();

		}

		e1_hiter.reset();

		e0 = 0;
		//HARMONIC ENERGY ::
		for (; !e1_hiter.end(); ++e1_hiter)
		{
			Edge * e = *e1_hiter;
			Vertex *v1, *v2;
			e->get_vertices(v1, v2);

			Vertex *v1_new = mesh.idVertex(v1->id());
			Vertex *v2_new = mesh.idVertex(v2->id());

			temp = mesh.vertexEdge(v1_new, v2_new);

			e0 += ((v1->point() - v2->point()).norm2())*edgekuv[temp];//* (v1->point() - v2->point()).norm();

		}
		cout << "Harmonic energy is :" << e0 << endl;
		cout << "Energy Difference is: " << e0 - E0 << endl;

		if (fabs(e0 - E0)<1e-3)
			break;
		E0 = e0;

	}

	// Write out the resultant obj file
	int vObjID = 1;
	std::map<int, int> vidToObjID;

	std::ofstream os("e:\\face_model\\result.obj");

	SolidVertexIterator iter(&newMesh);

	for (; !iter.end(); ++iter)
	{
		Vertex *v = *iter;
		Point p = v->point();
		os << "v " << p[0] << " " << p[1] << " " << p[2] << std::endl;
		vidToObjID[v->id()] = vObjID++;
	}
	os << "# " << (unsigned int)newMesh.numVertices() << " vertices" << std::endl;

	float u = 0.0, v = 0.0;
	for (iter.reset(); !iter.end(); ++iter)
	{
		Vertex *vv = *iter;
		std::string key("uv");
		std::string s = Trait::getTraitValue(vv->string(), key);
		if (s.length() > 0)
		{
			sscanf_s(s.c_str(), "%f %f", &u, &v);
		}
		os << "vt " << u << " " << v << std::endl;
	}
	os << "# " << (unsigned int)newMesh.numVertices() << " texture coordinates" << std::endl;

	SolidFaceIterator fiter2(&newMesh);
	for (; !fiter2.end(); ++fiter2)
	{
		Face *f = *fiter2;
		FaceVertexIterator viter(f);
		os << "f ";
		for (; !viter.end(); ++viter)
		{
			Vertex *v = *viter;
			os << vidToObjID[v->id()] << "/" << vidToObjID[v->id()] << " ";
		}
		os << std::endl;
	}
	os.close();

	return 0;
}