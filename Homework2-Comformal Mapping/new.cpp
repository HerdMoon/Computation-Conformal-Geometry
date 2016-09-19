#include "stdafx.h"

#include "OBJFileReader.h"
#include "Solid.h"
#include "SolidDelegate.h"
#include "iterators.h"
#include <math.h>
#include <map>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>

using namespace MeshLib;
using namespace std;

int main()
{
	const double delta_d = 0.001;
	Solid mesh;
	Solid target_mesh;
	SolidDelegate temp_mesh; //soliddelegate类封装了很多调整solid边与面关系的算法
	OBJFileReader obj_file;
	ifstream in("E:\\face_model\\3.obj");
	obj_file.readToSolid(&mesh, in);
	
	//初始化Harmonic Edge weight.
	map<Edge *, double> edge_kuv;
	for (SolidEdgeIterator e_iter(&mesh); !e_iter.end(); ++e_iter)
	{
		Edge * e = *e_iter;
		edge_kuv[e] = e->kuv();
	}

	//Weingarten Map
	for (SolidVertexIterator v_iter(&mesh); !v_iter.end(); ++v_iter)
	{
		Vertex *v = temp_mesh.createVertex(&target_mesh, target_mesh.numVertices() + 1);
		v->point() = (*v_iter)->normal();
		v->id() = (*v_iter)->id();
	}
	for (SolidFaceIterator f_iter(&mesh); !f_iter.end(); ++f_iter)
	{
		Face *f = *f_iter;
		HalfEdge *he = f->halfedge();
		int vertices[3];
		vertices[0] = he->source()->id();
		vertices[1] = he->he_next()->source()->id();
		vertices[2] = he->he_next()->he_next()->source()->id();
		int face[3] = { vertices[0], vertices[1], vertices[2] };
		temp_mesh.createFace(&target_mesh, face, target_mesh.numFaces() + 1);
	}
	
	//首次计算能量E0

	double E0;
	E0 = 0;

	for (SolidEdgeIterator e_iter(&mesh); !e_iter.end(); ++e_iter)
	{
		Edge * e = * e_iter;
		Vertex * v1;
		Vertex * v2;
		e->get_vertices(v1, v2);
		E0+= (v1->point() - v2->point()).norm2();		
	}

	Point center;
	//大循环.调整f(v)
	double newE = E0;
	Point temp_p;
	SolidVertexIterator l_iter(&mesh);
	SolidVertexIterator sv_iter(&target_mesh);
	SolidEdgeIterator ne_iter(&target_mesh);
	while (true)
	{
		target_mesh.UpdateNormals();
		l_iter.reset();
		//计算Laplace
		for (;!l_iter.end();++l_iter)
		{
			Vertex *v     = *l_iter;
			Vertex *new_v = target_mesh.idVertex(v->id());
			temp_p[0] = 0;
			temp_p[1] = 0;
			temp_p[2] = 0;
			for (VertexVertexIterator vv_iter(v); !vv_iter.end(); ++vv_iter)
			{
				Vertex *vv = *vv_iter;
				Vertex *newvv = target_mesh.idVertex(vv->id());
				temp_p[0] += newvv->point()[0] - new_v->point()[0];
				temp_p[1] += newvv->point()[1] - new_v->point()[1];
				temp_p[2] += newvv->point()[2] - new_v->point()[2];
			}
			temp_p /= temp_p.norm();
			new_v->Laplacian() = temp_p;
		}

		sv_iter.reset();
		for (; !sv_iter.end(); ++sv_iter)
		{
			Vertex *v = *sv_iter;
			Point D = -(v->Laplacian() - (v->point() * (v->Laplacian() * v->point()))) * delta_d;
			v->point() -= D;
			v->point() /= v->point().norm();
		}

		center[0] = 0;
		center[1] = 0;
		center[2] = 0;
		for (SolidVertexIterator iter_centre(&target_mesh); !iter_centre.end(); ++iter_centre)
		{
			Vertex *v = *iter_centre;
			center = center + v->point()*v->v_farea();
		}
		center = center / target_mesh.numVertices();
		for (SolidVertexIterator pv(&target_mesh); !pv.end(); ++pv)
		{
			Vertex *pV = *pv;
			pV->point() -= center;
			pV->point() /= pV->point().norm();

		}
		printf("%f\n", E0);
		ne_iter.reset();
		newE = 0; 
		double e0 = 0;
		for (; !ne_iter.end(); ++ne_iter)
		{
			Edge * e = *ne_iter;
			Vertex *v1, *v2;
			e->get_vertices(v1, v2);
			e0 = (v1->point() - v2->point()).norm2();
			newE = newE + e0;
		}
		if (fabs(E0 - newE) < 1e-5)
			break;
		E0 = newE;
	}
	ofstream out("E:\\face_model\\result.obj");
	obj_file.writeToSolid(&target_mesh, out);
}