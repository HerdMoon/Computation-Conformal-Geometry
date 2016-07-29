// hw1.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <list>
#include <vector>
#include <math.h>
#include <Mesh\BaseMesh.h>
#include <Mesh\Face.h>
#include <Mesh\Vertex.h>
#include <Mesh\iterators.h>
#include <Geometry\Point.h>
#include <string>

using namespace MeshLib;
using namespace std;

typedef CHalfEdge * tHalfEdge;
typedef vector<CPoint> vec_pt;



int main()
{
	CBaseMesh<CVertex, CEdge, CFace, CHalfEdge> * BM = new(CBaseMesh<CVertex,CEdge,CFace,CHalfEdge>);
	BM->read_obj("e:\\model\\block.obj");
	int i;
	double sum = 0;
	for (i = 1;i <= BM->numVertices();i++)
	{
		vec_pt con_norm;
		con_norm.clear();

		VertexOutHalfedgeIterator<CVertex, CEdge, CFace, CHalfEdge>  it(BM, BM->idVertex(i));
		while (!it.end())
		{
			tHalfEdge now_he = it.value();
			CVertex  * v1 = now_he->target();
			CVertex  * v2 = now_he->source();
			CPoint pt1 = -v1->point() + v2->point();


			con_norm.push_back(pt1/pt1.norm());
			it++;
		}
		double  r = 0;
		for (int j = 0; j < con_norm.size() - 1; j++)
		{
			r = r + acos(con_norm[j] * con_norm[j + 1]);
		}
		r = r + acos(con_norm[0] * con_norm[con_norm.size()-1]);

		sum = sum + 2 * 3.1415926 - r;
	}
	printf("%f\n", sum);
	printf("\n%d", BM->numVertices()-BM->numEdges()+BM->numFaces());
    return 0;
}

