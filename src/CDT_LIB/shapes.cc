#include "stdafx.h"
/* 
 * Poly2Tri Copyright (c) 2009-2010, Poly2Tri Contributors
 * http://code.google.com/p/poly2tri/
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice,
 *   this list of conditions and the following disclaimer.
 * * Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 * * Neither the name of Poly2Tri nor the names of its contributors may be
 *   used to endorse or promote products derived from this software without specific
 *   prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#include "shapes.h"
#include <iostream>


TINTriangle::TINTriangle(TINPoint& a, TINPoint& b, TINPoint& c)
{
  points_[0] = &a; points_[1] = &b; points_[2] = &c;
  neighbors_[0] = NULL; neighbors_[1] = NULL; neighbors_[2] = NULL;
  constrained_edge[0] = constrained_edge[1] = constrained_edge[2] = false;
  delaunay_edge[0] = delaunay_edge[1] = delaunay_edge[2] = false;
  interior_ = false;
}

// Update neighbor pointers
void TINTriangle::MarkNeighbor(TINPoint* p1, TINPoint* p2, TINTriangle* t)
{
  if ((p1 == points_[2] && p2 == points_[1]) || (p1 == points_[1] && p2 == points_[2]))
    neighbors_[0] = t;
  else if ((p1 == points_[0] && p2 == points_[2]) || (p1 == points_[2] && p2 == points_[0]))
    neighbors_[1] = t;
  else if ((p1 == points_[0] && p2 == points_[1]) || (p1 == points_[1] && p2 == points_[0]))
    neighbors_[2] = t;
  else
    assert(0);
}

// Exhaustive search to update neighbor pointers
void TINTriangle::MarkNeighbor(TINTriangle& t)
{
  if (t.Contains(points_[1], points_[2])) {
    neighbors_[0] = &t;
    t.MarkNeighbor(points_[1], points_[2], this);
  } else if (t.Contains(points_[0], points_[2])) {
    neighbors_[1] = &t;
    t.MarkNeighbor(points_[0], points_[2], this);
  } else if (t.Contains(points_[0], points_[1])) {
    neighbors_[2] = &t;
    t.MarkNeighbor(points_[0], points_[1], this);
  }
}

/**
 * Clears all references to all other triangles and points
 */
void TINTriangle::Clear()
{
    TINTriangle *t;
    for( int i=0; i<3; i++ )
    {
        t = neighbors_[i];
        if( t != NULL )
        {
            t->ClearNeighbor( this );
        }
    }
    ClearNeighbors();
    points_[0]=points_[1]=points_[2] = NULL;
}

void TINTriangle::ClearNeighbor(TINTriangle *triangle )
{
    if( neighbors_[0] == triangle )
    {
        neighbors_[0] = NULL;
    }
    else if( neighbors_[1] == triangle )
    {
        neighbors_[1] = NULL;            
    }
    else
    {
        neighbors_[2] = NULL;
    }
}
    
void TINTriangle::ClearNeighbors()
{
  neighbors_[0] = NULL;
  neighbors_[1] = NULL;
  neighbors_[2] = NULL;
}

void TINTriangle::ClearDelunayEdges()
{
  delaunay_edge[0] = delaunay_edge[1] = delaunay_edge[2] = false;
}

TINPoint* TINTriangle::OppositePoint(TINTriangle& t, TINPoint& p)
{
  TINPoint *cw = t.PointCW(p);
  double x = cw->xyz[0];
  double y = cw->xyz[1];
  x = p.xyz[0];
  y = p.xyz[1];
  return PointCW(*cw);
}

// Legalized triangle by rotating clockwise around point(0)
void TINTriangle::Legalize(TINPoint& point)
{
  points_[1] = points_[0];
  points_[0] = points_[2];
  points_[2] = &point;
}

// Legalize triagnle by rotating clockwise around oPoint
void TINTriangle::Legalize(TINPoint& opoint, TINPoint& npoint)
{
  if (&opoint == points_[0]) {
    points_[1] = points_[0];
    points_[0] = points_[2];
    points_[2] = &npoint;
  } else if (&opoint == points_[1]) {
    points_[2] = points_[1];
    points_[1] = points_[0];
    points_[0] = &npoint;
  } else if (&opoint == points_[2]) {
    points_[0] = points_[2];
    points_[2] = points_[1];
    points_[1] = &npoint;
  } else {
    assert(0);
  }
}

int TINTriangle::Index(const TINPoint* p)
{
  if (p == points_[0]) {
    return 0;
  } else if (p == points_[1]) {
    return 1;
  } else if (p == points_[2]) {
    return 2;
  }
  assert(0);
}

int TINTriangle::EdgeIndex(const TINPoint* p1, const TINPoint* p2)
{
  if (points_[0] == p1) {
    if (points_[1] == p2) {
      return 2;
    } else if (points_[2] == p2) {
      return 1;
    }
  } else if (points_[1] == p1) {
    if (points_[2] == p2) {
      return 0;
    } else if (points_[0] == p2) {
      return 2;
    }
  } else if (points_[2] == p1) {
    if (points_[0] == p2) {
      return 1;
    } else if (points_[1] == p2) {
      return 0;
    }
  }
  return -1;
}

void TINTriangle::MarkConstrainedEdge(const int index)
{
  constrained_edge[index] = true;
}

void TINTriangle::MarkConstrainedEdge(TINEdge& edge)
{
  MarkConstrainedEdge(edge.p, edge.q);
}

// Mark edge as constrained
void TINTriangle::MarkConstrainedEdge(TINPoint* p, TINPoint* q)
{
  if ((q == points_[0] && p == points_[1]) || (q == points_[1] && p == points_[0])) {
    constrained_edge[2] = true;
  } else if ((q == points_[0] && p == points_[2]) || (q == points_[2] && p == points_[0])) {
    constrained_edge[1] = true;
  } else if ((q == points_[1] && p == points_[2]) || (q == points_[2] && p == points_[1])) {
    constrained_edge[0] = true;
  }
}

// The point counter-clockwise to given point
TINPoint* TINTriangle::PointCW(TINPoint& point)
{
  if (&point == points_[0]) {
    return points_[2];
  } else if (&point == points_[1]) {
    return points_[0];
  } else if (&point == points_[2]) {
    return points_[1];
  }
  assert(0);
}

// The point counter-clockwise to given point
TINPoint* TINTriangle::PointCCW(TINPoint& point)
{
  if (&point == points_[0]) {
    return points_[1];
  } else if (&point == points_[1]) {
    return points_[2];
  } else if (&point == points_[2]) {
    return points_[0];
  }
  assert(0);
}

// The neighbor clockwise to given point
TINTriangle* TINTriangle::NeighborCW(TINPoint& point)
{
  if (&point == points_[0]) {
    return neighbors_[1];
  } else if (&point == points_[1]) {
    return neighbors_[2];
  }
  return neighbors_[0];
}

// The neighbor counter-clockwise to given point
TINTriangle* TINTriangle::NeighborCCW(TINPoint& point)
{
  if (&point == points_[0]) {
    return neighbors_[2];
  } else if (&point == points_[1]) {
    return neighbors_[0];
  }
  return neighbors_[1];
}

bool TINTriangle::GetConstrainedEdgeCCW(TINPoint& p)
{
  if (&p == points_[0]) {
    return constrained_edge[2];
  } else if (&p == points_[1]) {
    return constrained_edge[0];
  }
  return constrained_edge[1];
}

bool TINTriangle::GetConstrainedEdgeCW(TINPoint& p)
{
  if (&p == points_[0]) {
    return constrained_edge[1];
  } else if (&p == points_[1]) {
    return constrained_edge[2];
  }
  return constrained_edge[0];
}

void TINTriangle::SetConstrainedEdgeCCW(TINPoint& p, bool ce)
{
  if (&p == points_[0]) {
    constrained_edge[2] = ce;
  } else if (&p == points_[1]) {
    constrained_edge[0] = ce;
  } else {
    constrained_edge[1] = ce;
  }
}

void TINTriangle::SetConstrainedEdgeCW(TINPoint& p, bool ce)
{
  if (&p == points_[0]) {
    constrained_edge[1] = ce;
  } else if (&p == points_[1]) {
    constrained_edge[2] = ce;
  } else {
    constrained_edge[0] = ce;
  }
}

bool TINTriangle::GetDelunayEdgeCCW(TINPoint& p)
{
  if (&p == points_[0]) {
    return delaunay_edge[2];
  } else if (&p == points_[1]) {
    return delaunay_edge[0];
  }
  return delaunay_edge[1];
}

bool TINTriangle::GetDelunayEdgeCW(TINPoint& p)
{
  if (&p == points_[0]) {
    return delaunay_edge[1];
  } else if (&p == points_[1]) {
    return delaunay_edge[2];
  }
  return delaunay_edge[0];
}

void TINTriangle::SetDelunayEdgeCCW(TINPoint& p, bool e)
{
  if (&p == points_[0]) {
    delaunay_edge[2] = e;
  } else if (&p == points_[1]) {
    delaunay_edge[0] = e;
  } else {
    delaunay_edge[1] = e;
  }
}

void TINTriangle::SetDelunayEdgeCW(TINPoint& p, bool e)
{
  if (&p == points_[0]) {
    delaunay_edge[1] = e;
  } else if (&p == points_[1]) {
    delaunay_edge[2] = e;
  } else {
    delaunay_edge[0] = e;
  }
}

// The neighbor across to given point
TINTriangle& TINTriangle::NeighborAcross(TINPoint& opoint)
{
  if (&opoint == points_[0]) {
    return *neighbors_[0];
  } else if (&opoint == points_[1]) {
    return *neighbors_[1];
  }
  return *neighbors_[2];
}

void TINTriangle::DebugPrint()
{
  using namespace std;
  cout << points_[0]->xyz[0] << "," << points_[0]->xyz[1] << " ";
  cout << points_[1]->xyz[0] << "," << points_[1]->xyz[1] << " ";
  cout << points_[2]->xyz[0] << "," << points_[2]->xyz[1] << endl;
}

