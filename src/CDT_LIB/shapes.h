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

// Include guard
#ifndef SHAPES_H
#define SHAPES_H

#include <vector>
#include <cstddef>
#include <assert.h>
#include <cmath>
#include <Eigen\Dense>
using Eigen::Vector3d;

//namespace p2t {

struct TINEdge;

struct TINPoint {

  //double x, y, z;
	Vector3d xyz;

  /// Default constructor does nothing (for performance).
  TINPoint()
  {
    xyz.setZero();
  }

  /// The edges this point constitutes an upper ending point
  std::vector<TINEdge*> edge_list;

  /// Construct using coordinates.
  TINPoint(double x, double y, double z) : xyz(Vector3d(x, y, z)) {}

  /// Set this point to all zeros.
  void set_zero()
  {
	  xyz.setZero();
  }

  /// Set this point to some specified coordinates.
  void set(double x_, double y_, double z_)
  {
    xyz[0] = x_;
	xyz[1] = y_;
	xyz[2] = z_;
  }

  /// Negate this point.
  /*TINPoint operator -() const
  {
    TINPoint v;
    v.set(-xyz[0], -xyz[1], -xyz[2]);
    return v;
  }

  /// Add a point to this point.
  void operator +=(const TINPoint& v)
  {
    xyz[0] += v.xyz[0] ;
	xyz[1] += v.xyz[1] ;
	xyz[2] += v.xyz[2] ;
  }

  /// Subtract a point from this point.
  void operator -=(const TINPoint& v)
  {
    xyz[0] -= v.xyz[0] ;
	xyz[1] -= v.xyz[1] ;
	xyz[2] -= v.xyz[2] ;
  }*/

  /// Multiply this point by a scalar.
//   void operator *=(double a)
//   {
//     xyz[0] *= a;
// 	xyz[1] *= a;
// 	xyz[2] *= a;
//   }

  /// Get the length of this point (the norm).
  double Length() const
  {
    return sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1]);
  }

  /// Convert this point into a unit point. Returns the Length.
//   double Normalize()
//   {
//     double len = Length();
//     xyz[0] /= len;
//     xyz[1] /= len;
//     return len;
//   }

};

// Represents a simple polygon's edge
struct TINEdge {

  TINPoint* p, *q;

  /// Constructor
  TINEdge(TINPoint& p1, TINPoint& p2) : p(&p1), q(&p2)
  {
    if (p1.xyz[1] > p2.xyz[1] ) {
      q = &p1;
      p = &p2;
    } else if (p1.xyz[1] == p2.xyz[1] ) {
      if (p1.xyz[0] > p2.xyz[0] ) {
        q = &p1;
        p = &p2;
      } else if (p1.xyz[0] == p2.xyz[0] ) {
        // Repeat points
        assert(false);
      }
    }

    q->edge_list.push_back(this);
  }
};

// TINTriangle-based data structures are know to have better performance than quad-edge structures
// See: J. Shewchuk, "TINTriangle: Engineering a 2D Quality Mesh Generator and Delaunay Triangulator"
//      "Triangulations in CGAL"
struct TINTriangle {
public:

/// Constructor
TINTriangle(TINPoint& a, TINPoint& b, TINPoint& c);

/// Flags to determine if an edge is a Constrained edge
bool constrained_edge[3];
/// Flags to determine if an edge is a Delauney edge
bool delaunay_edge[3];

TINPoint* GetPoint(const int& index);
TINPoint* PointCW(TINPoint& point);
TINPoint* PointCCW(TINPoint& point);
TINPoint* OppositePoint(TINTriangle& t, TINPoint& p);

TINTriangle* GetNeighbor(const int& index);
void MarkNeighbor(TINPoint* p1, TINPoint* p2, TINTriangle* t);
void MarkNeighbor(TINTriangle& t);

void MarkConstrainedEdge(const int index);
void MarkConstrainedEdge(TINEdge& edge);
void MarkConstrainedEdge(TINPoint* p, TINPoint* q);

int Index(const TINPoint* p);
int EdgeIndex(const TINPoint* p1, const TINPoint* p2);

TINTriangle* NeighborCW(TINPoint& point);
TINTriangle* NeighborCCW(TINPoint& point);
bool GetConstrainedEdgeCCW(TINPoint& p);
bool GetConstrainedEdgeCW(TINPoint& p);
void SetConstrainedEdgeCCW(TINPoint& p, bool ce);
void SetConstrainedEdgeCW(TINPoint& p, bool ce);
bool GetDelunayEdgeCCW(TINPoint& p);
bool GetDelunayEdgeCW(TINPoint& p);
void SetDelunayEdgeCCW(TINPoint& p, bool e);
void SetDelunayEdgeCW(TINPoint& p, bool e);

bool Contains(TINPoint* p);
bool Contains(const TINEdge& e);
bool Contains(TINPoint* p, TINPoint* q);
void Legalize(TINPoint& point);
void Legalize(TINPoint& opoint, TINPoint& npoint);
/**
 * Clears all references to all other triangles and points
 */
void Clear();
void ClearNeighbor(TINTriangle *triangle );
void ClearNeighbors();
void ClearDelunayEdges();

inline bool IsInterior();
inline void IsInterior(bool b);

TINTriangle& NeighborAcross(TINPoint& opoint);

void DebugPrint();

private:

/// TINTriangle points
TINPoint* points_[3];
/// Neighbor list
TINTriangle* neighbors_[3];

/// Has this triangle been marked as an interior triangle?
bool interior_;
};

inline bool cmp(const TINPoint* a, const TINPoint* b)
{
  if (a->xyz[1] < b->xyz[1]) {
    return true;
  } else if (a->xyz[1] == b->xyz[1]) {
    // Make sure q is point with greater xyz[0] value
    if (a->xyz[0] < b->xyz[0]) {
      return true;
    }
  }
  return false;
}

/// Add two points_ component-wise.
inline TINPoint operator +(const TINPoint& a, const TINPoint& b)
{
  return TINPoint(a.xyz[0] + b.xyz[0] , a.xyz[1] + b.xyz[1] , a.xyz[2] + b.xyz[2] );
}

/// Subtract two points_ component-wise.
inline TINPoint operator -(const TINPoint& a, const TINPoint& b)
{
  return TINPoint(a.xyz[0] - b.xyz[0] , a.xyz[1] - b.xyz[1] , a.xyz[2] - b.xyz[2] );
}

/// Multiply point by scalar
inline TINPoint operator *(double s, const TINPoint& a)
{
  return TINPoint(s * a.xyz[0] , s * a.xyz[1] , s * a.xyz[2] );
}

inline bool operator ==(const TINPoint& a, const TINPoint& b)
{
  return a.xyz[0] == b.xyz[0] && a.xyz[1] == b.xyz[1] ;
}

inline bool operator !=(const TINPoint& a, const TINPoint& b)
{
  return !(a.xyz[0] == b.xyz[0] ) && !(a.xyz[1] == b.xyz[1] );
}

/// Peform the dot product on two vectors.
inline double Dot(const TINPoint& a, const TINPoint& b)
{
  return a.xyz[0] * b.xyz[0] + a.xyz[1] * b.xyz[1] ;
}


inline TINPoint* TINTriangle::GetPoint(const int& index)
{
  return points_[index];
}

inline TINTriangle* TINTriangle::GetNeighbor(const int& index)
{
  return neighbors_[index];
}

inline bool TINTriangle::Contains(TINPoint* p)
{
  return p == points_[0] || p == points_[1] || p == points_[2];
}

inline bool TINTriangle::Contains(const TINEdge& e)
{
  return Contains(e.p) && Contains(e.q);
}

inline bool TINTriangle::Contains(TINPoint* p, TINPoint* q)
{
  return Contains(p) && Contains(q);
}

inline bool TINTriangle::IsInterior()
{
  return interior_;
}

inline void TINTriangle::IsInterior(bool b)
{
  interior_ = b;
}

//}

#endif


