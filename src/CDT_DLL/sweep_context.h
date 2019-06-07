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

#ifndef SWEEP_CONTEXT_H
#define SWEEP_CONTEXT_H
#include "stdafx.h"
#include <list>
#include <vector>
#include <cstddef>

namespace p2t {

// Inital triangle factor, seed triangle will extend 30% of
// PointSet width to both left and right.
const double kAlpha = 0.3;

//struct TINPoint;
//class TINTriangle;
struct Node;
//struct TINEdge;
class AdvancingFront;

class SweepContext {
public:

/// Constructor
SweepContext(std::vector<TINPoint*> polyline);
/// Destructor
~SweepContext();

void set_head(TINPoint* p1);

TINPoint* head();

void set_tail(TINPoint* p1);

TINPoint* tail();

int point_count();

Node& LocateNode(TINPoint& point);

void RemoveNode(Node* node);

void CreateAdvancingFront(std::vector<Node*> nodes);

/// Try to map a node to all sides of this triangle that don't have a neighbor
void MapTriangleToNodes(TINTriangle& t);

void AddToMap(TINTriangle* triangle);

TINPoint* GetPoint(const int& index);

TINPoint* GetPoints();

void RemoveFromMap(TINTriangle* triangle);

void AddHole(std::vector<TINPoint*> polyline);

void AddPoint(TINPoint* point);

AdvancingFront* front();

void MeshClean(TINTriangle& triangle);

std::vector<TINTriangle*> GetTriangles();
std::list<TINTriangle*> GetMap();

std::vector<TINEdge*> edge_list;

struct Basin {
  Node* left_node;
  Node* bottom_node;
  Node* right_node;
  double width;
  bool left_highest;

  Basin() : left_node(NULL), bottom_node(NULL), right_node(NULL), width(0.0), left_highest(false)
  {
  }

  void Clear()
  {
    left_node = NULL;
    bottom_node = NULL;
    right_node = NULL;
    width = 0.0;
    left_highest = false;
  }
};

struct EdgeEvent {
  TINEdge* constrained_edge;
  bool right;

  EdgeEvent() : constrained_edge(NULL), right(false)
  {
  }
};

Basin basin;
EdgeEvent edge_event;

private:

friend class Sweep;

std::vector<TINTriangle*> triangles_;
std::list<TINTriangle*> map_;
std::vector<TINPoint*> points_;

// Advancing front
AdvancingFront* front_;
// head point used with advancing front
TINPoint* head_;
// tail point used with advancing front
TINPoint* tail_;

Node *af_head_, *af_middle_, *af_tail_;

void InitTriangulation();
void InitEdges(std::vector<TINPoint*> polyline);

};

inline AdvancingFront* SweepContext::front()
{
  return front_;
}

inline int SweepContext::point_count()
{
  return points_.size();
}

inline void SweepContext::set_head(TINPoint* p1)
{
  head_ = p1;
}

inline TINPoint* SweepContext::head()
{
  return head_;
}

inline void SweepContext::set_tail(TINPoint* p1)
{
  tail_ = p1;
}

inline TINPoint* SweepContext::tail()
{
  return tail_;
}

}

#endif
