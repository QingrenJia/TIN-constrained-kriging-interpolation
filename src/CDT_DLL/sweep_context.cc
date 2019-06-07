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
#include "sweep_context.h"
#include <algorithm>
#include "advancing_front.h"

namespace p2t {

SweepContext::SweepContext(std::vector<TINPoint*> polyline) :
  front_(0),
  head_(0),
  tail_(0),
  af_head_(0),
  af_middle_(0),
  af_tail_(0)
{
  basin = Basin();
  edge_event = EdgeEvent();

  points_ = polyline;

  InitEdges(points_);
}

void SweepContext::AddHole(std::vector<TINPoint*> polyline)
{
  InitEdges(polyline);
  for(unsigned int i = 0; i < polyline.size(); i++) {
    points_.push_back(polyline[i]);
  }
}

void SweepContext::AddPoint(TINPoint* point) {
  points_.push_back(point);
}

std::vector<TINTriangle*> SweepContext::GetTriangles()
{
  return triangles_;
}

std::list<TINTriangle*> SweepContext::GetMap()
{
  return map_;
}

void SweepContext::InitTriangulation()
{
  double xmax(points_[0]->xyz[0]), xmin(points_[0]->xyz[0]);
  double ymax(points_[0]->xyz[1]), ymin(points_[0]->xyz[1]);

  // Calculate bounds.
  for (unsigned int i = 0; i < points_.size(); i++) {
    TINPoint& p = *points_[i];
    if (p.xyz[0] > xmax)
      xmax = p.xyz[0];
    if (p.xyz[0] < xmin)
      xmin = p.xyz[0];
    if (p.xyz[1] > ymax)
      ymax = p.xyz[1];
    if (p.xyz[1] < ymin)
      ymin = p.xyz[1];
  }

  double dx = kAlpha * (xmax - xmin);
  double dy = kAlpha * (ymax - ymin);
  head_ = new TINPoint(xmax + dx, ymin - dy, 0.0);
  tail_ = new TINPoint(xmin - dx, ymin - dy, 0.0);

  // Sort points along y-axis
  std::sort(points_.begin(), points_.end(), cmp);

}

void SweepContext::InitEdges(std::vector<TINPoint*> polyline)
{
  int num_points = polyline.size();
  for (int i = 0; i < num_points; i++) {
    int j = i < num_points - 1 ? i + 1 : 0;
    edge_list.push_back(new TINEdge(*polyline[i], *polyline[j]));
  }
}

TINPoint* SweepContext::GetPoint(const int& index)
{
  return points_[index];
}

void SweepContext::AddToMap(TINTriangle* triangle)
{
  map_.push_back(triangle);
}

Node& SweepContext::LocateNode(TINPoint& point)
{
  // TODO implement search tree
  return *front_->LocateNode(point.xyz[0]);
}

void SweepContext::CreateAdvancingFront(std::vector<Node*> nodes)
{

  (void) nodes;
  // Initial triangle
  TINTriangle* triangle = new TINTriangle(*points_[0], *tail_, *head_);

  map_.push_back(triangle);

  af_head_ = new Node(*triangle->GetPoint(1), *triangle);
  af_middle_ = new Node(*triangle->GetPoint(0), *triangle);
  af_tail_ = new Node(*triangle->GetPoint(2));
  front_ = new AdvancingFront(*af_head_, *af_tail_);

  // TODO: More intuitive if head is middles next and not previous?
  //       so swap head and tail
  af_head_->next = af_middle_;
  af_middle_->next = af_tail_;
  af_middle_->prev = af_head_;
  af_tail_->prev = af_middle_;
}

void SweepContext::RemoveNode(Node* node)
{
  delete node;
}

void SweepContext::MapTriangleToNodes(TINTriangle& t)
{
  for (int i = 0; i < 3; i++) {
    if (!t.GetNeighbor(i)) {
      Node* n = front_->LocatePoint(t.PointCW(*t.GetPoint(i)));
      if (n)
        n->triangle = &t;
    }
  }
}

void SweepContext::RemoveFromMap(TINTriangle* triangle)
{
  map_.remove(triangle);
}

void SweepContext::MeshClean(TINTriangle& triangle)
{
  std::vector<TINTriangle *> triangles;
  triangles.push_back(&triangle);

  while(!triangles.empty()){
	TINTriangle *t = triangles.back();
	triangles.pop_back();

    if (t != NULL && !t->IsInterior()) {
      t->IsInterior(true);
      triangles_.push_back(t);
      for (int i = 0; i < 3; i++) {
        if (!t->constrained_edge[i])
          triangles.push_back(t->GetNeighbor(i));
      }
    }
  }
}

SweepContext::~SweepContext()
{

    // Clean up memory

    delete head_;
    delete tail_;
    delete front_;
    delete af_head_;
    delete af_middle_;
    delete af_tail_;

    typedef std::list<TINTriangle*> type_list;

    for(type_list::iterator iter = map_.begin(); iter != map_.end(); ++iter) {
        TINTriangle* ptr = *iter;
        delete ptr;
    }

     for(unsigned int i = 0; i < edge_list.size(); i++) {
        delete edge_list[i];
    }

}

}
