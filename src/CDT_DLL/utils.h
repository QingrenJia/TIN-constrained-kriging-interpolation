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
 
#ifndef UTILS_H
#define UTILS_H

// Otherwise #defines like M_PI are undeclared under Visual Studio
#define _USE_MATH_DEFINES

#include <exception>
#include <math.h>

namespace p2t {

const double PI_3div4 = 3 * M_PI / 4;
const double PI_div2 = 1.57079632679489661923;
const double EPSILON = 1e-12;

enum Orientation { CW, CCW, COLLINEAR };

/**
 * Forumla to calculate signed area<br>
 * Positive if CCW<br>
 * Negative if CW<br>
 * 0 if collinear<br>
 * <pre>
 * A[P1,P2,P3]  =  (x1*y2 - y1*x2) + (x2*y3 - y2*x3) + (x3*y1 - y3*x1)
 *              =  (x1-x3)*(y2-y3) - (y1-y3)*(x2-x3)
 * </pre>
 */
Orientation Orient2d(TINPoint& pa, TINPoint& pb, TINPoint& pc)
{
  double detleft = (pa.xyz[0] - pc.xyz[0]) * (pb.xyz[1] - pc.xyz[1]);
  double detright = (pa.xyz[1] - pc.xyz[1]) * (pb.xyz[0] - pc.xyz[0]);
  double val = detleft - detright;
  if (val > -EPSILON && val < EPSILON) {
    return COLLINEAR;
  } else if (val > 0) {
    return CCW;
  }
  return CW;
}

/*
bool InScanArea(TINPoint& pa, TINPoint& pb, TINPoint& pc, TINPoint& pd)
{
  double pdx = pd.xyz[0];
  double pdy = pd.xyz[1];
  double adx = pa.xyz[0] - pdx;
  double ady = pa.xyz[1] - pdy;
  double bdx = pb.xyz[0] - pdx;
  double bdy = pb.xyz[1] - pdy;

  double adxbdy = adx * bdy;
  double bdxady = bdx * ady;
  double oabd = adxbdy - bdxady;

  if (oabd <= EPSILON) {
    return false;
  }

  double cdx = pc.xyz[0] - pdx;
  double cdy = pc.xyz[1] - pdy;

  double cdxady = cdx * ady;
  double adxcdy = adx * cdy;
  double ocad = cdxady - adxcdy;

  if (ocad <= EPSILON) {
    return false;
  }

  return true;
}

*/

bool InScanArea(TINPoint& pa, TINPoint& pb, TINPoint& pc, TINPoint& pd)
{
  double oadb = (pa.xyz[0] - pb.xyz[0])*(pd.xyz[1] - pb.xyz[1]) - (pd.xyz[0] - pb.xyz[0])*(pa.xyz[1] - pb.xyz[1]);
  if (oadb >= -EPSILON) {
    return false;
  }

  double oadc = (pa.xyz[0] - pc.xyz[0])*(pd.xyz[1] - pc.xyz[1]) - (pd.xyz[0] - pc.xyz[0])*(pa.xyz[1] - pc.xyz[1]);
  if (oadc <= EPSILON) {
    return false;
  }
  return true;
}

}

#endif

