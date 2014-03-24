/*************************************************************************\

  Copyright 2010 The University of North Carolina at Chapel Hill.
  All Rights Reserved.

  Permission to use, copy, modify and distribute this software and its
  documentation for educational, research and non-profit purposes, without
   fee, and without a written agreement is hereby granted, provided that the
  above copyright notice and the following three paragraphs appear in all
  copies.

  IN NO EVENT SHALL THE UNIVERSITY OF NORTH CAROLINA AT CHAPEL HILL BE
  LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR
  CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE
  USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY
  OF NORTH CAROLINA HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH
  DAMAGES.

  THE UNIVERSITY OF NORTH CAROLINA SPECIFICALLY DISCLAIM ANY
  WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE SOFTWARE
  PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
  NORTH CAROLINA HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT,
  UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

  The authors may be contacted via:

  US Mail:             GAMMA Research Group at UNC
                       Department of Computer Science
                       Sitterson Hall, CB #3175
                       University of N. Carolina
                       Chapel Hill, NC 27599-3175

  Phone:               (919)962-1749

  EMail:              geom@cs.unc.edu; tang_m@zju.edu.cn


\**************************************************************************/

#pragma once

#include "vec3f.h"


inline void
vmin(vec3f &a, const vec3f &b)
{
	a.set_value(
		std::min(a[0], b[0]),
		std::min(a[1], b[1]),
		std::min(a[2], b[2]));
}

inline void
vmax(vec3f &a, const vec3f &b)
{
	a.set_value(
		std::max(a[0], b[0]),
		std::max(a[1], b[1]),
		std::max(a[2], b[2]));
}


class aabb {
public:
	vec3f _min;
	vec3f _max;

	FORCEINLINE aabb() {
		empty();
	}

	FORCEINLINE aabb(const vec3f &v) {
		_min = _max = v;
	}

	FORCEINLINE aabb(const vec3f &a, const vec3f &b) {
		_min = a;
		_max = a;
		vmin(_min, b);
		vmax(_max, b);
	}

	FORCEINLINE bool overlaps(const aabb& b) const
	{
		if (_min[0] > b._max[0]) return false;
		if (_min[1] > b._max[1]) return false;
		if (_min[2] > b._max[2]) return false;

		if (_max[0] < b._min[0]) return false;
		if (_max[1] < b._min[1]) return false;
		if (_max[2] < b._min[2]) return false;

		return true;
	}

	FORCEINLINE bool overlaps(const aabb &b, aabb &ret) const
	{
		if (!overlaps(b))
			return false;

		ret._min.set_value(
			std::max(_min[0],  b._min[0]),
			std::max(_min[1],  b._min[1]),
			std::max(_min[2],  b._min[2]));

		ret._max.set_value(
			std::min(_max[0], b._max[0]),
			std::min(_max[1], b._max[1]),
			std::min(_max[2], b._max[2]));

		return true;
	}

	FORCEINLINE bool inside(const vec3f &p) const
	{
		if (p[0] < _min[0] || p[0] > _max[0]) return false;
		if (p[1] < _min[1] || p[1] > _max[1]) return false;
		if (p[2] < _min[2] || p[2] > _max[2]) return false;

		return true;
	}

	FORCEINLINE aabb &operator += (const vec3f &p)
	{
		vmin(_min, p);
		vmax(_max, p);
		return *this;
	}

	FORCEINLINE aabb &operator += (const aabb &b)
	{
		vmin(_min, b._min);
		vmax(_max, b._max);
		return *this;
	}

	FORCEINLINE aabb operator + ( const aabb &v) const
		{ aabb rt(*this); return rt += v; }

	FORCEINLINE float width()  const { return _max[0] - _min[0]; }
	FORCEINLINE float height() const { return _max[1] - _min[1]; }
	FORCEINLINE float depth()  const { return _max[2] - _min[2]; }
	FORCEINLINE vec3f center() const { return (_min+_max)*0.5; }
	FORCEINLINE float volume() const { return width()*height()*depth(); }

	FORCEINLINE void empty() {
		_max = vec3f(-std::numeric_limits< float >::max(), -std::numeric_limits< float >::max(), -std::numeric_limits< float >::max());
		_min = vec3f(std::numeric_limits< float >::max(), std::numeric_limits< float >::max(), std::numeric_limits< float >::max());
	}
};