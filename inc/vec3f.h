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

#include <limits>
#include <cmath>
#include "forceline.h"

template<typename T>
class vec3 {
	union {
		struct {
		T x, y, z;
		};
		struct {
		T v[3];
		};
	};
public:

	FORCEINLINE vec3 ()
	{x=0; y=0; z=0;}

	FORCEINLINE vec3(const vec3 &v)
	{
		x = v.x;
		y = v.y;
		z = v.z;
	}

	FORCEINLINE vec3(const T *v)
	{
		x = v[0];
		y = v[1];
		z = v[2];
	}

	FORCEINLINE vec3(T x, T y, T z)
	{
		this->x = x;
		this->y = y;
		this->z = z;
	}

	FORCEINLINE T operator [] ( int i ) const {return v[i];}

	FORCEINLINE vec3 &operator += (const vec3 &v) {
		x += v.x;
		y += v.y;
		z += v.z;
		return *this;
	}

	FORCEINLINE vec3 &operator -= (const vec3 &v) {
		x -= v.x;
		y -= v.y;
		z -= v.z;
		return *this;
	}

	FORCEINLINE void negate() {
		x = -x;
		y = -y;
		z = -z;
	}

	FORCEINLINE vec3 operator - () const {
		return vec3(-x, -y, -z);
	}

	FORCEINLINE vec3 operator+ (const vec3 &v) const
	{
		return vec3(x+v.x, y+v.y, z+v.z);
	}

	FORCEINLINE vec3 operator- (const vec3 &v) const
	{
		return vec3(x-v.x, y-v.y, z-v.z);
	}

	FORCEINLINE vec3 operator *(T t) const
	{
		return vec3(x*t, y*t, z*t);
	}

     // cross product
     FORCEINLINE const vec3 cross(const vec3 &vec) const
     {
          return vec3(y*vec.z - z*vec.y, z*vec.x - x*vec.z, x*vec.y - y*vec.x);
     }

	 FORCEINLINE T dot(const vec3 &vec) const {
		 return x*vec.x+y*vec.y+z*vec.z;
	 }

	 FORCEINLINE void normalize() 
	 { 
		 T sum = x*x+y*y+z*z;
		 if (sum > std::numeric_limits<T>::epsilon()*std::numeric_limits<T>::epsilon()) {
			 T base = T(1.0/sqrt(sum));
			 x *= base;
			 y *= base;
			 z *= base;
		 }
	 }

	 FORCEINLINE T length() {
		 return T(sqrt(x*x + y*y + z*z));
	 }

	FORCEINLINE vec3 & set_value( const T &vx, const T &vy, const T &vz)
	{ x = vx; y = vy; z = vz; return *this; }

	FORCEINLINE bool equal_abs(const vec3 &other) {
		return x == other.x && y == other.y && z == other.z;
	}

	FORCEINLINE T square_norm() const {
		return x*x+y*y+z*z;
	}
};

typedef vec3<float> vec3f;

#include <vector>
typedef std::vector<vec3f> vec3f_list;
