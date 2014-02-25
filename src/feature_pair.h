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
#include <assert.h>
#include "forceline.h"

class ee_pair {
	unsigned int _id[2];

public:
	FORCEINLINE ee_pair(unsigned int id1, unsigned int id2) {
		if (id1 > id2) {
			_id[0] = id1;
			_id[1] = id2;
		} else {
			_id[0] = id2;
			_id[1] = id1;
		}
	}

	FORCEINLINE void get_param(unsigned int &id1, unsigned int &id2) {
		id1 = _id[0];
		id2 = _id[1];
	}

	bool operator == (const ee_pair &other) const {
		return	(_id[0] == other._id[0] && _id[1] == other._id[1]);
	}

	bool operator < (const ee_pair &other) const {
		if (_id[0] == other._id[0])
			return _id[1] < other._id[1];
		else
			return _id[0] < other._id[0];
	}
};

class vf_pair {
	unsigned int _id[2];

public:
	FORCEINLINE vf_pair(unsigned int f, unsigned int v) {
		_id[0] = f;
		_id[1] = v;
	}

	FORCEINLINE void get_param(unsigned int &f, unsigned int &v) {
		f = _id[0];
		v = _id[1];
	}
	bool operator == (const vf_pair &other) const {
		return	(_id[0] == other._id[0] && _id[1] == other._id[1]);
	}

	bool operator < (const vf_pair &other) const {
		if (_id[0] == other._id[0])
			return _id[1] < other._id[1];
		else
			return _id[0] < other._id[0];
	}
};

#include <vector>
typedef std::vector<ee_pair> ee_list;
typedef std::vector<vf_pair> vf_list;
