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

#include "DeformModel.h"
#include "aabb.h"

#include "timing.h"
extern CBVHTimer tm;

#include "ccdAPI.h"
extern ccdEETestCallback *cbFuncEE;
extern ccdVFTestCallback *cbFuncVF;

#ifndef swapI
#define swapI(a, b) {\
	unsigned int tmp = a;\
	a = b;\
	b = tmp;\
}
#endif

extern float
Intersect_VF(const vec3f &ta0, const vec3f &tb0, const vec3f &tc0,
			 const vec3f &ta1, const vec3f &tb1, const vec3f &tc1,
			 const vec3f &q0, const vec3f &q1,
			 vec3f &qi, vec3f &baryc);
extern float
Intersect_EE(const vec3f &ta0, const vec3f &tb0, const vec3f &tc0, const vec3f &td0,
			 const vec3f &ta1, const vec3f &tb1, const vec3f &tc1, const vec3f &td1,
			 vec3f &qi);

inline vec3f norm(vec3f &a, vec3f &b, vec3f &c, vec3f &d)
{
	return (a-b).cross(c-d);
}

inline vec3f norm(vec3f &p1, vec3f &p2, vec3f &p3)
{
	return (p2-p1).cross(p3-p1);
}

inline bool side(vec3f &a, vec3f &b, vec3f &c, vec3f &p)
{
	return norm(a, b, c).dot(p-a) > 0;
}

extern unsigned int g_frame;

inline bool
check_abcd(vec3f &a0, vec3f &b0, vec3f &c0, vec3f &d0,
					vec3f &a1, vec3f &b1, vec3f &c1, vec3f &d1)
{
	vec3f n0 = norm(a0, b0, c0);
	vec3f n1 = norm(a1, b1, c1);
	vec3f delta = norm(a1-a0, b1-b0, c1-c0); //((b1-b0)-(a1-a0)).cross((c1-c0)-(a1-a0));
	vec3f nX = (n0+n1-delta)*0.5;

	vec3f pa0 = d0-a0;
	vec3f pa1 = d1-a1;

	float A = n0.dot(pa0);
	float B = n1.dot(pa1);
	float C = nX.dot(pa0);
	float D = nX.dot(pa1);
	float E = n1.dot(pa0);
	float F = n0.dot(pa1);

	if (A > 0 && B > 0 && (2*C+F) > 0 && (2*D+E) > 0)
		return false;

	if (A < 0 && B < 0 && (2*C+F) < 0 && (2*D+E) < 0)
		return false;

	return true;
}

bool
DeformModel::check_vf(unsigned int fid1, unsigned int vid2)
{
	unsigned v0 = _tris[fid1].id0();
	unsigned v1 = _tris[fid1].id1();
	unsigned v2 = _tris[fid1].id2();

	vec3f &a0 = _prev_vtxs[v0];
	vec3f &b0 = _prev_vtxs[v1];
	vec3f &c0 = _prev_vtxs[v2];
	vec3f &p0 = _prev_vtxs[vid2];

	vec3f &a1 = _cur_vtxs[v0];
	vec3f &b1 = _cur_vtxs[v1];
	vec3f &c1 = _cur_vtxs[v2];
	vec3f &p1 = _cur_vtxs[vid2];

	return check_abcd(a0, b0, c0, p0, a1, b1, c1, p1);
}

bool
DeformModel::check_ee(unsigned int e1, unsigned int e2)
{
	unsigned v0 = _edges[e1].vid(0);
	unsigned v1 = _edges[e1].vid(1);
	unsigned w0 = _edges[e2].vid(0);
	unsigned w1 = _edges[e2].vid(1);

	vec3f &a0 = _prev_vtxs[v0];
	vec3f &b0 = _prev_vtxs[v1];
	vec3f &c0 = _prev_vtxs[w0];
	vec3f &d0 = _prev_vtxs[w1];

	vec3f &a1 = _cur_vtxs[v0];
	vec3f &b1 = _cur_vtxs[v1];
	vec3f &c1 = _cur_vtxs[w0];
	vec3f &d1 = _cur_vtxs[w1];

	return check_abcd(a0, b0, c0, d0, a1, b1, c1, d1);
}

float
DeformModel::intersect_vf(unsigned int fid1, unsigned int vid2, unsigned int fid2)
{
	if (!_fac_boxes[fid1].overlaps(_vtx_boxes[vid2]))
		return -1.f;

	_num_ccd_tests++;

	for (id_list::iterator it1=_vtx_fids[vid2].begin(); it1!=_vtx_fids[vid2].end(); it1++) {
			unsigned int fid = *it1;
			if (!Covertex_F(fid1, fid)) {
				if (fid == fid2) {
					if (check_vf(fid1, vid2) == false)
						return -1.f;
					else
						return do_vf(fid1, vid2);
				} else
					return -1.f;
			}
	}
	return -1.f;
}

float
DeformModel::intersect_vf(unsigned int fid, unsigned int vid)
{
	if (!_fac_boxes[fid].overlaps(_vtx_boxes[vid]))
		return -1.f;

	_num_ccd_tests++;
	if (check_vf(fid, vid) == false)
		return -1.f;
	else
		return do_vf(fid, vid);
}

float
DeformModel::do_vf(unsigned int fid, unsigned int vid)
{
	_num_vf_test++;

	vec3f qi, baryc;
	unsigned v0 = _tris[fid].id0();
	unsigned v1 = _tris[fid].id1();
	unsigned v2 = _tris[fid].id2();

	//tm.startTiming(8);
	float ret = Intersect_VF(
		_prev_vtxs[v0],  _prev_vtxs[v1], _prev_vtxs[v2],
		_cur_vtxs[v0],   _cur_vtxs[v1],  _cur_vtxs[v2],
		_prev_vtxs[vid], _cur_vtxs[vid], qi, baryc);
	//tm.endTiming(8);

	if (ret> -0.5) {
		_num_lp_tests++;
		_num_vf_true++;
		if (cbFuncVF)
			(*cbFuncVF)(vid, fid, ret);
	}

	return ret;
}

float
DeformModel::intersect_ee(unsigned int e1, unsigned int e2, unsigned int f1, unsigned int f2)
{
	if (!_edg_boxes[e1].overlaps(_edg_boxes[e2]))
		return -1.f;

	unsigned int e[2];
	unsigned int f[2];

	if (e1 > e2) {
		e[0] = e1, e[1] = e2;
		f[0] = f1, f[1] = f2;
	} else {
		e[0] = e2, e[1] = e1;
		f[0] = f2, f[1] = f1;
	}

	for (int i=0; i<2; i++)
		for (int j=0; j<2; j++) {
			unsigned int ff1 = _edges[e[0]].fid(i);
			unsigned int ff2 = _edges[e[1]].fid(j);

			if (ff1 == -1 || ff2 == -1)
				continue;

			if (!Covertex_F(ff1, ff2)) {
				if (ff1 == f[0] && ff2 == f[1]) {
						if (check_ee(e1, e2) == false)
							return -1.f;
						else
							return do_ee(e1, e2);
				}else
					return -1.f;
			}
		}

	return -1.f;
}

float
DeformModel::intersect_ee(unsigned int e1, unsigned int e2)
{
	if (!_edg_boxes[e1].overlaps(_edg_boxes[e2]))
		return -1.f;

	_num_ccd_tests++;

	if (check_ee(e1, e2) == false)
		return -1.f;
	else
		return do_ee(e1, e2);
}

float
DeformModel::do_ee(unsigned int e1, unsigned int e2)
{
	_num_ee_test++;

	vec3f qi;
	unsigned v0 = _edges[e1].vid(0);
	unsigned v1 = _edges[e1].vid(1);
	unsigned w0 = _edges[e2].vid(0);
	unsigned w1 = _edges[e2].vid(1);

	//tm.startTiming(7);
	float ret = Intersect_EE(
		_prev_vtxs[v0], _prev_vtxs[v1],
		_prev_vtxs[w0], _prev_vtxs[w1],
		_cur_vtxs[v0], _cur_vtxs[v1],
		_cur_vtxs[w0], _cur_vtxs[w1],
		qi);
	//tm.endTiming(7);

	if (ret> -0.5) {
		_num_lp_tests++;
		_num_ee_true++;

		if (cbFuncEE) {
			(*cbFuncEE)(v0, v1, w0, w1, ret);
		}
	}

	return ret;
}

void
DeformModel::test_feature_0(unsigned id1, unsigned int id2)
{
#if defined(FOR_BART) || defined(FOR_DRAGON)
	// 6 VF test
	for (int i=0; i<3; i++) {
		intersect_vf(id1, _tris[id2].id(i));
		intersect_vf(id2, _tris[id1].id(i));
	}

	// 9 EE test
	for (int i=0; i<3; i++)
		for (int j=0; j<3; j++) {
		unsigned int e0 = _tri_edges[id1].id(i);
		unsigned int e1 = _tri_edges[id2].id(j);
		
		intersect_ee(e0, e1);
	}

	return;
#endif

	// 6 VF test
	for (int i=0; i<3; i++) {
		intersect_vf(id1, _tris[id2].id(i), id2);
		intersect_vf(id2, _tris[id1].id(i), id1);
	}

	// 9 EE test
	for (int i=0; i<3; i++)
		for (int j=0; j<3; j++) {
		unsigned int e0 = _tri_edges[id1].id(i);
		unsigned int e1 = _tri_edges[id2].id(j);
		
		intersect_ee(e0, e1, id1, id2);
	}
}
