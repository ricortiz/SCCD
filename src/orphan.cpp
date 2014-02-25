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


#include "timing.h"
#include <algorithm>

#include "tri_pair.h"
extern vector<adjacent_pair> adj_2_list;
extern vector<adjacent_pair> adj_1_list;

#include "feature_pair.h"
static ee_list ee_keeper;
static vf_list vf_keeper;

#ifndef swapI
#define swapI(a, b) {\
	unsigned int tmp = a;\
	a = b;\
	b = tmp;\
}
#endif

void
DeformModel::do_orphans()
{
	for (ee_list::iterator it=ee_keeper.begin(); it!=ee_keeper.end(); it++) {
		unsigned int e1, e2;
		(*it).get_param(e1, e2);
		intersect_ee(e1, e2);
	}
	for (vf_list::iterator it=vf_keeper.begin(); it!=vf_keeper.end(); it++) {
		unsigned int v, f;
		(*it).get_param(f, v);
		intersect_vf(f, v);
	}
}

void
DeformModel::get_orphans()
{
	unsigned int id1, id2, st1, st2;
	char status;

	for (vector<adjacent_pair>::iterator it=adj_1_list.begin(); it!=adj_1_list.end(); it++) {
		(*it).get_param(id1, id2, st1, st2, status);
		get_feature_1(id1, id2, st1, st2);
	}
	sort(vf_keeper.begin(), vf_keeper.end());
	vf_keeper.erase(unique(vf_keeper.begin(), vf_keeper.end()), vf_keeper.end());

	sort(ee_keeper.begin(), ee_keeper.end());
	ee_keeper.erase(unique(ee_keeper.begin(), ee_keeper.end()), ee_keeper.end());

	for (vector<adjacent_pair>::iterator it=adj_2_list.begin(); it!=adj_2_list.end(); it++) {
		(*it).get_param(id1, id2, st1, st2, status);
		get_feature_2(id1, id2, st1, st2);
	}

	sort(vf_keeper.begin(), vf_keeper.end());
	vf_keeper.erase(unique(vf_keeper.begin(), vf_keeper.end()), vf_keeper.end());

	sort(ee_keeper.begin(), ee_keeper.end());
	ee_keeper.erase(unique(ee_keeper.begin(), ee_keeper.end()), ee_keeper.end());
}

bool
DeformModel::test_orphan_ee(unsigned int e1, unsigned int e2)
{
	for (int i=0; i<2; i++) {
		unsigned int f1 = _edges[e1].fid(i);
		if (f1 == -1) continue;

		for (int j=0; j<2; j++) {
			unsigned int f2 = _edges[e2].fid(j);
			if (f2 == -1) continue;

			if (!Covertex_F(f1, f2))
				return false;
		}
	}

	return true;
}

bool
DeformModel::test_orphan_vf(unsigned int f, unsigned int v)
{
	for (id_list::iterator it=_vtx_fids[v].begin(); it!=_vtx_fids[v].end(); it++) {
		unsigned int f2 = *it;
		if (!Covertex_F(f2, f))
			return false;
	}

	return true;
}

void
DeformModel::insert_ee(unsigned int e1, unsigned int e2)
{
	if (test_orphan_ee(e1, e2)) {
		ee_keeper.push_back(ee_pair(e1, e2));
	}
}

void
DeformModel::insert_vf(unsigned int f, unsigned int v)
{
	if (test_orphan_vf(f, v)) {
		vf_keeper.push_back(vf_pair(f, v));
	}
}

void
DeformModel::get_feature_2(
		unsigned id1, unsigned int id2, unsigned int st1, unsigned int st2)
{
	unsigned int e0 = _tri_edges[id1].id(st1);
	unsigned int e00 = _tri_edges[id1].id((st1+2)%3);
	unsigned int e1 = _tri_edges[id2].id(st2);
	unsigned int e11 = _tri_edges[id2].id((st2+2)%3);
	if (Covertex_E(e0, e1)) {
		swapI(e1, e11);
	}

	unsigned int fid = id2;
	unsigned int vid = _tris[id1].id(st1);
	insert_vf(fid, vid);
	insert_ee(e0, e1);

	// 2 VF test
	fid = id1;
	vid = _tris[id2].id(st2);
	insert_vf(fid, vid);
	insert_ee(e00, e11);
}

void
DeformModel::get_feature_1(unsigned id1, unsigned int id2, unsigned int st1, unsigned int st2)
{
	{
		// 2 VF test
		unsigned int fid = id1;
		for (int i=0; i<3; i++) {
			if (i == st2) continue; // skip one

			unsigned int vid = _tris[id2].id(i);
			insert_vf(fid, vid);
		}

		// 2 EE test
		unsigned int e1 = _tri_edges[id2].id((st2+1)%3);
		for (int i=0; i<2; i++) {
			unsigned int idx = (i == 0) ? 0 : 2;
			unsigned int e0 = _tri_edges[id1].id((st1+idx)%3);

			insert_ee(e0, e1);
		}
	}

	{
		// 2 VF test
		unsigned int fid = id2;
		for (int i=0; i<3; i++) {
			if (i == st1) continue; // skip one

			unsigned int vid = _tris[id1].id(i);
			insert_vf(fid, vid);
		}

		// 3 EE test
		unsigned int e0 = _tri_edges[id1].id((st1+1)%3);
		for (int i=0; i<3; i++) {
			unsigned int e1 = _tri_edges[id2].id(i);

			insert_ee(e0, e1);
		}
	}
}
