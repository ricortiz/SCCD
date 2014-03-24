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

#include <math.h>
#include <stdio.h>
#include <assert.h>

#include "timing.h"
extern CBVHTimer tm;

#include <vector>
#include <algorithm>
#include <list>

using namespace std;

#include "tri_pair.h"
vector<adjacent_pair> adj_1_list;
vector<adjacent_pair> adj_2_list;
non_adjacent_pair_list non_adj_list;

#include "DeformModel.h"
#include "DeformBVH.h"


#include <omp.h>

DeformModel::~DeformModel()
{
	Clean();
}

void
DeformModel::Init()
{
	_num_vtx = 0;
	_num_frame = 0;

	_vtxs = NULL;
	_cur_vtxs = NULL;
	_prev_vtxs = NULL;

	_cur_nrms = NULL;
	_prev_nrms = NULL;
	_cur_flags = _prev_flags = NULL;

	_vtx_fids = NULL;

	_num_tri = 0;
	_tris = NULL;
	_tri_edges = NULL;
	_tri_nrms = NULL;
	_old_tri_nrms = NULL;
	_tri_flags = NULL;

	_num_edge = 0;
	_edges = NULL;

	_tree = NULL;
	_tri_centers = NULL;
	_tri_boxes = NULL;

	_num_box_tests = 0;
	_num_tri_tests = 0;
	_num_cov_tests = 0;
	_num_lp_tests = 0;
	_num_ccd_true = 0;
	_num_ccd_tests = 0;

	_num_vf_test = 0;
	_num_ee_test = 0;
	_num_vf_true = 0;
	_num_ee_true = 0;

	_num_parts = 0;
	_parts = NULL;
}

void
DeformModel::Clean()
{
	if (_vtxs) delete [] _vtxs;
	if (_vtx_fids) delete [] _vtx_fids;
	if (_cur_vtxs) delete [] _cur_vtxs;
	if (_prev_vtxs) delete [] _prev_vtxs;

	if (_cur_nrms) delete [] _cur_nrms;
	if (_prev_nrms) delete [] _prev_nrms;
	if (_cur_flags) delete [] _cur_flags;
	if (_prev_flags) delete [] _prev_flags;

	if (_vtx_boxes) delete [] _vtx_boxes;

	if (_edges) delete [] _edges;
	if (_edg_boxes) delete [] _edg_boxes;

	if (_tris) delete [] _tris;
	if (_tri_edges) delete [] _tri_edges;
	if (_fac_boxes) delete [] _fac_boxes;
	
	if (_tri_nrms) delete [] _tri_nrms;
	if (_old_tri_nrms) delete [] _old_tri_nrms;
	if (_tri_flags) delete [] _tri_flags;

	if (_parts) delete [] _parts;
}

void
DeformModel::Display()
{
}

bool
DeformModel::Deform(float t, float circle)
{
	if (t >= circle)
		return false;

	t = t/circle*_num_frame;
	unsigned int prev_frame = ((unsigned int)t)%_num_frame;
	unsigned int next_frame = (prev_frame+1)%_num_frame;

	if (next_frame < prev_frame)
		return false;

	tm.startTiming(9);
	tm.startTiming(11);
	UpdateVert(prev_frame, next_frame, t-prev_frame);
	UpdateBoxes();

	tm.endTiming(11);
	tm.endTiming(9);

	return true;
}

#define swap(a, b) {\
	vec3f *tmp = a;\
	a = b;\
	b = tmp;\
}

inline vec3f interp(const vec3f &p1, const vec3f &p2, float t)
{
	return p1*(1-t)+p2*t;
}


inline vec3f update(vec3f &v1, vec3f &v2, vec3f &v3)
{
	vec3f s = (v2-v1);
	return s.cross(v3-v1);
}

inline vec3f
update(tri3f &tri, vec3f *vtxs)
{
	vec3f &v1 = vtxs[tri.id0()];
	vec3f &v2 = vtxs[tri.id1()];
	vec3f &v3 = vtxs[tri.id2()];

	return update(v1, v2, v3);
}

void
DeformModel::UpdateVert(vec3f_list &vtxs)
{
	swap(_prev_vtxs, _cur_vtxs);

	for (int i=0; i<_num_vtx; i++) {
		_cur_vtxs[i] = vtxs[i];
	}
}

void
DeformModel::UpdateVert(unsigned int prev, unsigned int next, float t)
{
	vec3f *prev_pt = _vtxs+prev*_num_vtx;
	vec3f *next_pt = _vtxs+next*_num_vtx;

	swap(_prev_vtxs, _cur_vtxs);

	for (int i=0; i<_num_vtx; i++) {
		_cur_vtxs[i] = interp(prev_pt[i], next_pt[i], t);
	}
}

void
DeformModel::UpdateBoxes()
{
	for (int i=0; i<_num_vtx; i++) {
		_vtx_boxes[i] = BOX(_cur_vtxs[i]) + _prev_vtxs[i];
	}

	for (int i=0; i<_num_edge; i++) {
		unsigned int id0 = _edges[i].vid(0);
		unsigned int id1 = _edges[i].vid(1);

		_edg_boxes[i] = _vtx_boxes[id0] + _vtx_boxes[id1];
	}

	for (int i=0; i<_num_tri; i++) {
		unsigned int id0 = _tris[i].id0();
		unsigned int id1 = _tri_edges[i].id(1);

		_fac_boxes[i] = _vtx_boxes[id0] + _edg_boxes[id1];
	}
}

void
DeformModel::BuildBVH(bool ccd)
{
	_tree = new DeformBVHTree(this, ccd);
	_tree->refit();
}

void
DeformModel::RebuildBVH(bool ccd)
{
	delete _tree;
	_tree = new DeformBVHTree(this, ccd);
	_tree->refit();
}

float DeformModel::RefitBVH(bool ccd)
{
	return _tree->refit();
}

void DeformModel::ResetCounter()
{
	// reset results
	_num_box_tests = 0;
	_num_tri_tests = 0;
	_num_ccd_tests = 0;
	_num_cov_tests = 0;
	_num_lp_tests = 0;
	_num_ccd_true = 0;

	_num_vf_true = 0;
	_num_ee_true = 0;
	_num_vf_test = 0;
	_num_ee_test = 0;
}

void DeformModel::SelfCollide(bool ccd)
{
	if (_tree == NULL)
		return;

	non_adj_list.clear();

	tm.startTiming(1);
	_tree->self_collide();
	tm.endTiming(1);

	do_pairs();
}

void
DeformModel::do_pairs()
{
	tm.startTiming(4);
	do_non_adj_pairs();
	tm.endTiming(4);

	tm.startTiming(2);
	do_orphans();
	tm.endTiming(2);
}

void
DeformModel::do_non_adj_pairs()
{
	for (vector<non_adjacent_pair>::iterator it=non_adj_list.begin(); it != non_adj_list.end(); it++)
	{
		unsigned int id1, id2;
		(*it).get_param(id1, id2);
		test_feature_0(id1, id2);
	}
}

unsigned int
DeformModel::Covertex_F(unsigned int id1, unsigned int id2, unsigned int &st1, unsigned int &st2)
{
	unsigned int keeps[4];
	unsigned int num = 0;
	for (int i=0; i<3; i++)
		for (int j=0; j<3; j++) {
			if (_tris[id1].id(i) == _tris[id2].id(j)) {
				if (num < 2) {
					keeps[num*2] = i;
					keeps[num*2+1] = j;
				}
				num++;
			}
		}

	assert(num <= 3);
	if (num == 1) {
		st1 = keeps[0];
		st2 = keeps[1];
	} else
	if (num == 2) {
		for (int i=0; i<3; i++) {
			if (i != keeps[0] && i != keeps[2]) {
				st1 = i;
			}
			if (i != keeps[1] && i != keeps[3]) {
				st2 = i;
			}
		}
	} else {
		st1 = st2 = 0;
	}

	return num;
}

//##########################################################
#include <stdio.h>
#include <iostream>

DeformModel::DeformModel(vec3f_list &vtxs, tri_list &tris)
{
	Build(vtxs, tris);
}

void
DeformModel::Build(vec3f_list &vtxs, tri_list &tris)
{
	Init();

	_num_vtx = vtxs.size();
	_cur_vtxs = new vec3f[_num_vtx];
	_cur_vtxs = new vec3f[_num_vtx];
	_prev_vtxs = new vec3f[_num_vtx];
	_vtx_boxes = new BOX[_num_vtx];
	_vtx_fids = new id_list[_num_vtx];
	for (int i=0; i<_num_vtx; i++) {
		_cur_vtxs[i] = _prev_vtxs[i] = vtxs[i];
	}
	cout << "Vtx # = " << _num_vtx << endl;

	list<edge2f> edgelist_temp;
	vector<tri3f> trilist_temp;

	for (int i=0; i<tris.size(); i++) {
		unsigned int id0 = tris[i].id0(), id1 = tris[i].id1(), id2 = tris[i].id2(), fid = i;

		trilist_temp.push_back(tris[i]);
		edgelist_temp.push_back(edge2f(id0, id1, fid));
		edgelist_temp.push_back(edge2f(id1, id2, fid));
		edgelist_temp.push_back(edge2f(id2, id0, fid));
	}
	edgelist_temp.sort();


	list<edge2f> edge_unqie;
	for (list<edge2f>::iterator it=edgelist_temp.begin(); it!=edgelist_temp.end(); it++) {
		if (!edge_unqie.empty() && *it == edge_unqie.back()) { // find duplicated with other fid
			unsigned int fid = (*it).fid(0);
			assert(fid != -1);
			edge_unqie.back().set_fid2(fid);
		} else
			edge_unqie.push_back(*it);
	}

	edgelist_temp.clear();
	vector<edge2f> edge_array;

	_num_edge = (unsigned int)edge_unqie.size();
	_edges = new edge2f[_num_edge];
	_edg_boxes = new BOX[_num_edge];

	unsigned int t=0;
	for (list<edge2f>::iterator it=edge_unqie.begin(); it != edge_unqie.end(); it++)
	{
		_edges[t++] = *it;
		edge_array.push_back(*it);
	}

	// copy over temp list to static array
	_num_tri = (unsigned int)trilist_temp.size();
	_tris = new tri3f[_num_tri];
	_tri_nrms = new vec3f[_num_tri];
	_old_tri_nrms = NULL;

	_cur_nrms = new vec3f[_num_tri];
	_prev_nrms = new vec3f[_num_tri];
	_cur_flags = new unsigned int[_num_tri];
	_prev_flags = new unsigned int [_num_tri];
	for (int i=0; i<_num_tri; i++)
		_cur_flags[i] = _prev_flags[i] = -1;

	_tri_edges = new tri3e[_num_tri];
	_fac_boxes = new BOX[_num_tri];
	_tri_flags = new char[_num_tri];

	vector <edge2f>::iterator first = edge_array.begin();
	vector <edge2f>::iterator last = edge_array.end();
	for (t = 0; t < _num_tri; t++) {
		_tris[t] = trilist_temp[t];

		vector <edge2f>::iterator it1 = lower_bound(first, last, edge2f(_tris[t].id0(), _tris[t].id1(), 0));
		vector <edge2f>::iterator it2 = lower_bound(first, last, edge2f(_tris[t].id1(), _tris[t].id2(), 0));
		vector <edge2f>::iterator it3 = lower_bound(first, last, edge2f(_tris[t].id2(), _tris[t].id0(), 0));

		_tri_edges[t].set(it1-first, it2-first, it3-first);
	}

	cout << "Edge # = " << _num_edge << endl;
	cout << "Tri # = " << _num_tri << endl;
	UpdateBoxes();

	// build _vtx_fids
	for (unsigned t = 0; t < _num_tri; t++)
		for (int i=0; i<3; i++) {
			unsigned int vid = _tris[t].id(i);
			_vtx_fids[vid].push_back(t);
		}

	BufferAdjacent();
}


#include "tri_pair.h"
extern vector<adjacent_pair> adj_2_list;
extern vector<adjacent_pair> adj_1_list;
#define swapI(a, b) {\
	unsigned int tmp = a;\
	a = b;\
	b = tmp;\
}

char
DeformModel::get_status_1(unsigned int id1, unsigned int id2, unsigned int st1, unsigned int st2)
{
	unsigned int e0 = _tri_edges[id1].id(st1);
	unsigned int e00 = _tri_edges[id1].id((st1+2)%3);

	unsigned int fid1 = _edges[e0].fid(0);
	unsigned int f1_mates = (fid1 == id1) ? _edges[e0].fid(1) : fid1;

	unsigned int e1 = _tri_edges[id2].id(st2);
	unsigned int e11 = _tri_edges[id2].id((st2+2)%3);
	if (Covertex_E(e0, e1)) {
		swapI(e1, e11);
	}
	unsigned int fid2 = _edges[e11].fid(0);
	unsigned int f2_mates = (fid2 == id2) ? _edges[e11].fid(1) : fid2;

	if (f1_mates != -1)
		if (f2_mates != -1)
			return 3; // just skip it
		else
			return 2; //group 2 only
	else
		if (f2_mates != -1)
			return 1; //group 1 only
		else
			return 0; //both group
}

char
DeformModel::get_status_2(unsigned int id1, unsigned int id2, unsigned int st1, unsigned int st2)
{
	unsigned int edg1 = _tri_edges[id1].id((st1+1)%3);
	unsigned int fid1 = _edges[edg1].fid(0);
	unsigned int f1_mates = (fid1 == id1) ? _edges[edg1].fid(1) : fid1;

	unsigned int edg2 = _tri_edges[id2].id((st2+1)%3);
	unsigned int fid2 = _edges[edg2].fid(0);
	unsigned int f2_mates = (fid2 == id2) ? _edges[edg2].fid(1) : fid2;


	if (f1_mates != -1)
		if (f2_mates != -1)
			return 3; //just skip it ...
		else
			return 2; //group 2 only
	else
		if (f2_mates != -1)
			return 1; //group 1 only
		else
			return 0; //both group
}

void
DeformModel::BufferAdjacent()
{
	adj_2_list.clear();
	for (unsigned int i=0; i<_num_edge; i++) {
		unsigned int id1 = _edges[i].fid(0);
		unsigned int id2 = _edges[i].fid(1);

		if (id1 == -1 || id2 == -1) continue;
		if (id1 < id2) swapI(id1, id2);

		unsigned int st1 = 0, st2 = 0;
		unsigned int cov = Covertex_F(id1, id2, st1, st2);
		assert(cov == 2);

		char status = get_status_2(id1, id2, st1, st2);
		//if (status == 3) continue;

		adj_2_list.push_back(adjacent_pair(id1, id2, st1, st2, status));
	}

	adj_1_list.clear();
	for (unsigned int i=0; i<_num_vtx; i++) {
		for (id_list::iterator it1=_vtx_fids[i].begin(); it1!=_vtx_fids[i].end(); it1++) {
			for (id_list::iterator it2=it1; it2!=_vtx_fids[i].end(); it2++) {
				if (it1 == it2) continue;

				unsigned int id1 = *it1;
				unsigned int id2 = *it2;
				if (id1 < id2) swapI(id1, id2);

				unsigned int st1 = 0, st2 = 0;
				unsigned int cov = Covertex_F(id1, id2, st1, st2);
				if (cov == 2) continue;

				char status = get_status_1(id1, id2, st1, st2);
				//if (status == 3) continue;

				adj_1_list.push_back(adjacent_pair(id1, id2, st1, st2, status));
			}
		}
	}

	sort(adj_1_list.begin(), adj_1_list.end());
	adj_1_list.erase(unique(adj_1_list.begin(), adj_1_list.end()), adj_1_list.end());

	get_orphans();
}