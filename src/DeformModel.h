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
#include <stdlib.h>

#include "vec3f.h"

#include <vector>
#include <set>
#include "feature.h"

using namespace std;
typedef vector<unsigned int> id_list;

class DeformBVHTree;

#include "box.h"

class DeformModel {
	unsigned int _num_vtx;
	unsigned int _num_frame;

	vec3f *_vtxs;
	vec3f *_cur_vtxs;
	vec3f *_prev_vtxs;

	unsigned int*_cur_flags;
	unsigned int*_prev_flags;
	vec3f *_cur_nrms;
	vec3f *_prev_nrms;

	id_list *_vtx_fids;
	BOX *_vtx_boxes;

	unsigned int _num_tri;
	tri3f *_tris;

	tri3e *_tri_edges;
	BOX *_fac_boxes;

	vec3f *_tri_nrms;
	vec3f *_old_tri_nrms;

	char *_tri_flags;

	unsigned int _num_edge;
	edge2f *_edges;
	BOX *_edg_boxes;

	// for building BVH
	DeformBVHTree *_tree;

	vec3f *_tri_centers;
	BOX *_tri_boxes;

	// for collide
	unsigned int _num_box_tests;
	unsigned int _num_tri_tests;
	unsigned int _num_ccd_tests;
	unsigned int _num_cov_tests;
	unsigned int _num_lp_tests;
	unsigned int _num_ccd_true;

	unsigned int _num_vf_test;
	unsigned int _num_ee_test;
	unsigned int _num_vf_true;
	unsigned int _num_ee_true;

	unsigned int _num_parts;
	unsigned int *_parts;

	void do_pairs();
	void do_non_adj_pairs();

	void BufferAdjacent();
	char get_status_1(unsigned int id1, unsigned int id2, unsigned int st1, unsigned int st2);
	char get_status_2(unsigned int id1, unsigned int id2, unsigned int st1, unsigned int st2);

	void Build(vec3f_list &vtxs, tri_list &tris);

public:
	DeformModel(vec3f_list &vtxs, tri_list &tris);
	~DeformModel();

	unsigned int GetFrames() { return _num_frame; }
	void Display();
	bool Deform(float t, float circle);

	void Init();
	void Clean();

	void UpdateVert(unsigned int prev, unsigned int next, float t);
	void UpdateVert(vec3f_list &vtxs);
	void UpdateBoxes();

	void BuildBVH(bool ccd);
	void RebuildBVH(bool ccd);
	float RefitBVH(bool ccd);

	void ResetCounter();
	void SelfCollide(bool ccd);

	FORCEINLINE int NumTri() { return _num_tri; }
	FORCEINLINE int NumBoxTest() { return _num_box_tests; }
	FORCEINLINE int NumTriTest() { return _num_tri_tests; }
	FORCEINLINE int NumContact() { return 0; }
	FORCEINLINE int NumCCDTest() { return _num_ccd_tests; }
	FORCEINLINE int NumCovTest() { return _num_cov_tests; }
	FORCEINLINE int NumLpTest() { return _num_lp_tests; }
	FORCEINLINE int NumCCDTrue() { return _num_ccd_true; }

	FORCEINLINE int NumVFTest() { return _num_vf_test; }
	FORCEINLINE int NumEETest() { return _num_ee_test; }
	FORCEINLINE int NumVFTrue() { return _num_vf_true; }
	FORCEINLINE int NumEETrue() { return _num_ee_true; }

	FORCEINLINE bool Covertex_E(unsigned int e1, unsigned int e2) {
		for (int i=0; i<2; i++)
			for (int j=0; j<2; j++)
				if (_edges[e1].vid(i) == _edges[e2].vid(j))
					return true;

		return false;
	}

	FORCEINLINE bool Coedge_F(unsigned int id1, unsigned int id2) {
		unsigned int num = 0;
		for (int i=0; i<3; i++)
			for (int j=0; j<3; j++) {
				if (_tris[id1].id(i) == _tris[id2].id(j)) {
					num++;
				}
			}

		return num > 1;
	}

	FORCEINLINE bool Covertex_F(unsigned int id1, unsigned int id2) {
		for (int i=0; i<3; i++)
			for (int j=0; j<3; j++) {
				if (_tris[id1].id(i) == _tris[id2].id(j))
					return true;
			}

		return false;
	}

	unsigned int Covertex_F(unsigned int id1, unsigned int id2, unsigned int &st1, unsigned int &st2);

	friend class DeformBVHTree;
	friend class DeformBVHNode;

	float intersect_vf(unsigned int fid1, unsigned int vid2, unsigned int fid2);
	bool check_vf(unsigned int fid, unsigned int vid);

	float intersect_ee(unsigned int e1, unsigned int e2, unsigned int fid1, unsigned int fid2);
	bool check_ee(unsigned int e1, unsigned int e2);

	float intersect_vf(unsigned int fid, unsigned int vid);
	float intersect_ee(unsigned int e1, unsigned int e2);
	float do_vf(unsigned int fid, unsigned int vid);
	float do_ee(unsigned int e1, unsigned int e2);

	void test_feature_0(unsigned id1, unsigned int id2);

	/// for orphan set
	void do_orphans();
	void get_orphans();
	void get_feature_1(unsigned int, unsigned int, unsigned int, unsigned int);
	void get_feature_2(unsigned int, unsigned int, unsigned int, unsigned int);
	void insert_ee(unsigned int, unsigned int);
	void insert_vf(unsigned int, unsigned int);
	bool test_orphan_ee(unsigned int e1, unsigned int e2);
	bool test_orphan_vf(unsigned int f, unsigned int v);
};