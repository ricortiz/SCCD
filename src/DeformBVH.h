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

#include "box.h"

class DeformBVHNode;
class DeformModel;

class bvh_front_list;
class non_adjacent_pair_list;

class DeformBVHNode {
	BOX _box;

	unsigned int _id;

	DeformBVHNode *_parent;
	DeformBVHNode *_left;
	DeformBVHNode *_right;

public:
	DeformBVHNode();
	DeformBVHNode(DeformBVHNode *, unsigned int);
	DeformBVHNode(DeformBVHNode *, unsigned int *, unsigned int, DeformModel *);

	~DeformBVHNode();

	void getChildren(DeformBVHNode *&n1, DeformBVHNode *&n2, DeformBVHNode *&n3, DeformBVHNode *&n4);
	void mergeBox(DeformBVHNode *n1, DeformBVHNode *n2, DeformBVHNode *n3, DeformBVHNode *n4);

	void collide(DeformBVHNode *);
	void self_collide();

	void refit();
	bool find(unsigned int);

	FORCEINLINE DeformBVHNode *getLeftChild() { return _left; }
	FORCEINLINE DeformBVHNode *getRightChild() { return _right; }
	FORCEINLINE DeformBVHNode *getParent() { return _parent; }

	FORCEINLINE int getTriID() { return _id; }
	FORCEINLINE bool isLeaf() { return _left == NULL; }
	FORCEINLINE bool isRoot() { return _parent == NULL;}

friend class DeformBVHTree;
};

class DeformBVHTree {
	DeformModel		*_mdl;
	DeformBVHNode	*_root;
	unsigned int	_part;
	unsigned int *idx_buffer;
	unsigned int *_tri_idxes;

public:
	DeformBVHTree(DeformModel *, bool, unsigned int = -1);
	~DeformBVHTree();

	void Construct(DeformModel *, bool);

	float refit(bool = true);
	float refit1(bool = true);

	void collide(DeformBVHTree *);
	void self_collide();
	void self_collide_norecursive();
	void parallel_self_collide();
	void check_connection();

	void do_task_1();
	void do_task_2();

	BOX box();
	FORCEINLINE DeformBVHNode *getRoot() { return _root; }

friend class DeformBVHNode;
friend class DeformModel;
};
