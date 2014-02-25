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

#include <stdlib.h>
#include <assert.h>

#include "DeformBVH.h"
#include "DeformModel.h"
#include "timing.h"
extern CBVHTimer tm;

static DeformModel *s_mdl1, *s_mdl2;

#include "tri_pair.h"
extern non_adjacent_pair_list non_adj_list;

static DeformModel *s_mdl;

void
DeformBVHNode::getChildren(DeformBVHNode *&n1, DeformBVHNode *&n2, DeformBVHNode *&n3, DeformBVHNode *&n4)
{
	n1 = getLeftChild()->getLeftChild();
	n2 = getLeftChild()->getRightChild();
	n3 = getRightChild()->getLeftChild();
	n4 = getRightChild()->getRightChild();
}

void
DeformBVHNode::mergeBox(DeformBVHNode *n1, DeformBVHNode *n2, DeformBVHNode *n3, DeformBVHNode *n4)
{
	getLeftChild()->_box = n1->_box + n2->_box;
	getRightChild()->_box = n3->_box + n4->_box;
	_box = getLeftChild()->_box + getRightChild()->_box;
}

float
DeformBVHTree::refit(bool openmp)
{
	s_mdl = _mdl;

	getRoot()->refit();

	return 0.f;
}

void
DeformBVHTree::collide(DeformBVHTree *other)
{
	s_mdl1 = _mdl;
	s_mdl2 = other->_mdl;

	getRoot()->collide(other->getRoot());
}

void
DeformBVHTree::self_collide()
{
	s_mdl1 = _mdl;
	s_mdl2 = _mdl;

	getRoot()->self_collide();
}

BOX
DeformBVHTree::box()
{
	return getRoot()->_box;
}

inline vec3f norm(vec3f &p1, vec3f &p2, vec3f &p3)
{
	vec3f s = p2-p1;
	vec3f t = p3-p1;
	vec3f n = s.cross(t);
	return n;
}

void
DeformBVHNode::refit()
{
	if (isLeaf()) {
		_box = s_mdl->_fac_boxes[getTriID()];
	} else {
		getLeftChild()->refit();
		getRightChild()->refit();

		_box = getLeftChild()->_box + getRightChild()->_box;
	}
}

bool
DeformBVHNode::find(unsigned int id)
{
	if (isLeaf())
		return getTriID() == id;

	if (getLeftChild()->find(id))
		return true;

	if (getRightChild()->find(id))
		return true;

	return false;
}

void
DeformBVHNode::self_collide()
{
	if (isLeaf())
		return;

	getLeftChild()->self_collide();
	getRightChild()->self_collide();
	getLeftChild()->collide(getRightChild());
}

void
DeformBVHNode::collide(DeformBVHNode *other)
{
	if (isLeaf() && other->isLeaf()) {
		bool cov = s_mdl1->Covertex_F(getTriID(), other->getTriID());

		if (!cov) {
			s_mdl1->_num_box_tests++;
			if (!_box.overlaps(other->_box))
				return;

			s_mdl1->_num_tri_tests++;
			non_adj_list.push_back(non_adjacent_pair(getTriID(), other->getTriID()));
		} else {
			s_mdl1->_num_cov_tests++;
		}

		return;
	}

	s_mdl1->_num_box_tests++;
	if (!_box.overlaps(other->_box)) {
		return;
	}

	if (isLeaf()) {
		collide(other->getLeftChild());
		collide(other->getRightChild());
	} else {
		getLeftChild()->collide(other);
		getRightChild()->collide(other);
	}
}
