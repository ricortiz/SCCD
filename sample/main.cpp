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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#pragma warning(disable: 4996)

#include "ccdAPI.h"

char DATA_PATH[512]="d:\\data\\cloth_ball.plys\\cloth_ball";
int NUM_FRAME=94;
float MODEL_SCALE = 1.f;
int SLICES = 5;

extern void BuildSession(char *fname, unsigned int num_frame, float ply_scale,
				  unsigned int &vtx_num, vec3f *&all_vtxs, vec3f_list &vtxs, tri_list &tris);

inline vec3f interp(const vec3f &p1, const vec3f &p2, float t)
{
	return p1*(1-t)+p2*t;
}

void
GetCurrentVtxs(vec3f *all_vtxs, unsigned int vtx_num, unsigned int prev, unsigned int next, float t, vec3f_list &vtxs)
{
	vec3f *prev_pt = all_vtxs+prev*vtx_num;
	vec3f *next_pt = all_vtxs+next*vtx_num;

	for (int i=0; i<vtx_num; i++) {
		vtxs[i] = interp(prev_pt[i], next_pt[i], t);
	}
}

void EECallback(unsigned int e1_v1, unsigned e1_v2,
				unsigned int e2_v1, unsigned int e2_v2, float t) {
	printf("EE result: e1(%d, %d), e2(%d, %d) @ t=%f\n", e1_v1, e1_v2, e2_v1, e2_v2, t);
}
void VFCallback(unsigned int vid, unsigned int fid, float t) {
	printf("VF result: v=%d, f=%d @ t=%f\n", vid, fid, t);
}

int main(int argc, char **argv)
{
	unsigned int vtx_num;
	vec3f *all_vtxs;
	vec3f_list vtxs;
	tri_list tris;

	BuildSession(DATA_PATH, NUM_FRAME, MODEL_SCALE,
		vtx_num, all_vtxs, vtxs, tris);

	ccdInitModel(vtxs, tris);

	float circle = SLICES*NUM_FRAME;
	for (int i=0; i<circle; i++) {
		float t= i/circle*NUM_FRAME;
		unsigned int prev_frame = ((int)t)%NUM_FRAME;
		unsigned int next_frame = (prev_frame+1)%NUM_FRAME;

		if (next_frame < prev_frame)
			break;

		GetCurrentVtxs(all_vtxs, vtx_num, prev_frame, next_frame, t-prev_frame, vtxs);
		ccdUpdateVtxs(vtxs);

//		printf("i = %d\n", i);
		if (i == 300) {
			ccdSetEECallback(EECallback);
			ccdSetVFCallback(VFCallback);
		} else {
			ccdSetEECallback((ccdEETestCallback*)NULL);
			ccdSetVFCallback((ccdVFTestCallback*)NULL);
		}

		ccdChecking(true);
	}

	ccdReport();
	ccdQuitModel();

	return 0;
}
