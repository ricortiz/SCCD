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
#include "ccdAPI.h"
#include <stdio.h>
#include "timing.h"
CBVHTimer tm;

static DeformModel *mdl;
static double g_total = 0;

ccdEETestCallback *cbFuncEE;
ccdVFTestCallback *cbFuncVF;

void ccdSetEECallback(ccdEETestCallback *funcEE)
{
	cbFuncEE = funcEE;
}

void ccdSetVFCallback(ccdVFTestCallback *funcVF)
{
	cbFuncVF = funcVF;
}

void ccdInitModel(vec3f_list &vtxs, tri_list &tris)
{
	mdl = new DeformModel(vtxs, tris);
	mdl->BuildBVH(true);
}

void ccdUpdateVtxs(vec3f_list &vtxs)
{
	double ttmp = omp_get_wtime();
	tm.startTiming(9);
	tm.startTiming(11);
	mdl->UpdateVert(vtxs);
	mdl->UpdateBoxes();
	tm.endTiming(11);
	tm.endTiming(9);
	g_total += omp_get_wtime() - ttmp;
}

void ccdChecking(bool refit)
{
	double ttmp = omp_get_wtime();
	tm.startTiming(9);
	tm.startTiming(10);

	if (!refit) {
		mdl->RebuildBVH(true);
	} else {
		mdl->RefitBVH(true);
	}

	tm.endTiming(10);
	tm.endTiming(9);

	mdl->ResetCounter();

	tm.startTiming(0);
	mdl->SelfCollide(true);
	tm.endTiming(0);

	tm.incRecord(mdl->NumBoxTest(), mdl->NumTriTest(), mdl->NumContact(), mdl->NumCCDTest(), mdl->NumCovTest(), mdl->NumLpTest(), mdl->NumCCDTrue());
	tm.incRecord(mdl->NumVFTest(), mdl->NumEETest(), mdl->NumVFTrue(), mdl->NumEETrue());
	g_total += omp_get_wtime() - ttmp;
}

void ccdQuitModel()
{
	delete mdl;
}

void ccdReport()
{
	tm.report();
	printf("total: %g seconds\n", g_total);
}
