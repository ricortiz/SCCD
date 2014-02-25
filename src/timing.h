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
#include <omp.h>

# define	TIMING_BEGIN \
	{double tmp_timing_start = omp_get_wtime();

# define	TIMING_END(message) \
	{double tmp_timing_finish = omp_get_wtime();\
	double  tmp_timing_duration = tmp_timing_finish - tmp_timing_start;\
	printf("%s: %2.5f seconds\n", (message), tmp_timing_duration);}}

# define TIMER_NUM     15

class CBVHTimer {
	double _rets[TIMER_NUM];
	double _tmps[TIMER_NUM];

	int _count;
	long long _box_test, _tri_test, _contact, _ccd_test, _cov_test, _lp_test, _ccd_true;
	long long _vf_test, _vf_true, _ee_test, _ee_true;

public:
	CBVHTimer() {
		resetTiming();
	}

	void updatTiming() {
			_count++;
	}

	void startTiming(int i) {
		assert(i >= 0 && i<TIMER_NUM);
		_tmps[i] = omp_get_wtime();
	}

	double endTiming(int i) {
		assert(i>=0 && i<TIMER_NUM);
		double dur = omp_get_wtime()-_tmps[i];
		_rets[i] += dur;
		return dur;
	}

	void resetTiming();
	void incRecord(int box_test, int tri_test, int contact, int ccd_test, int cov_test, int lp_test, int ccd_true);
	void incRecord(int vf_test, int ee_test, int vf_true, int ee_true); 
	void report();
};

