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

  US Mail:       GAMMA Research Group at UNC
                       Department of Computer Science
                       Sitterson Hall, CB #3175
                       University of N. Carolina
                       Chapel Hill, NC 27599-3175

  Phone:            (919)962-1749

  EMail:              geom@cs.unc.edu; tang_m@zju.edu.cn


\**************************************************************************/

This is the Self-CCD collision detection package.

------------------------------------------------------------------------------
1. Installation Instructions:
------------------------------------------------------------------------------

Version 1.0 of Self-CCD is developed by Visual Studio 2005, so currently it
 can only be used on Windows platform.
The package is provided as static lib files (dccd.lib) and head files.
For external calling, the following directories will be used:
       inc\                  #head files
       lib\release\          #release version of dccd.lib
       lib\debug\            #debug version of dccd.lib
       make\                 #project files for Visual Stdio 2005

For rebuilding or running the demos:

Open make\Self-ccd.sln, then build and run it.

To run the demos, please download data files from following links:
download http://gamma.cs.unc.edu/SR/benchmarks/cloth_ball.plys.zip (100 MB, 
104,915,287 bytes) and put the extracted ply files into d:\data\cloth_ball.plys\

After build, run dccd-test.exe for testing.

For example, following information is output on my PC:
==========================================
CCD Information:
intersection time: 90.8777
    collecting nonadjacent: 54.7118
    processing nonadjacent: 36.0546
    processing adjacent: 0.110595
BVH and update time: 25.9455
    refit BVH: 8.07839
    update vtx, boxes for V/E/F: 17.8665
all vf tests: 886300, true vf tests: 103277, ratio: 0.116526
all ee tests: 3111810, true ee tests: 539010, ratio: 0.173214
box test: 1383337716
total: 116.824 seconds
=============================================

------------------------------------------------------------------------------
2. Directory Contents:
------------------------------------------------------------------------------

dccd\
        inc\                            #header files for external calling
        lib\                            #dccd.lib will be generated
        make\                           #Visual Studio 2005 project files
        sample\                         #test program calling dccd.lib
        lib\                            #source file for building dccd.lib

------------------------------------------------------------------------------
3. Who will need it?
------------------------------------------------------------------------------
I assume the reader is familiar with continuous collision detection.
In fact, Self-CCD is just designed for performing continuous collision detection
between deformable model. It checks both inter-object & intra-object collisions.

Features:

Procedural Representative Triangles
Orphan Set
Non-penetration Filters

Please refer following papers for detail algorithms:
a) Min Tang, Dinesh Manocha, Ruofeng Tong, Fast Continuous Collision Detection using 
 Deforming Non-Penetration Filters, In Proceedings of ACM SIGGRAPH Symposium on
 Interactive 3D Graphics and Games (i3D 2010), Washington DC, Feb. 19-21, 2010: 7-13.
b) Min Tang, Sean Curtis, Sung-Eui Yoon, Dinesh Manocha. ICCD: Interactive Continuous
 Collision Detection between Deformable Models using Connectivity-Based Culling,
 IEEE Transactions on Visualization and Computer Graphics, 2009, 15(4): 544-557.
c) Curtis, S., Tamstorf, R., and Manocha, D. 2008. Fast collision detection for 
  deformable models using representative-triangles. In Proceedings of the 2008 
  Symposium on interactive 3D Graphics and Games (Redwood City, California, February
  15 - 17, 2008). I3D 2008. ACM, New York, NY, 61-69. 

------------------------------------------------------------------------------
4. How to use it?
------------------------------------------------------------------------------

All the API are listed in inc/ccdAPI.h. An example is shown in sample/main.cpp:

a. Providing the initial configuration of the deformable model:

ccdInitModel(vtxs, tris);

Here vtxs is a list of vertices, and tris is a list of triangles (each is defined by
three indices)

b. As the model is deforming, update the model by providing current vertices:

ccdUpdateVtxs(vtxs);

c. Checking collisions:

ccdChecking(true); // true for refitting the BVH, false for rebuilding the BVH

d. Get timing and counting information:

ccdReport();

e. Clear memory:

ccdQuitModel();

f. Set callback functions for EE/VF tests:

ccdSetEECallback(EECallback);
ccdSetVFCallback(VFCallback);

They will be called then true EE/VF tests are founded.

------------------------------------------------------------------------------
5. bug report
------------------------------------------------------------------------------

We would be interested in knowing more about your application as well as any
bugs you may encounter in the collision detection library. You can
report them by sending e-mail to geom@cs.unc.edu or tang_m@zju.edu.cn.

------------------------------------------------------------------------------
6. release notes
------------------------------------------------------------------------------

Release 1.0: 2010/6/16

