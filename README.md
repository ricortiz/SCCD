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

------------------------------------------------------------------------------
  The authors may be contacted via:

  US Mail:       GAMMA Research Group at UNC
                       Department of Computer Science
                       Sitterson Hall, CB #3175
                       University of N. Carolina
                       Chapel Hill, NC 27599-3175

  Phone:            (919)962-1749

  EMail:              geom@cs.unc.edu; tang_m@zju.edu.cn

------------------------------------------------------------------------------

## This is the Self-CCD collision detection package.

### Building Instructions:
 
 Use cmake to build the library. You can obtain cmake from [here](http://www.cmake.org). CMake
 allows you to generate projects files for Visual Studio(.sln files) and Makefiles for unix based build systems.

 Follow a quick step-by-step specific for SCCD. In the example we are running windows and want Visual Studio X project files.
 
 **Step 1**

  First you open the CMake GUI (you probably got an icon on your desktop, so just double click on this icon)
  Specify the path for the SCCD source code. This is the folder location on your hard drive where you choice to check-out SCCD from the git repository.
  Specify a path where you want the binaries to be. This folder will contain the generated Visual Studio project files.

 **Step 2**

  Next click on the configure button. CMake now ask you want type of generator you want. This is where you pick what kind of makefiles you want in the build-folder you specified in the previous step. Select the genarator corresponding to the Visual Studio version installed in your computer.
 
 **Step 3**

  CMake might ask you if it should create the build folder you specified. Just accept.
  
  Now you have to wait for cmake to configure your project. 
  Now you should see some Cache values marked with a red-color.

 **Step 4**

  The red color simply mean that CMake asks you to decide whether you want to override/set any of the shown values.  For convenience I have made a quick explanation of some of the most magical cache values in the table below. If you do not care about what they are used for simply skip the table.

  
 **Step 5**

  Press the configure button once more, and wait while CMake process everything a second time around.
  Now all the cache values should be colored in gray, this implies that everything is good to go and we can start creating the content of the build folder.
  Finally press the OK-button. Now CMake will create the content of the build folder and automatically close the CMake application.
  After having run CMake you simple go to the build folder you specified and look for the file SCCD.sln, double click on this file and you are up and running, ready to code:-)
  
####To run the demos...

Please download data files from following links:
download http://gamma.cs.unc.edu/SR/benchmarks/cloth_ball.plys.zip (100 MB, 
104,915,287 bytes) and put the extracted ply files into d:\data\cloth_ball.plys\

After build, run dccd-test.exe for testing.

For example, following information is output on my PC:

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

###Directory Contents:

    SDCC\
      inc\                            #header files
      src\                            #source files
      sample\                         #test program
      CMakeLists.txt                  #cmake file


###Who will need it?

I assume the reader is familiar with continuous collision detection.
In fact, Self-CCD is just designed for performing continuous collision detection
between deformable model. It checks both inter-object & intra-object collisions.

Features:

Procedural Representative Triangles
Orphan Set
Non-penetration Filters

Please refer following papers for detail algorithms:

    a) Min Tang, Dinesh Manocha, Ruofeng Tong, Fast Continuous Collision Detection using Deforming Non-Penetration Filters, In Proceedings of ACM SIGGRAPH Symposium on Interactive 3D Graphics and Games (i3D 2010), Washington DC, Feb. 19-21, 2010: 7-13.
    b) Min Tang, Sean Curtis, Sung-Eui Yoon, Dinesh Manocha. ICCD: Interactive Continuous Collision Detection between Deformable Models using Connectivity-Based Culling, IEEE Transactions on Visualization and Computer Graphics, 2009, 15(4): 544-557.
    c) Curtis, S., Tamstorf, R., and Manocha, D. 2008. Fast collision detection for deformable models using representative-triangles. In Proceedings of the 2008 Symposium on interactive 3D Graphics and Games (Redwood City, California, February 15 - 17, 2008). I3D 2008. ACM, New York, NY, 61-69. 

###How to use it?

All the API are listed in inc/ccdAPI.h. An example is shown in sample/main.cpp:

 a. Providing the initial configuration of the deformable model:

    ccdInitModel(vtxs, tris);

  Here vtxs is a list of vertices, and tris is a list of triangles (each is defined by three indices)

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

###bug report

We would be interested in knowing more about your application as well as any
bugs you may encounter in the collision detection library. You can
report them by sending e-mail to geom@cs.unc.edu or tang_m@zju.edu.cn.

###release notes

Release 1.0: 2010/6/16

