#ifndef _MAINPAGE_H_
#define _MAINPAGE_H_

/** @mainpage DENDRO: A Parallel Geometric Multigrid Library for Finite Elements on Octree Meshes 
*
* @authors <A href="http://www.cc.gatech.edu/~rahulss">Rahul S. Sampath</A>
* @authors <A href="http://www.seas.upenn.edu/~hsundar">Hari Sundar</A>
* @authors <A href="http://www.seas.upenn.edu/~adavani">Santi S. Adavani</A>
* @authors <A href="http://www.cc.gatech.edu/~ilashuk">Ilya Lashuk</A>
* @authors <A href="http://www.cc.gatech.edu/~gbiros">George Biros</A>
*
* \image html hangingNodes.bmp
* 
Dendro (from the greek word, \htmlonly &delta;&epsilon;&nu;&delta;&rho;&omicron; \endhtmlonly, for tree) is a suite of parallel algorithms for the discretization and solution of partial differential equations that require discretization of second-order elliptic operators. It supports trilinear finite element discretizations constructed using octees. The package comprises of four main modules: a templated parallel utilities module ('par'), a bottom-up octree generation and 2:1 balancing module ('oct'), a meshing module ('oda'), a geometric multigrid module ('omg'). It supports the PETSc objects 'Mat' and 'Vec' and provides interfaces to PETSc's linear and non-linear solvers. Dendro can be best viewed as an extension of PETSc's DA and DMMG modules that supports octree discretizations.

* This package uses the following libraries:
* - <A href="http://www-unix.mcs.anl.gov/petsc/petsc-2/">PETSc (Version 2.3.3)</A>  PETSc (Argonne National Laboratories) is a suite of data structures and routines for the scalable (parallel) solution of scientific applications modeled by partial differential equations.
* - <A href="http://www-unix.mcs.anl.gov/mpi/">MPI</A> MPI is a library specification for message-passing, proposed as a standard by a broadly based committee of vendors, implementors, and users.
* - <A href="http://www.cppreference.com/cppstl.html">C++ STL</A> The C++ STL (Standard Template Library) is a generic collection of class templates and algorithms that allow programmers to easily implement standard data structures like queues, lists, and stacks. 
*

Dendro was initially developed at the Computational science and engineering laboratory <A href="http://www.seas.upenn.edu/~csela">(CSELa@Penn)</A> at the University of Pennsylvania. It is currently maintained by the Computational science and engineering laboratory <A href="http://www.cc.gatech.edu/csela">(CSELa@GT)</A> at the Georgia Institute of Technology. It was supported by grants from the U.S. Department of Energy and the U.S. National Science Foundation and TeraGrid resources provided by the Pittsburgh Supercomputing Center, the National Center for Supercomputing Applications and the Texas Advanced Computing Center. 

Copyright Notification

Copyright (C) 2008 Rahul S. Sampath, Hari Sundar, Santi S. Adavani, Ilya Lashuk and George Biros   
                                                                         
This program is a free software; you can redistribute it and/or modify  
it under the terms of the GNU General Public License as published by  
the Free Software Foundation; either version 2 of the License, or     
(at your option) any later version.                                   
                                                                         
This program is distributed in the hope that it will be useful,       
but WITHOUT ANY WARRANTY; without even the implied warranty of        
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         
GNU General Public License for more details.                                                             

The GNU General Public License can be found <A href="http://www.cc.gatech.edu/csela/dendro/License.gpl"> here </A>.          

Other Links
- <A href="http://www.cc.gatech.edu/csela/dendro/download-dendro.html"> Download Dendro </A>
- <A href="http://www.cc.gatech.edu/csela/dendro/dendro-papers.html"> Related Publications </A>
- <A href="http://www.cc.gatech.edu/csela/dendro/Manual.pdf"> Dendro Manual </A>
- <A href="http://www.cc.gatech.edu/csela/dendro/FAQ.html"> FAQ </A> 

To contact us please send an email to dendro.maintenance@gmail.com. Please join the google group: <A href="http://groups.google.com/group/dendro-users"> Dendro-Users </A>. We will post updates about Dendro to this group. 

*/

#endif /*_MAINPAGE_H_*/


