/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "integrator" library , Copyright (C) 2007 CMU, 2009 MIT, 2018 USC     *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code author: Jernej Barbic                                            *
 * http://www.jernejbarbic.com/vega                                      *
 *                                                                       *
 * Research: Jernej Barbic, Hongyi Xu, Yijing Li,                        *
 *           Danyong Zhao, Bohan Wang,                                   *
 *           Fun Shing Sin, Daniel Schroeder,                            *
 *           Doug L. James, Jovan Popovic                                *
 *                                                                       *
 * Funding: National Science Foundation, Link Foundation,                *
 *          Singapore-MIT GAMBIT Game Lab,                               *
 *          Zumberge Research and Innovation Fund at USC,                *
 *          Sloan Foundation, Okawa Foundation,                          *
 *          USC Annenberg Foundation                                     *
 *                                                                       *
 * This library is free software; you can redistribute it and/or         *
 * modify it under the terms of the BSD-style license that is            *
 * included with this library in the file LICENSE.txt                    *
 *                                                                       *
 * This library is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the file     *
 * LICENSE.TXT for more details.                                         *
 *                                                                       *
 *************************************************************************/

/*
  A class to timestep large sparse dynamics using standard Euler, or symplectic
  Euler.

  standard Euler:
  x_{n+1} = x_n + h * v_n
  v_{n+1} = v_n + h * (F_n / m)

  symplectic Euler
  v_{n+1} = v_n + h * (F_n / m)
  x_{n+1} = x_n + h * v_{n+1}
*/

#ifndef _EULERSPARSE_H_
#define _EULERSPARSE_H_

#include "integratorBaseSparse.h"

class PardisoSolver;
class SPOOLESSolver;
class CGSolver;

class EulerSparse : public IntegratorBaseSparse {
public:
  // constrainedDOFs is an integer array of degrees of freedom that are to be
  // fixed to zero (e.g., to permanently fix a vertex in a deformable
  // simulation) constrainedDOFs are 0-indexed (separate DOFs for x,y,z), and
  // must be pre-sorted (ascending) dampingMatrix is optional and provides
  // damping (in addition to mass damping)
  EulerSparse(int r, double timestep, SparseMatrix *massMatrix,
              ForceModel *forceModel, int symplectic = 0,
              int numConstrainedDOFs = 0, int *constrainedDOFs = NULL,
              double dampingMassCoef = 0.0, int numSolverThreads = 1);

  virtual ~EulerSparse();

  // sets q, and (optionally) qvel
  // returns 0
  virtual int SetState(double *q, double *qvel = NULL);

  virtual int DoTimestep();

protected:
  int symplectic;
  SparseMatrix *systemMatrix;
  double *bufferConstrained;

  PardisoSolver *pardisoSolver = nullptr;
  SPOOLESSolver *spoolesSolver = nullptr;
  CGSolver *jacobiPreconditionedCGSolver = nullptr;
};

#endif
