//* This file is developed from MOOSE Framework
//* https://www.mooseframework.org
//* by Dr. Yinkai Lei at NETL (yinkai.lei@netl.doe.gov)
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "KernelValue.h"
#include "JvarMapInterface.h"
#include "DerivativeMaterialInterface.h"

/**
 * The forward declaration is so that we can declare the validParams() function
 * before we actually define the class... that way the definition isn't lost
 * at the bottom of the file.
 */

// Forward Declarations
class ChemFE;

/**
 * validParams returns the parameters that this Kernel accepts / needs
 * The actual body of the function MUST be in the .C file.
 */
template <>
InputParameters validParams<ChemFE>();

/**
 * Define the Kernel for a diffusion operator that looks like:
 *
 * M kappa (\rho \partial f/\partial c, test)
 *
 * where c is an order parameter.
 */
class ChemFE : public DerivativeMaterialInterface<JvarMapKernelInterface<Kernel>>
{
public:
  /**
   * This is the constructor declaration.  This class takes a
   * string and a InputParameters object, just like other
   * Kernel-derived classes.
   */
  ChemFE(const InputParameters & parameters);

protected:
  /**
   * Responsible for computing the residual at one quadrature point.
   * This function should always be defined in the .C file.
   */
  virtual Real computeQpResidual() override;

  /**
   * Responsible for computing the diagonal block of the preconditioning matrix.
   * This is essentially the partial derivative of the residual with respect to
   * the variable this kernel operates on ("x").
   *
   * Note that this can be an approximation or linearization.  In this case it's
   * not because the Jacobian of this operator is easy to calculate.
   *
   * This function should always be defined in the .C file.
   */
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:

  // coupled x_Ni variable
  const VariableValue & _x; 
  unsigned int _x_var;

  // equilibrium mole fraction of Ni in solid and pore, RT/Vm
  const MaterialProperty<Real> & _M;
  const MaterialProperty<Real> & _x_seq;
  const MaterialProperty<Real> & _x_peq;
  const MaterialProperty<Real> & _RT_Vm;
  const MaterialProperty<Real> & _dh_FE;
  const MaterialProperty<Real> & _d2h_FE;
  const MaterialProperty<Real> & _h_c;
  const MaterialProperty<Real> & _dh_c;

  // get Delta xs
  Real get_Delta_xs(const Real x, const Real xs, const Real xp, const Real h);
  Real get_dx(const Real xs, const Real xp);
  Real get_gf(const Real xs, const Real xp, const Real RTV);
};

