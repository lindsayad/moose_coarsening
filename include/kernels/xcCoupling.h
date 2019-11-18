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
class xcCoupling;

/**
 * validParams returns the parameters that this Kernel accepts / needs
 * The actual body of the function MUST be in the .C file.
 */
template <>
InputParameters validParams<xcCoupling>();

/**
 * Define the Kernel for a diffusion operator that looks like:
 *
 * (dh/dc (x_seq-x_peq) \nabla c, \nabla test)
 *
 * where c is an order parameter.
 */
class xcCoupling : public DerivativeMaterialInterface<JvarMapKernelInterface<Kernel>>
{
public:
  /**
   * This is the constructor declaration.  This class takes a
   * string and a InputParameters object, just like other
   * Kernel-derived classes.
   */
  xcCoupling(const InputParameters & parameters);

protected:
  /**
   * Responsible for computing the residual at one quadrature point.
   * This function should always be defined in the .C file.
   */
  virtual Real computeQpResidual() override;

  /**
   * Responsible for computing the diagonal block of the preconditioning matrix.
   * This is essentially the partial derivative of the residual with respect to
   * the variable this kernel operates on ("u").
   *
   * Note that this can be an approximation or linearization.  In this case it's
   * not because the Jacobian of this operator is easy to calculate.
   *
   * This function should always be defined in the .C file.
   */
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:
  // volume fraction order parameter c and its gradient
  const VariableValue & _c;
  const VariableGradient & _grad_c;
  unsigned int _c_var;

  // difference between the equilibrium molar fraction of x in solid and pore
  const MaterialProperty<Real> & _D_Ni;
  const MaterialProperty<Real> & _dD_Ni;
  const MaterialProperty<Real> & _dh_c;
  const MaterialProperty<Real> & _d2h_c;
  const MaterialProperty<Real> & _x_seq;
  const MaterialProperty<Real> & _x_peq;

};

