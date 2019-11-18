//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "maskedTimeDerivative.h"


registerMooseObject("coarseningApp", maskedTimeDerivative);

template <>
InputParameters
validParams<maskedTimeDerivative>()
{
  InputParameters params = validParams<TimeDerivative>();
  params.addRequiredCoupledVar("c_YSZ", "YSZ order parameter");
  return params;
}

maskedTimeDerivative::maskedTimeDerivative(const InputParameters & parameters)
  : TimeDerivative(parameters),
    _c_YSZ(coupledValue("c_YSZ"))
{
}

Real
maskedTimeDerivative::computeQpResidual()
{
  return (1.0-_c_YSZ[_qp]) * TimeDerivative::computeQpResidual();
}

Real
maskedTimeDerivative::computeQpJacobian()
{
  return (1.0-_c_YSZ[_qp]) * TimeDerivative::computeQpJacobian();
}
