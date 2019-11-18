//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "TimeDerivative.h"

// Forward Declarations
class maskedTimeDerivative;

template <>
InputParameters validParams<maskedTimeDerivative>();

class maskedTimeDerivative : public TimeDerivative
{
public:
  maskedTimeDerivative(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;

private:
  const VariableValue & _c_YSZ;
};
