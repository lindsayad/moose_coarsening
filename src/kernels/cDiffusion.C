//* This file is developed from MOOSE Framework
//* https://www.mooseframework.org
//* by Dr. Yinkai Lei at NETL (yinkai.lei@netl.doe.gov)
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "cDiffusion.h"

// register object

registerMooseObject("coarseningApp", cDiffusion);

template <>
InputParameters validParams<cDiffusion>(){
   InputParameters params = validParams<Diffusion>();
   params.addRequiredParam<MaterialPropertyName>("M_c", "Mobility of c");
   params.addRequiredParam<MaterialPropertyName>("kappa_c", "gradient coefficient of c");
   return params;
}

cDiffusion::cDiffusion(const InputParameters & parameters) 
   : Diffusion(parameters),
     _M(getMaterialProperty<Real>("M_c")),
     _kappa(getMaterialProperty<Real>("kappa_c"))
{
}

Real cDiffusion::computeQpResidual() {
    return _M[_qp] * _kappa[_qp] * Diffusion::computeQpResidual(); 
}

Real cDiffusion::computeQpJacobian(){
    return _M[_qp] * _kappa[_qp] * Diffusion::computeQpJacobian();
}

