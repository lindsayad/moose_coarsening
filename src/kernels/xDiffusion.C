//* This file is developed from MOOSE Framework
//* https://www.mooseframework.org
//* by Dr. Yinkai Lei at NETL (yinkai.lei@netl.doe.gov)
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "xDiffusion.h"

// register object

registerMooseObject("coarseningApp", xDiffusion);

template <>
InputParameters validParams<xDiffusion>(){
   InputParameters params = validParams<Kernel>();

   params.addRequiredCoupledVar(
          "c_variable", "the order parameter of Ni volume fraction.");
   params.addRequiredParam<MaterialPropertyName>("D_Ni", "diffusivity of Ni");
   return params;
}

xDiffusion::xDiffusion(const InputParameters & parameters) 
   : DerivativeMaterialInterface<JvarMapKernelInterface<Kernel>>(parameters),
     _c(coupledValue("c_variable")),
     _c_var(coupled("c_variable")),
     _D_Ni(getMaterialProperty<Real>("D_Ni")),
     _dD_Ni(getMaterialPropertyDerivative<Real>("D_Ni", getVar("c_variable", 0)->name()))
{
}

Real xDiffusion::computeQpResidual() {
    return _grad_test[_i][_qp] * (_D_Ni[_qp] * _grad_u[_qp]); 
}

Real xDiffusion::computeQpJacobian(){
    return _grad_test[_i][_qp] * (_D_Ni[_qp] * _grad_phi[_j][_qp]);
}

Real xDiffusion::computeQpOffDiagJacobian(unsigned int jvar){
    if(jvar == _c_var)
       return _grad_test[_i][_qp] * (_phi[_j][_qp] * _dD_Ni[_qp] * _grad_u[_qp]);
    else
       return 0.0;
}

