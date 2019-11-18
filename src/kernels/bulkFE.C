//* This file is developed from MOOSE Framework
//* https://www.mooseframework.org
//* by Dr. Yinkai Lei at NETL (yinkai.lei@netl.doe.gov)
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "bulkFE.h"

// register object

registerMooseObject("coarseningApp", bulkFE);

template <>
InputParameters validParams<bulkFE>(){
   InputParameters params = validParams<Kernel>();
   params.addRequiredParam<MaterialPropertyName>("M_c", "Mobility of c");
   params.addRequiredParam<MaterialPropertyName>("rho_c", "bulk FE parameter of c");
   params.addRequiredParam<MaterialPropertyName>("f_FE", "bulk Free Energy");

   return params;
}

bulkFE::bulkFE(const InputParameters & parameters) 
   : DerivativeMaterialInterface<JvarMapKernelInterface<Kernel>>(parameters),
     _M(getMaterialProperty<Real>("M_c")),
     _rho(getMaterialProperty<Real>("rho_c")),
     _df_FE(getMaterialPropertyDerivative<Real>("f_FE", _var.name())),
     _d2f_FE(getMaterialPropertyDerivative<Real>("f_FE", _var.name(), _var.name()))
{
}

Real bulkFE::computeQpResidual() {
    return _M[_qp] * _rho[_qp] * _test[_i][_qp] * _df_FE[_qp]; 
}

Real bulkFE::computeQpJacobian(){
    return _M[_qp] * _rho[_qp] * _test[_i][_qp] * _d2f_FE[_qp] * _phi[_j][_qp];
}

