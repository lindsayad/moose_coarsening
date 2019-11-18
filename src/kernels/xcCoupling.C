//* This file is developed from MOOSE Framework
//* https://www.mooseframework.org
//* by Dr. Yinkai Lei at NETL (yinkai.lei@netl.doe.gov)
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "xcCoupling.h"

// register object

registerMooseObject("coarseningApp", xcCoupling);

template <>
InputParameters validParams<xcCoupling>(){
   InputParameters params = validParams<Kernel>();

   params.addRequiredCoupledVar(
          "c_variable", "the order parameter of Ni volume fraction.");
   params.addRequiredParam<MaterialPropertyName>("D_Ni", "diffusivity of Ni");
   params.addRequiredParam<MaterialPropertyName>("h_c", "interpolation function of x");
   params.addRequiredParam<MaterialPropertyName>("x_seq", "Equilibrium mole fraction in Ni");
   params.addRequiredParam<MaterialPropertyName>("x_peq", "Equilibrium mole fraction in Pore");
   return params;
}

xcCoupling::xcCoupling(const InputParameters & parameters) 
   : DerivativeMaterialInterface<JvarMapKernelInterface<Kernel>>(parameters),
     _c(coupledValue("c_variable")),
     _grad_c(coupledGradient("c_variable")),
     _c_var(coupled("c_variable")),
     _D_Ni(getMaterialProperty<Real>("D_Ni")),
     _dD_Ni(getMaterialPropertyDerivative<Real>("D_Ni", getVar("c_variable", 0)->name())),
     _dh_c(getMaterialPropertyDerivative<Real>("h_c", getVar("c_variable", 0)->name())),
     _d2h_c(getMaterialPropertyDerivative<Real>("h_c", getVar("c_variable", 0)->name(), getVar("c_variable", 0)->name())),
     _x_seq(getMaterialProperty<Real>("x_seq")),
     _x_peq(getMaterialProperty<Real>("x_peq"))
{
}

Real xcCoupling::computeQpResidual() {
    return -_grad_test[_i][_qp] * (_D_Ni[_qp] * _dh_c[_qp] * (_x_seq[_qp] - _x_peq[_qp]) * _grad_c[_qp]); 
}

Real xcCoupling::computeQpJacobian(){
    return 0.0;
}

Real xcCoupling::computeQpOffDiagJacobian(unsigned int jvar){
    if(jvar == _c_var)
       return -_grad_test[_i][_qp] * (_x_seq[_qp] - _x_peq[_qp]) * ( _D_Ni[_qp] * _phi[_j][_qp] * _d2h_c[_qp] * _grad_c[_qp] 
                                                    + _D_Ni[_qp] * _dh_c[_qp] * _grad_phi[_j][_qp]
                                                    + _phi[_j][_qp] * _dD_Ni[_qp] * _dh_c[_qp] * _grad_c[_qp] );
    else
       return 0.0;
}


