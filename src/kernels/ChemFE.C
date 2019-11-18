//* This file is developed from MOOSE Framework
//* https://www.mooseframework.org
//* by Dr. Yinkai Lei at NETL (yinkai.lei@netl.doe.gov)
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ChemFE.h"

// register object

registerMooseObject("coarseningApp", ChemFE);

template <>
InputParameters validParams<ChemFE>(){
   InputParameters params = validParams<Kernel>();

   params.addRequiredCoupledVar("x_variable", "the order parameter of Ni mole fraction.");
   params.addRequiredParam<MaterialPropertyName>("M_c", "Mobility of c");
   params.addRequiredParam<MaterialPropertyName>("x_seq", "Equilibrium mole fraction in Ni");
   params.addRequiredParam<MaterialPropertyName>("x_peq", "Equilibrium mole fraction in Pore");
   params.addRequiredParam<MaterialPropertyName>("RTVm", "RT/Vm, Vm is the molar volume of Ni");
   params.addRequiredParam<MaterialPropertyName>("h_FE", "Interpolation of free energy");
   params.addRequiredParam<MaterialPropertyName>("h_c", "Interpolation of mole fraction");

   return params;
}

ChemFE::ChemFE(const InputParameters & parameters) 
   : DerivativeMaterialInterface<JvarMapKernelInterface<Kernel>>(parameters),
     _x(coupledValue("x_variable")),
     _x_var(coupled("x_variable")),
     _M(getMaterialProperty<Real>("M_c")),
     _x_seq(getMaterialProperty<Real>("x_seq")),
     _x_peq(getMaterialProperty<Real>("x_peq")),
     _RT_Vm(getMaterialProperty<Real>("RTVm")),
     _dh_FE(getMaterialPropertyDerivative<Real>("h_FE", _var.name())),
     _d2h_FE(getMaterialPropertyDerivative<Real>("h_FE", _var.name(), _var.name())),
     _h_c(getMaterialProperty<Real>("h_c")),
     _dh_c(getMaterialPropertyDerivative<Real>("h_c", _var.name()))
{
}

Real ChemFE::computeQpResidual() {
    Real _gf=get_gf(_x_seq[_qp], _x_peq[_qp], _RT_Vm[_qp]);
    return _M[_qp] * _test[_i][_qp] * _dh_FE[_qp] * _gf * get_Delta_xs(_x[_qp], _x_seq[_qp], _x_peq[_qp], _h_c[_qp]);
}

Real ChemFE::computeQpJacobian(){
    Real _dxeq=get_dx(_x_seq[_qp], _x_peq[_qp]);
    Real _gf=get_gf(_x_seq[_qp], _x_peq[_qp], _RT_Vm[_qp]);
    return _M[_qp] * _test[_i][_qp] * _phi[_j][_qp] * _gf 
              * (_d2h_FE[_qp] * get_Delta_xs(_x[_qp], _x_seq[_qp], _x_peq[_qp], _h_c[_qp])
                 - _dh_FE[_qp] * _dh_c[_qp] * _dxeq );
}

Real ChemFE::computeQpOffDiagJacobian(unsigned int jvar){
    Real _gf=get_gf(_x_seq[_qp], _x_peq[_qp], _RT_Vm[_qp]);
    if(jvar == _x_var)
       return _M[_qp] * _test[_i][_qp] * _phi[_j][_qp] * _gf * _dh_FE[_qp];
    else
       return 0.0;
}

Real ChemFE::get_Delta_xs(const Real x, const Real xs, const Real xp, const Real h){
    return x - (xp + h * ( xs - xp ) );
}

Real ChemFE::get_dx(const Real xs, const Real xp){
    return xs - xp;
}

Real ChemFE::get_gf(const Real xs, const Real xp, const Real RTV){
    return ( xp - xs ) / xs / ( 1.0 - xs ) * RTV;
}
