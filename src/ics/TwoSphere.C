#include "TwoSphere.h"

// register object

registerMooseObject("coarseningApp", TwoSphere);

template<>
InputParameters validParams<TwoSphere>(){
   InputParameters params = validParams<InitialCondition>();
   params.addRequiredParam<RealVectorValue>("center1", "center of the first sphere");
   params.addRequiredParam<RealVectorValue>("center2", "center of the second sphere");
   params.addRequiredParam<Real>("radius1", "radius of the first sphere");
   params.addRequiredParam<Real>("radius2", "radius of the second sphere");

   return params;
}

TwoSphere::TwoSphere(const InputParameters & parameters)
   : InitialCondition(parameters),
     _c1(getParam<RealVectorValue>("center1")),
     _c2(getParam<RealVectorValue>("center2")),
     _r1(getParam<Real>("radius1")),
     _r2(getParam<Real>("radius2"))
{
}

Real TwoSphere::value(const Point & p){
   Real r12=(p(0)-_c1(0))*(p(0)-_c1(0))
           +(p(1)-_c1(1))*(p(1)-_c1(1))
           +(p(2)-_c1(2))*(p(2)-_c1(2));
   Real r22=(p(0)-_c2(0))*(p(0)-_c2(0))
           +(p(1)-_c2(1))*(p(1)-_c2(1))
           +(p(2)-_c2(2))*(p(2)-_c2(2));
   if(r12<=_r1*_r1 or r22<=_r2*_r2) return 1.0;
   else return 0.0;
}
