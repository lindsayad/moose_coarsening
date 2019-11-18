#include "Plane.h"

// register object

registerMooseObject("coarseningApp", Plane);

template<>
InputParameters validParams<Plane>(){
   InputParameters params = validParams<InitialCondition>();
   params.addRequiredParam<Real>("x", "position of the plane in x-direction");

   return params;
}

Plane::Plane(const InputParameters & parameters)
   : InitialCondition(parameters),
     _x(getParam<Real>("x"))
{
}

Real Plane::value(const Point & p){
   if(p(0) > _x) return 1.0;
   else return 0.0;
}
