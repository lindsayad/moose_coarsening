#pragma once
#include "InitialCondition.h"

class TwoSphere;

template <>
InputParameters validParams<TwoSphere>();

class TwoSphere : public InitialCondition{
public:
   TwoSphere(const InputParameters & parameters);

   virtual Real value(const Point & p) override;

private: 
   RealVectorValue _c1, _c2;
   Real _r1, _r2;   
};
