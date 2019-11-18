#pragma once
#include "InitialCondition.h"

class Plane;

template <>
InputParameters validParams<Plane>();

class Plane : public InitialCondition{
public:
   Plane(const InputParameters & parameters);

   virtual Real value(const Point & p) override;

private:
   Real _x;
};
