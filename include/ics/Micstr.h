#pragma once
#include "InitialCondition.h"
#include <string>
#include <iostream>
#include <cmath>
#include <netcdf.h>

class Micstr;

template <>
InputParameters validParams<Micstr>();

class Micstr : public InitialCondition{
public:
   Micstr(const InputParameters & parameters);

   ~Micstr();

   virtual Real value(const Point & p) override;

private:
   
   void nc_error(const int e);
   
   std::string _file_name;
   
   int *** _micstr=NULL;
   size_t _nx, _ny, _nz;
   double _dx, _dy, _dz;
};
