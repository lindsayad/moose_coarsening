#include "Micstr.h"

// register object

registerMooseObject("coarseningApp", Micstr);

template<>
InputParameters validParams<Micstr>(){
   InputParameters params = validParams<InitialCondition>();
   params.addRequiredParam<std::string>("file_name", "name of the .nc file");
   return params;
}

Micstr::Micstr(const InputParameters & parameters)
   : InitialCondition(parameters),
     _file_name(getParam<std::string>("file_name"))
{

    int ncid, varid, dimid;
    int retval;
    retval = nc_open(_file_name.c_str(), NC_NOWRITE, &ncid);
    if(retval != NC_NOERR) nc_error(retval);

    // get nx, ny, nz
    retval = nc_inq_dimid(ncid, "x", &dimid);
    if(retval != NC_NOERR) nc_error(retval);
    retval = nc_inq_dimlen(ncid, dimid, &_nx);
    if(retval != NC_NOERR) nc_error(retval);
    retval = nc_inq_dimid(ncid, "y", &dimid);
    if(retval != NC_NOERR) nc_error(retval);
    retval = nc_inq_dimlen(ncid, dimid, &_ny);
    if(retval != NC_NOERR) nc_error(retval);
    retval = nc_inq_dimid(ncid, "z", &dimid);
    if(retval != NC_NOERR) nc_error(retval);
    retval = nc_inq_dimlen(ncid, dimid, &_nz);
    if(retval != NC_NOERR) nc_error(retval);
  
    // get dx, dy, dz
    retval = nc_inq_varid(ncid, "dx", &varid);
    if(retval != NC_NOERR) nc_error(retval);
    retval = nc_get_var(ncid, varid, &_dx);
    if(retval != NC_NOERR) nc_error(retval);
    retval = nc_inq_varid(ncid, "dy", &varid);
    if(retval != NC_NOERR) nc_error(retval);
    retval = nc_get_var(ncid, varid, &_dy);
    if(retval != NC_NOERR) nc_error(retval);
    retval = nc_inq_varid(ncid, "dz", &varid);
    if(retval != NC_NOERR) nc_error(retval);
    retval = nc_get_var(ncid, varid, &_dz);
    if(retval != NC_NOERR) nc_error(retval);
    
    // get microstr
    int *m = new int[_nx*_ny*_nz];
    retval = nc_inq_varid(ncid, "microstr", &varid);
    if(retval != NC_NOERR) nc_error(retval);
    retval = nc_get_var(ncid, varid, m);
    if(retval != NC_NOERR) nc_error(retval);
         
    _micstr = new int **[_nx];
    int cnt=0;
    for(size_t i=0; i<_nx; i++){
       _micstr[i] = new int *[_ny];
       for(size_t j=0; j<_ny; j++){
          _micstr[i][j] = new int[_nz];
          for(size_t k=0; k<_nz; k++){
              _micstr[i][j][k]=m[cnt];
              cnt++;
          }
       }
    }

    delete[] m;
}

Micstr::~Micstr(){
    for(size_t i=0; i<_nx; i++){
        for(size_t j=0; j<_ny; j++){
            delete[] _micstr[i][j];
        }
        delete[] _micstr[i];
    }      
    delete[] _micstr;
}

Real Micstr::value(const Point & p){
     size_t xx=p(0)/_dx;
     if(xx>_nx-1) xx=_nx-1;
     size_t yy=p(1)/_dy;
     if(yy>_ny-1) yy=_ny-1;
     size_t zz=p(2)/_dz;
     if(zz>_nz-1) zz=_nz-1;

     Real a((double) _micstr[xx][yy][zz]);
     return a; 
}

void Micstr::nc_error(const int e){
   std::cout << "Error: " << nc_strerror(e) << std::endl;
}
