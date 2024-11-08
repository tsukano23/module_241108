#include <cassert>
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <iostream>
#include <iomanip>
#include <limits>

#include "mbconfig.h"
#include "seabedprop.h"

/* seabedprop-----------------------------------------------------------------------*/
seabedprop::seabedprop(void)
{
    NO_OP;
}

seabedprop::~seabedprop(void)
{
    NO_OP;
}

void
seabedprop::setValue(doublereal& pz_seabed, doublereal& pg)
{   
    z_seabed = pz_seabed;
    g = pg;
}

void
seabedprop::get(doublereal& pz_seabed, doublereal& pg) const
{
    pz_seabed   = z_seabed;
    pg          = g;
}

/* ---------------------------------------------------------------------------------*/

/* seabedpropowner------------------------------------------------------------------*/
seabedpropowner::seabedpropowner(void)
{
    NO_OP;
}

seabedpropowner::~seabedpropowner(void)
{
    NO_OP;   
}

void
seabedpropowner::setSeabedprop(doublereal& z_seabed, doublereal& g)
{
    pSeabedprop.setValue(z_seabed,g);
}

void
seabedpropowner::get(doublereal& z_seabed, doublereal& g) const
{
    pSeabedprop.get(z_seabed,g);
}
/* ---------------------------------------------------------------------------------*/