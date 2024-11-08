#include "dataman.h"

/* seabedprop-----------------------------------------------------------------------*/

#ifndef SEABEDPROP_H
#define SEABEDPROP_H

class seabedprop
{
private:
    doublereal z_seabed;                        //海底面のZ座標
    doublereal g;
public:
    seabedprop(void);
    virtual ~seabedprop(void);
    virtual void setValue(doublereal& pz_seabed, doublereal& pg);
    virtual void get(doublereal& pz_seabed, doublereal& pg) const;
};

#endif // SEABEDPROP_H
/* ---------------------------------------------------------------------------------*/
/* seabedpropowner------------------------------------------------------------------*/

#ifndef SEABEDPROPOWNER_H
#define SEABEDPROPOWNER_H

class seabedpropowner
{
protected:
    seabedprop pSeabedprop; 
public:
    seabedpropowner(void);
    virtual ~seabedpropowner(void);

    virtual void setSeabedprop(doublereal& z_seabed, doublereal& g);
    virtual void get(doublereal& z_seabed, doublereal& g) const;
};

#endif // SEABEDPROPOWBER_H

/* ---------------------------------------------------------------------------------*/