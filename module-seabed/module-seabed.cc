#include "mbconfig.h"

#include <cassert>
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <iostream>
#include <iomanip>
#include <limits>

#include "module-seabed.h"
#include "seabedprop.h"

/* ----------------------------- Seabed start --------------------------------------*/

/*=======================================================================================
 * Constructor and Destructor
 *=======================================================================================*/
//constructor
Seabed::Seabed (
	unsigned uLabel,
	const DofOwner *pDO,
	DataManager* pDM,
	MBDynParser& HP
)
: Elem(uLabel, flag(0)), UserDefinedElem(uLabel, pDO), seabedpropowner()
{
	// help message or no arg error
	if (HP.IsKeyWord("help")) {
		silent_cout(
			"help message\n"
			"==== Module: Seabed ====\n"
			"- Note: \n"
			"\tTest, \n"
			"- Usage: \n"
			"\tSeabed, z_node, z_seabed;\n"
			<< std::endl);
		
		if (!HP.IsArg()) {
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	doublereal z_seabed = HP.GetReal();
	doublereal g = HP.GetReal();
	pSeabedprop.setValue(z_seabed,g);

	//output flag
	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));
	//export log file
	pDM->GetLogFile()
		<< "Seabed: " << uLabel
		<< std::endl;
}

//destructor
Seabed::~Seabed (void)
{
	NO_OP;
}


/*=======================================================================================
 * Intial Assembly Process
 *=======================================================================================*/
//set number of DOF
unsigned int
Seabed::iGetInitialNumDof(void) const
{
	return 0;
}

//set initial value
void
Seabed::SetInitialValue(VectorHandler& XCurr)
{
	return;
}

//set initial assembly matrix dimension
void 
Seabed::InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

//calculate residual vector for initial assembly analysis
SubVectorHandler& 
Seabed::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	WorkVec.ResizeReset(0);
	return WorkVec;
}

//calculate Jaconbian for initial assembly analysis
VariableSubMatrixHandler&
Seabed::InitialAssJac(
	VariableSubMatrixHandler& WorkMat, 
	const VectorHandler& XCurr)
{
	WorkMat.SetNullMatrix();
	return WorkMat;
}


/*=======================================================================================
 * Initial Value Problem
 *=======================================================================================*/
//set number of DOF
unsigned int
Seabed::iGetNumDof(void) const
{
	return 0;
}

//set DOF type
DofOrder::Order
Seabed::GetDofType(unsigned int i) const
{
	return DofOrder::DIFFERENTIAL;
}

//set initial value
void
Seabed::SetValue(
	DataManager *pDM,
	VectorHandler& X,
	VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	return;
}

//set matrix dimension
void
Seabed::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;	
}

//calculate residual vector
SubVectorHandler& 
Seabed::AssRes(
	SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	WorkVec.ResizeReset(0);
	return WorkVec;
}

//calculate Jacobian matrix
VariableSubMatrixHandler& 
Seabed::AssJac(
	VariableSubMatrixHandler& WorkMat,
	doublereal dCoef, 
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	WorkMat.SetNullMatrix();
	return WorkMat;
}


/*=======================================================================================
 * Private Data
 *=======================================================================================*/
//set number of private data
unsigned int
Seabed::iGetNumPrivData(void) const
{
	return 0;
}


/*=======================================================================================
 * Configure runtime processing
 *=======================================================================================*/
//describe update function
void 
Seabed::Update(const VectorHandler& XCurr, const VectorHandler& XPrimeCurr)
{
	return;
}

//process before each iteration
void
Seabed::BeforePredict(VectorHandler& /* X */ ,
					   VectorHandler& /* XP */ ,
					   VectorHandler& /* XPrev */ ,
					   VectorHandler& /* XPPrev */ ) const
{
	return;
}

//process after each iteration
void
Seabed::AfterPredict(VectorHandler& X, VectorHandler& XP)
{
	return;
}

//process after convergence (each time step)
void
Seabed::AfterConvergence(const VectorHandler& X, const VectorHandler& XP)
{
	return;
}


/*=======================================================================================
 * Output
 *=======================================================================================*/
//output file 
void
Seabed::Output(OutputHandler& OH) const
{
	if (bToBeOutput()) {
		if (OH.UseText(OutputHandler::LOADABLE)) {
			OH.Loadable() << GetLabel()
				<< std::endl;
		}
	}
}


/*=======================================================================================
 * etc
 *=======================================================================================*/

//output restart file
std::ostream&
Seabed::Restart(std::ostream& out) const
{
   	return out << "# Seabed (" << GetLabel() << "): not implemented yet" << std::endl;
}
/* ----------------------------- seabed end -------------------------------------- */


/*=======================================================================================
 *  Module init function
 *=======================================================================================*/
extern "C"
int module_init(const char *module_name, void *pdm, void *php)
{
	bool UDEset = true;

	UserDefinedElemRead *rf = new UDERead<Seabed>;
	if (!SetUDE("Seabed", rf)) {
		delete rf;
		return false;
	}

	if (!UDEset) {
		silent_cerr("Seabed: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);
		return -1;
	}

	return 0;
}
