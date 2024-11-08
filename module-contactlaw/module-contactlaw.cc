#include "mbconfig.h"

#include <cassert>
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <iostream>
#include <iomanip>
#include <limits>

#include "module-contactlaw.h"


/* -----------------------------contactlaw  start --------------------------------------*/
/*=======================================================================================
 * Constructor and Destructor
 *=======================================================================================*/
//constructor
Contactlaw::Contactlaw (
	unsigned uLabel,
	const DofOwner *pDO,
	DataManager* pDM,
	MBDynParser& HP
)
: Elem(uLabel, flag(0)), UserDefinedElem(uLabel, pDO)
{
	// help message or no arg error
	if (HP.IsKeyWord("help")) {
		silent_cout(
			"help message\n"
			"==== Module: Contactlaw ====\n"
			"- Note: \n"
			"\ttest\n"
			"- Usage: \n"
			"\tContactlaw;\n"
			<< std::endl);
		if (!HP.IsArg()) {
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	// read Seabed object
	unsigned int uElemLabel = (unsigned int)HP.GetInt();
	pSeabed = dynamic_cast<Seabed *>(pDM->pFindElem(Elem::LOADABLE, uElemLabel));

	// read m
	if (!HP.IsKeyWord("m")) {
		silent_cerr("Contactlaw(" << GetLabel() << "): keyword \"m\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	m = HP.GetReal();

	// read n_node
	if (!HP.IsKeyWord("n_node")) {
		silent_cerr("Contactlaw(" << GetLabel() << "): keyword \"n_node\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	n_node = HP.GetReal();
	//output flag
	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));
	//export log file
	pDM->GetLogFile()
		<< "Contactlaw: " << uLabel
		<< " " << pSeabed->GetLabel()
		<< std::endl;
}

//destructor
Contactlaw::~Contactlaw (void)
{
	NO_OP;
}


/*=======================================================================================
 * Intial Assembly Process
 *=======================================================================================*/
//set number of DOF
unsigned int
Contactlaw::iGetInitialNumDof(void) const
{
	return 0;
}

//set initial value
void
Contactlaw::SetInitialValue(VectorHandler& XCurr)
{
	return;
}

//set initial assembly matrix dimension
void 
Contactlaw::InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

//calculate residual vector for initial assembly analysis
SubVectorHandler& 
Contactlaw::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	WorkVec.ResizeReset(0);
	return WorkVec;
}

//calculate Jaconbian for initial assembly analysis
VariableSubMatrixHandler&
Contactlaw::InitialAssJac(
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
Contactlaw::iGetNumDof(void) const
{
	return 0;
}

//set DOF type
DofOrder::Order
Contactlaw::GetDofType(unsigned int i) const
{

	return DofOrder::DIFFERENTIAL;
}

//set initial value
void
Contactlaw::SetValue(
	DataManager *pDM,
	VectorHandler& X,
	VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	return;
}

//print explanation of variables and equations
std::ostream&
Contactlaw::DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const
{
}
std::ostream&
Contactlaw::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{
}


//set matrix dimension
void
Contactlaw::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const //@@自由度について2節点（それぞれ3自由度）を扱いたいので6自由度用意すべきでは
{
	*piNumRows = 1;
	*piNumCols = 1;	
}

//calculate residual vector
SubVectorHandler& 
Contactlaw::AssRes(
	SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	//obtain vector dimension
	integer iNumRows;
	integer iNumCols;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	//set index from global state vector
	const integer iPositionIndex = pNode->iGetFirstPositionIndex();
	const integer iMomentumIndex = pNode->iGetFirstMomentumIndex();
	WorkVec.PutRowIndex(1, iMomentumIndex+3);

	//compare z coordinates
	doublereal z_seabed,g;
	pSeabed->get(z_seabed,g);		
	doublereal D = XCurr(iPositionIndex+3) - z_seabed;

	//calculate refrectionforece
	doublereal s;						
	doublereal F;
	//Vec3 z = [0.0,0.0,1.0]; 								

	if(D > 0.0){
		s = 0.0;
	}else{
		s = 1.0;
		
	}
	F = s*m*g/n_node;

	//set value
	WorkVec.PutCoef(1,F);
	return WorkVec;
}
/*-------------------------------------------------------------------------------------*/

//calculate Jacobian matrix
VariableSubMatrixHandler& 
Contactlaw::AssJac(
	VariableSubMatrixHandler& WorkMat,
	doublereal dCoef, 
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	//obtain vector dimension
	integer iNumRows;
	integer iNumCols;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	//set index from global matrix
  	const integer iPositionIndex = pNode->iGetFirstPositionIndex();
	const integer iMomentumIndex = pNode->iGetFirstMomentumIndex();
	WM.PutRowIndex(1, iMomentumIndex+3);
	WM.PutColIndex(1, iPositionIndex+3);

	//compare z coordinates
	doublereal z_seabed,g;
	pSeabed->get(z_seabed,g);		
	doublereal D = XCurr(iPositionIndex+3) - z_seabed;

	//calculate refrectionforece					
	doublereal ds_dx;
	doublereal dF_dx;
	//Vec3 z = [0.0,0.0,1.0]; 								

	if(D > 0.0){
		ds_dx = 0.0;
	}else{
		ds_dx = 1.0;
		
	}
	dF_dx = ds_dx*m*g / n_node;


	// set value
	WM.PutCoef( 1,  1, -dF_dx*dCoef );
	return WorkMat;	
}


/*=======================================================================================
 * Private Data
 *=======================================================================================*/
//set number of private data
unsigned int
Contactlaw::iGetNumPrivData(void) const
{
	return 0;
}


/*=======================================================================================
 * Configure runtime processing
 *=======================================================================================*/
//describe update function
void 
Contactlaw::Update(const VectorHandler& XCurr, const VectorHandler& XPrimeCurr)
{
	return;
}
//process before each iteration
void
Contactlaw::BeforePredict(VectorHandler& /* X */ ,
					VectorHandler& /* XP */ ,
					VectorHandler& /* XPrev */ ,
					VectorHandler& /* XPPrev */ ) const
{
	return;
}
//process after each iteration
void
Contactlaw::AfterPredict(VectorHandler& X, VectorHandler& XP)
{
	return;
}
//process after convergence (each time step)
void
Contactlaw::AfterConvergence(const VectorHandler& X, const VectorHandler& XP)
{
	return;
}


/*=======================================================================================
 * Output
 *=======================================================================================*/
//output file 
void
Contactlaw::Output(OutputHandler& OH) const
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
//print information of connected nodes
int
Contactlaw::iGetNumConnectedNodes(void) const
{
	return 0;
}
void
Contactlaw::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	return;
}
//output restart file
std::ostream&
Contactlaw::Restart(std::ostream& out) const
{
   	return out << "# Contactlaw (" << GetLabel() << "): not implemented yet" << std::endl;
}

/* ----------------------------- Contactlaw end -------------------------------------- */


/*=======================================================================================
 *  Module init function
 *=======================================================================================*/
extern "C"
int module_init(const char *module_name, void *pdm, void *php)
{
	bool UDEset = true;

	UserDefinedElemRead *rf = new UDERead<Contactlaw>;
	if (!SetUDE("Contactlaw", rf)) {
		delete rf;
		return false;
	}

	if (!UDEset) {
		silent_cerr("Contactlaw: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);
		return -1;
	}

	return 0;
}
