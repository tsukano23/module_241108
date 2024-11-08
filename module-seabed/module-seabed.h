#ifndef MODULE_SEABED_H
#define MODULE_SEABED_H

#include "dataman.h"
#include "userelem.h"
#include"seabedprop.h"

class Seabed
: virtual public Elem, public UserDefinedElem, public seabedpropowner
{
public:
	/*===================================================================
	 * Constructor and Destructor
	 *===================================================================*/
	//constructor
	Seabed(unsigned uLabel, const DofOwner *pDO,
		DataManager* pDM, MBDynParser& HP);
	//destructor
	virtual ~Seabed(void);


	/*===================================================================
	 * Intial Assembly Process
	 *===================================================================*/
	/*-------------------------------------------------------------------
	 * Configure private DOF for intial assembly process
	 *-------------------------------------------------------------------*/
	//set number of DOF
	virtual unsigned int iGetInitialNumDof(void) const;
	//set initial value
	virtual void SetInitialValue(VectorHandler& XCurr);

	/*-------------------------------------------------------------------
	 * Configure contribution for intial assembly iteration
	 *-------------------------------------------------------------------*/
	//set initial assembly matrix dimension
	virtual void 
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
	//calculate residual vector for initial assembly analysis
   	SubVectorHandler& 
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);
	//calculate Jaconbian for initial assembly analysis
   	VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat, 
		      const VectorHandler& XCurr);


	/*===================================================================
	 * Initial Value Problem
	 *===================================================================*/
	/*-------------------------------------------------------------------
	 * Configure private DOF
	 *-------------------------------------------------------------------*/
	//set number of DOF
	virtual unsigned int iGetNumDof(void) const;
	//set DOF type
	virtual DofOrder::Order GetDofType(unsigned int i) const;
	//set initial value
	void SetValue(DataManager *pDM, VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph);
	

	//set matrix dimension
	virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
	//calculate residual vector, b
	SubVectorHandler& 
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr, 
		const VectorHandler& XPrimeCurr);
	//calculate Jacobian matrix, A
	VariableSubMatrixHandler& 
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef, 
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);


	/*===================================================================
	 * Private Data
	 *===================================================================*/
	/*-------------------------------------------------------------------
	 * Configure private data
	 *-------------------------------------------------------------------*/
	//set number of private data
	virtual unsigned int iGetNumPrivData(void) const;

	/*-------------------------------------------------------------------
	 * Configure runtime processing
	 *-------------------------------------------------------------------*/
	//describe update function
	virtual void 
	Update(const VectorHandler& XCurr, const VectorHandler& XPrimeCurr);
	//process before each iteration
	virtual void
	BeforePredict(VectorHandler& /* X */ ,
					VectorHandler& /* XP */ ,
					VectorHandler& /* XPrev */ ,
					VectorHandler& /* XPPrev */ ) const;
	//process after each iteration
	virtual void
	AfterPredict(VectorHandler& X, VectorHandler& XP);
	//process after convergence (each time step)
	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP);


	/*===================================================================
	 * Output
	 *===================================================================*/
	//output file 
	virtual void Output(OutputHandler& OH) const;


	/*===================================================================
	 * etc
	 *===================================================================*/
	//print information of connected nodes
	//virtual int iGetNumConnectedNodes(void) const;
	//virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
	//output restart file
	virtual std::ostream& Restart(std::ostream& out) const;
};

#endif // MODULE_SEABED_H
