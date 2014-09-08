#ifndef DcBranchBase_hpp_
#define DcBranchBase_hpp_

//#############################################################################
// This file is modified from COIN-OR's ALPS library example file
// AbcBranchBase.h
//#############################################################################

#include <string>
#include <vector>

class OsiSolverInterface;

class DcModel;
class DcNode;
class DcNodeDesc;
class DcBranchingObject;

//#############################################################################

/** Abstract branching decision base class

    In the abstract, an DcBranchDecision object is expected to be able to
    compare two possible branching choices.

    The #betterBranch() method is the crucial routine. It is expected to be
    able to compare two integer variables.
*/
class DcBranchDecision {
 public:
    /// Default Constructor
    DcBranchDecision ();

    /// Destructor
    virtual ~DcBranchDecision();

    /// Clone
    virtual DcBranchDecision * clone() const = 0;

    /// Initialize <i>e.g.</i> before starting to choose a branch at a node
    virtual void initialize(DcModel * model) = 0;

    /** \brief Compare two branching objects (current just integer variables).
	Return nonzero if branching using \p thisOne is better than
	branching using \p bestSoFar.

	If \p bestSoFar is NULL, the routine should return a nonzero value.
	This routine is used only after strong branching.

	It is now reccommended that bestBranch is used - see below.
	This has been left for compatibility.
    */
    virtual int
	betterBranch(int thisOne,
		     int bestSoFar,
		     double changeUp,
		     int numberInfeasibilitiesUp,
		     double changeDown,
		     int numberInfeasibilitiesDown) = 0 ;

    /** \brief Compare N branching objects. Return index of best
	and sets way of branching in chosen object.

	This routine is used only after strong branching.
	This is reccommended version as it can be more sophisticated
    */
    virtual int	bestBranch ( DcModel* model,
			     int* objects,
			     int numberObjects,
			     int numberUnsatisfied,
			     double * changeUp,
			     int * numberInfeasibilitiesUp,
			     double * changeDown,
			     int * numberInfeasibilitiesDown,
			     double objectiveValue );

 private:

    /// Assignment is illegal
    DcBranchDecision & operator=(const DcBranchDecision& rhs);

};

#endif
