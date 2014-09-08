#ifndef DcCutGenerator_h_
#define DcCutGenerator_h_

//#############################################################################
// This file is modified from COIN-OR's ALPS library example file
// AbcCutGenerator.h
//#############################################################################

#include "OsiSolverInterface.hpp"
#include "OsiCuts.hpp"

class DcModel;
class OsiRowCut;
class OsiRowCutDebugger;
class CglCutGenerator;

//#############################################################################

/** Interface between Dc and Cut Generation Library.

    \c DcCutGenerator is intended to provide an intelligent interface between
    Dc and the cutting plane algorithms in the CGL. A \c DcCutGenerator is
    bound to a \c CglCutGenerator and to an \c DcModel. It contains parameters
    which control when and how the \c generateCuts method of the
    \c CglCutGenerator will be called.

    The builtin decision criteria available to use when deciding whether to
    generate cuts are limited: every <i>X</i> nodes, when a solution is found,
    and when a subproblem is found to be infeasible. The idea is that the class
    will grow more intelligent with time.

    \todo Add a pointer to function member which will allow a client to install
        their own decision algorithm to decide whether or not to call the CGL
        \p generateCuts method. Create a default decision method that looks
	at the builtin criteria.

    \todo It strikes me as not good that generateCuts contains code specific to
	individual CGL algorithms. Another set of pointer to function members,
	so that the client can specify the cut generation method as well as
	pre- and post-generation methods? Taken a bit further, should this
	class contain a bunch of pointer to function members, one for each
	of the places where the cut generator might be referenced?
	Initialization, root node, search tree node, discovery of solution,
	and termination all come to mind. Initialization and termination would
	also be useful for instrumenting sbb.
*/

class DcCutGenerator  {

 public:

    /** \name Generate Cuts */
    //@{
    /** Generate cuts for the client model.

	Evaluate the state of the client model and decide whether to
	generate cuts. The generated cuts are inserted into and returned
	in the collection of cuts \p cs.

	If \p fullScan is true, the generator is obliged to call the CGL
	\c generateCuts routine.  Otherwise, it is free to make a local
	decision. The current implementation uses \c whenCutGenerator_
	to decide.

	The routine returns true if reoptimisation is needed (because the
	state of the solver interface has been modified).
    */
    bool generateCuts( OsiCuts &cs, bool fullScan);
    //@}


    /**@name Constructors and destructors */
    //@{
    /// Default constructor
    DcCutGenerator ();

    /// Normal constructor
    DcCutGenerator(DcModel * model,CglCutGenerator * generator,
		    int howOften=1, const char * name=NULL,
		    bool normal=true, bool atSolution=false,
		    bool infeasible=false);

    /// Copy constructor
    DcCutGenerator (const DcCutGenerator &);

    /// Assignment operator
    DcCutGenerator & operator=(const DcCutGenerator& rhs);

    /// Destructor
    ~DcCutGenerator ();
    //@}

    /**@name Gets and sets */
    //@{
    /** Set the client model.

	In addition to setting the client model, refreshModel also calls
	the \c refreshSolver method of the CglCutGenerator object.
    */
    void refreshModel(DcModel * model);

    /// return name of generator
    inline const char * cutGeneratorName() const
	{
	    return generatorName_;
	}

    /** Set the cut generation interval

	Set the number of nodes evaluated between calls to the Cgl object's
	\p generateCuts routine.

	If \p value is positive, cuts will always be generated at the specified
	interval.
	If \p value is negative, cuts will initially be generated at the
	specified interval, but Dc may adjust the value depending on the
	success of cuts produced by this generator.

	A value of -100 disables the generator, while a value of -99 means
	just at root.
    */
    void setHowOften(int value) ;

    /// Get the cut generation interval.
    inline int howOften() const
	{ return whenCutGenerator_; }

    /// Get whether the cut generator should be called in the normal place
    inline bool normal() const
	{ return normal_; }
    /// Set whether the cut generator should be called in the normal place
    inline void setNormal(bool value)
	{ normal_=value; }
    /// Get whether the cut generator should be called when a solution is found
    inline bool atSolution() const
	{ return atSolution_; }
    /// Set whether the cut generator should be called when a solution is found
    inline void setAtSolution(bool value)
	{ atSolution_=value; }
    /** Get whether the cut generator should be called when the subproblem is
	found to be infeasible.
    */
    inline bool whenInfeasible() const
	{ return whenInfeasible_; }
    /** Set whether the cut generator should be called when the subproblem is
	found to be infeasible.
    */
    inline void setWhenInfeasible(bool value)
	{ whenInfeasible_=value; }
    /// Get the \c CglCutGenerator bound to this \c DcCutGenerator.
    inline CglCutGenerator * generator() const
	{ return generator_; }
    //@}

 private:
    /// The client model
    DcModel *model_;

    // The CglCutGenerator object
    CglCutGenerator * generator_;

    /** Number of nodes between calls to the CglCutGenerator::generateCuts
	routine.
    */
    int whenCutGenerator_;

    /// Name of generator
    char * generatorName_;

    /// Whether to call the generator in the normal place
    bool normal_;

    /// Whether to call the generator when a new solution is found
    bool atSolution_;

    /// Whether to call generator when a subproblem is found to be infeasible
    bool whenInfeasible_;
};

#endif
