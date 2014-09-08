#ifndef DcTreeNode_hpp_
#define DcTreeNode_hpp_

//#############################################################################
// This file is modified from COIN-OR's ALPS library example file
// AbcTreeNode.h
//#############################################################################

#include <AlpsKnowledgeBroker.h>
#include <AlpsTreeNode.h>

#include "DcModel.hpp"
#include "DcNodeDesc.hpp"

//#############################################################################

enum DcVarStatus {
    DcVarFree,
    DcVarFixedToUB,
    DcVarFixedToLB
};

//#############################################################################

class DcTreeNode : public AlpsTreeNode {
 private:
    // NO: default constructor, copy constructor, assignment operator
    DcTreeNode(const DcTreeNode&);
    DcTreeNode& operator=(const DcTreeNode&);

 private:
    /** The index of the branching variable */
    int branchedOn_;

    /** The solution value (non-integral) of the branching variable. */
    double branchedOnVal_;

    /** Branching direction */
    int branchedDir_;

    /// Guessed satisfied Objective value
    double guessedObjectiveValue_;

    /// The number of objects unsatisfied at this node.
    int numberUnsatisfied_;

 public:
    DcTreeNode()
	:
	branchedOn_(-1),
	branchedOnVal_(ALPS_BND_MAX),
	branchedDir_(0),
	guessedObjectiveValue_(ALPS_OBJ_MAX),
	numberUnsatisfied_(0)
	{
	    desc_ = new DcNodeDesc(dynamic_cast<DcModel*>
                       (getKnowledgeBroker()->getModel()));
	}

    DcTreeNode(DcModel* m)
	:
	branchedOn_(-1),
	branchedOnVal_(ALPS_BND_MAX),
	branchedDir_(0),
	guessedObjectiveValue_(ALPS_OBJ_MAX),
	numberUnsatisfied_(0)
	{
	    desc_ = new DcNodeDesc(m);
	}

    DcTreeNode(DcNodeDesc*& desc)
	:
	branchedOn_(-1),
	branchedOnVal_(ALPS_BND_MAX),
	branchedDir_(0),
	guessedObjectiveValue_(ALPS_OBJ_MAX),
	numberUnsatisfied_(0)
	{
	    desc_ = desc;
	    desc = 0;
	    //At the time of registering node, that node hasn't set broker
	    //desc_->setModel(getKnowledgeBroker()->getDataPool()->getModel());
	}

    ~DcTreeNode()
	{
	}

    virtual AlpsTreeNode* createNewTreeNode(AlpsNodeDesc*& desc) const;

    /** Performing the bounding operation. */
    virtual int process(bool isRoot = false, bool rampUp = false);

    /** Select the branch variable.
	Return value:
	<ul>
	<li>  0: A branching object has been installed
	<li> -1: A monotone object was discovered
	<li> -2: An infeasible object was discovered
	</ul>
    */
    int chooseBranch (DcModel * model, bool& strongFound);

    /** Query/set the objective value (could be approximately or not exit)
	of the node. */
    ///@{
    inline double getObjValue() const { return quality_; }
    inline void setObjValue(const double objValue) { quality_ = objValue; }
    ///@}

    virtual std::vector< CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> >
	branch();

    /// Get the number of objects unsatisfied at this node.
    inline int numberUnsatisfied() const
	{ return numberUnsatisfied_; }

    /// Guessed objective value (for solution)
    inline double guessedObjectiveValue() const
	{ return guessedObjectiveValue_; }
    ///
    inline void setGuessedObjectiveValue(double value)
	{ guessedObjectiveValue_ = value; }
    ///
    void setBranchedOn(int b) { branchedOn_ = b; }
    ///
    void setBranchedDir(int d) { branchedDir_ = d; }
    ///
    void setBranchedOnValue(double b) { branchedOnVal_ = b; }
    ///
    int getBranchedOn() const { return branchedOn_; }
    ///
    int getBranchedDir() const { return branchedDir_; }
    ///
    double getBranchedOnValue() const { return branchedOnVal_; }
    ///
    virtual AlpsEncoded* encode() const;
    ///
    virtual AlpsKnowledge* decode(AlpsEncoded&) const;

};

#endif
