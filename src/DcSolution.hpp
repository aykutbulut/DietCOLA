#ifndef DcSolution_hpp_
#define DcSolution_hpp_

//#############################################################################
// This file is modified from COIN-OR's ALPS library example file
// AbcSolution.h
//#############################################################################

#include <AlpsSolution.h>

#include "DcModel.hpp"

/** This class holds a MIP feasible primal solution. */
class DcSolution : public AlpsSolution {
 private:
    int size_;
    double* value_;
    double objective_;

 public:
    DcSolution()
	:
	size_(0),
	value_(0),
	objective_()
	{}
    DcSolution(const int s, const double* val, const double obj)
	:
	size_(s)
	{
	    if (size_ >= 0) {
		value_ = new double [size_];
		memcpy(value_, val, sizeof(double) * size_);
	    }
	}

    ~DcSolution() {
	if (value_ != 0) {
	    delete [] value_;
	    value_ = 0;
	}
    }

    /** Get the objective value value */
    double getObjValue() const { return objective_; }

    virtual double getQuality() const { return getObjValue(); }

    /** Get the size of the solution */
    int getSize() const { return size_; }

    /** Get the column solution */
    const double* getColSolution() const
	{ return value_; }

    /** Get item i in the solution vector */
    double getColSolution(int i) const { return value_[i]; }

    /** Print out the solution.*/
    virtual void print(std::ostream& os) const;

    /** The method that encodes the solution into a encoded object. */
    virtual AlpsEncoded* encode() const;

    /** The method that decodes the solution from a encoded object. */
    virtual AlpsKnowledge* decode(AlpsEncoded&) const;
};

#endif
