//#############################################################################
// This file is modified from COIN-OR's ALPS library example file
// AbcSolution.cpp
//#############################################################################

#include <iomanip>
#include <iostream>
#include <set>

#include "DcModel.hpp"
#include "DcSolution.hpp"

//#############################################################################

void
DcSolution::print(std::ostream& os) const
{
    os <<std::setiosflags(std::ios::fixed|std::ios::showpoint)<<std::setw(14);

    os << "-------------------------" <<std::endl;
    for (int i = 0; i < size_; ++i) {
	if (fabs(value_[i]) > 1.0e-7) {
	    os << std::setw(6) << i << " " << value_[i] << std::endl;
	}
    }
    os << "-------------------------" <<std::endl;
    os <<std::resetiosflags(std::ios::fixed|std::ios::showpoint|std::ios::scientific);

}

/** The method that encodes the node into a encoded object. */
AlpsEncoded*
DcSolution::encode() const
{
  //  AlpsEncoded* encoded = new AlpsEncoded(typeid(*this).name());
  AlpsEncoded* encoded = new AlpsEncoded(AlpsKnowledgeTypeSolution);

  encoded->writeRep(size_);     // Base operand of `->' has non-pointer type
  encoded->writeRep(value_, size_);
  encoded->writeRep(objective_);

  return encoded;
}

//#############################################################################

/** The method that decodes the node from a encoded object. */
// Note: write and read sequence MUST same!
AlpsKnowledge*
DcSolution::decode(AlpsEncoded& encoded) const
{
  int s;
  double obj;
  double* val = 0;
  encoded.readRep(s);
  encoded.readRep(val, s);        // s must immediately before sol
  encoded.readRep(obj);
  DcSolution* sol = new DcSolution(s, val, obj);

  if (val != 0) {
      delete [] val;
      val = 0;
  }

  return sol;
}

//#############################################################################
