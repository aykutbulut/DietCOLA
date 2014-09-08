#ifndef DcHeuristic_hpp_
#define DcHeuristic_hpp_

//#############################################################################
// This file is modified from COIN-OR's ALPS library example file
// AbcHeuristic.h
//#############################################################################

#include <string>
#include <vector>
#include "CoinPackedMatrix.hpp"
#include "OsiCuts.hpp"

class OsiSolverInterface;
class DcModel;

//#############################################################################

/** Heuristic base class */
class DcHeuristic {
public:
  // Default Constructor
  DcHeuristic ();

  // Constructor with model - assumed before cuts
  DcHeuristic (DcModel & model);

  virtual ~DcHeuristic();

  /// update model (This is needed if cliques update matrix etc)
  virtual void setModel(DcModel * model);

  /// Clone
  virtual DcHeuristic * clone() const=0;

  /** returns 0 if no solution, 1 if valid solution
      with better objective value than one passed in
      Sets solution values if good, sets objective value
      This is called after cuts have been added - so can not add cuts
  */
  virtual int solution(double & objectiveValue,
		       double * newSolution)=0;

  /** returns 0 if no solution, 1 if valid solution, -1 if just
      returning an estimate of best possible solution
      with better objective value than one passed in
      Sets solution values if good, sets objective value (only if nonzero code)
      This is called at same time as cut generators - so can add cuts
      Default is do nothing
  */
  virtual int solution(double & objectiveValue,
		       double * newSolution,
		       OsiCuts & cs) {return 0;}

protected:

  /// Model
  DcModel * model_;
private:

  /// Illegal Assignment operator
  DcHeuristic & operator=(const DcHeuristic& rhs);

};

/** Rounding class
 */

class DcRounding : public DcHeuristic {
public:

  // Default Constructor
  DcRounding ();

  // Constructor with model - assumed before cuts
  DcRounding (DcModel & model);

  // Copy constructor
  DcRounding ( const DcRounding &);

  // Destructor
  ~DcRounding ();

  /// Clone
  virtual DcHeuristic * clone() const;

  /// update model (This is needed if cliques update matrix etc)
  virtual void setModel(DcModel * model);

  /** returns 0 if no solution, 1 if valid solution
      with better objective value than one passed in
      Sets solution values if good, sets objective value (only if good)
      This is called after cuts have been added - so can not add cuts
  */
  virtual int solution(double & objectiveValue,
		       double * newSolution);


  /// Set seed
  void setSeed(int value)
  { seed_ = value;}

protected:
  // Data

  // Original matrix by column
  CoinPackedMatrix matrix_;

  // Original matrix by
  CoinPackedMatrix matrixByRow_;

  // Seed for random stuff
  int seed_;

private:
  /// Illegal Assignment operator
  DcRounding & operator=(const DcRounding& rhs);
};


#endif

