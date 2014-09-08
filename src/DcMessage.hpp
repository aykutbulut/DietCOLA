#ifndef DcMessage_hpp_
#define DcMessage_hpp_

//#############################################################################
// This file is modified from COIN-OR's ALPS library example file
// AbcMessage.h
//#############################################################################

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

/** This deals with Dc messages (as against Clp messages etc).
    CoinMessageHandler.hpp is the general part of message handling.
    All it has are enum's for the various messages.
    DcMessage.cpp has text in various languages.

    It is trivial to use the .hpp and .cpp file as a basis for
    messages for other components.
 */

#include "CoinMessageHandler.hpp"

enum DC_Message
{
  DC_END_GOOD,
  DC_MAXNODES,
  DC_MAXTIME,
  DC_MAXSOLS,
  DC_SOLUTION,
  DC_END,
  DC_INFEAS,
  DC_STRONG,
  DC_SOLINDIVIDUAL,
  DC_INTEGERINCREMENT,
  DC_STATUS,
  DC_GAP,
  DC_ROUNDING,
  DC_ROOT,
  DC_GENERATOR,
  DC_BRANCH,
  DC_STRONGSOL,
  DC_NOINT,
  DC_VUB_PASS,
  DC_VUB_END,
  DC_NOTFEAS1,
  DC_NOTFEAS2,
  DC_NOTFEAS3,
  DC_CUTOFF_WARNING1,
  DC_CUTS,
  DC_BRANCHSOL,
  DC_DUMMY_END
};

class DcMessage : public CoinMessages {

public:

  /**@name Constructors etc */
  //@{
  /** Constructor */
  DcMessage(Language language=us_en);
  //@}

};

#endif
