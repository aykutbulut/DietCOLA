//#############################################################################
// This file is modified from COIN-OR's ALPS library example file
// AbcMessage.cpp
//#############################################################################

#include "DcMessage.hpp"
#include <cstring>

//#############################################################################

typedef struct {
    DC_Message internalNumber;
    int externalNumber;              // or continuation
    char detail;
    const char * message;
} Dc_message;

//#############################################################################

static Dc_message us_english[]=
{
    {DC_END_GOOD,1,1,"Search completed - best objective %g, took %d iterations and %d nodes"},
    {DC_MAXNODES,3,1,"Exiting on maximum nodes"},
    {DC_MAXTIME,20,1,"Exiting on maximum time"},
    {DC_MAXSOLS,19,1,"Exiting on maximum solutions"},
    {DC_SOLUTION,4,1,"Integer solution of %g found after %d iterations and %d nodes"},
    {DC_END,5,1,"Partial search took %d iterations and %d nodes"},
    {DC_INFEAS,6,1,"The LP relaxation is infeasible or too expensive"},
    {DC_STRONG,7,3,"Strong branching on %d (%d), down %g (%d) up %g (%d) value %d"},
    {DC_SOLINDIVIDUAL,8,2,"%d has value %g"},
    {DC_INTEGERINCREMENT,9,1,"Objective coefficients multiple of %g"},
    {DC_STATUS,10,1,"Process[%d]: after %d nodes, %d on tree, %g best solution, best possible %g"},
    {DC_GAP,11,1,"Exiting as integer gap of %g less than %g"},
    {DC_ROUNDING,12,1,"Integer solution of %g found by rounding after %d iterations and %d nodes"},
    {DC_ROOT,13,3,"At root node, %d cuts changed objective from %g to %g in %d passes"},
    {DC_GENERATOR,14,2,"Cut generator %d (%s) - %d row cuts (%d active), %d column cuts - new frequency is %d"},
    {DC_BRANCH,15,3,"Node %d Obj %g Unsat %d depth %d"},
    {DC_STRONGSOL,16,1,"Integer solution of %g found by strong branching after %d iterations and %d nodes"},
    {DC_NOINT,3007,0,"No integer variables - nothing to do"},
    {DC_VUB_PASS,17,1,"%d solved, %d variables fixed, %d tightened"},
    {DC_VUB_END,18,1,"After tightenVubs, %d variables fixed, %d tightened"},
    {DC_NOTFEAS1,21,2,"On closer inspection node is infeasible"},
    {DC_NOTFEAS2,22,2,"On closer inspection objective value of %g above cutoff of %g"},
    {DC_NOTFEAS3,23,2,"Allowing solution, even though largest row infeasibility is %g"},
    {DC_CUTOFF_WARNING1,23,1,"Cutoff set to %g - equivalent to best solution of %g"},
    {DC_CUTS,24,1, "At node %d, %d cuts changed objective from %g to %g in %d passes"},
    {DC_BRANCHSOL,25,1,"Integer solution of %g found by branching after %d iterations and %d nodes"},
    {DC_DUMMY_END, 999999, 0, ""}
};

//#############################################################################

/* Constructor */
DcMessage::DcMessage(Language language)
    :
    CoinMessages(sizeof(us_english) / sizeof(Dc_message))
{
    language_ = language;
    strcpy(source_, "Dc");
    Dc_message * message = us_english;

    while (message->internalNumber != DC_DUMMY_END) {
	CoinOneMessage oneMessage(message->externalNumber, message->detail,
				  message->message);
	addMessage(message->internalNumber, oneMessage);
	message++;
    }

    // now override any language ones

    switch (language) {

    default:
	message = NULL;
	break;
    }

    // replace if any found
    if (message) {
	while (message->internalNumber != DC_DUMMY_END) {
	    replaceMessage(message->internalNumber, message->message);
	    message++;
	}
    }
}
