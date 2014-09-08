//#############################################################################
// This file is modified from COIN-OR's ALPS library example file
// AbcParams.cpp
//#############################################################################

#include <AlpsParameterBase.h>
#include "DcParams.hpp"

using std::make_pair;

void
DcParams::createKeywordList() {

    // Create the list of keywords for parameter file reading
    //-------------------------------------------------------------------------
    // CharPar
    keys_.push_back(make_pair(std::string("Dc_cutDuringRampup"),
			      AlpsParameter(AlpsBoolPar, cutDuringRampup)));

    //-------------------------------------------------------------------------
    // BoolArrayPar

    //-------------------------------------------------------------------------
    // IntPar
    keys_.push_back(make_pair(std::string("Dc_statusInterval"),
			      AlpsParameter(AlpsIntPar, statusInterval)));
    keys_.push_back(make_pair(std::string("Dc_logLevel"),
			      AlpsParameter(AlpsIntPar, logLevel)));

    //-------------------------------------------------------------------------
    // DoublePar

    //-------------------------------------------------------------------------
    // StringPar

}

//#############################################################################

void
DcParams::setDefaultEntries() {
    //-------------------------------------------------------------------------
    // CharPar
    setEntry(cutDuringRampup, false);

    //-------------------------------------------------------------------------
    // IntPar
    setEntry(statusInterval, 50);
    setEntry(logLevel, 1);

    //-------------------------------------------------------------------------
    // DoublePar

    //-------------------------------------------------------------------------
    // StringPar

}
