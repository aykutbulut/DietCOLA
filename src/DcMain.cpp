//#############################################################################
// This file is modified from COIN-OR's ALPS library example file
// AbcMain.cpp
//#############################################################################

#include "AlpsConfig.h"

#include <iostream>

#include "CoinError.hpp"
#include "CoinTime.hpp"
//#include "OsiSolverInterface.hpp"
//#include "OsiClpSolverInterface.hpp"
#include <ColaModel.hpp>

#include "CglFlowCover.hpp"
#include "CglGomory.hpp"
#include "CglProbing.hpp"
#include "CglKnapsackCover.hpp"
#include "CglOddHole.hpp"

#ifdef COIN_HAS_MPI
#  include "AlpsKnowledgeBrokerMPI.h"
#else
#  include "AlpsKnowledgeBrokerSerial.h"
#endif

#include "DcHeuristic.hpp"
#include "DcSolution.hpp"
#include "DcTreeNode.hpp"
#include "DcModel.hpp"

#include <ConicConstraints.hpp>
#include <numeric>
#include <iomanip>
//#############################################################################
//#############################################################################

int main(int argc, char* argv[])
{

    try{
	// Declare application parameter, model and knowledge broker
      OsiConicSolverInterface * solver = new ColaModel();
      //solver1->options()->set_int_option(LOG_LEVEL, 0);
      DcModel model(*solver);
      //solver1.messageHandler()->setLogLevel(0);

#ifdef COIN_HAS_MPI
	AlpsKnowledgeBrokerMPI broker(argc, argv, model);
#else
	AlpsKnowledgeBrokerSerial broker(argc, argv, model);
#endif
	int DcLogLevel = model.DcPar()->entry(DcParams::logLevel);
	model.messageHandler()->setLogLevel(DcLogLevel);

	// Register model, solution, and tree node
	broker.registerClass(AlpsKnowledgeTypeModel, new DcModel);
	broker.registerClass(AlpsKnowledgeTypeSolution, new DcSolution);
	broker.registerClass(AlpsKnowledgeTypeNode, new DcTreeNode(&model));

	// Formulate the root node
	// NOTE: root will be deleted by ALPS
	AlpsTreeNode* root = new DcTreeNode(&model);

	// Search for solutions from give root.
	broker.rootSearch(root);

	// REPORT FEASIBILITY
	// get solution
	const double * full_sol = model.bestSolution();
	const ConicConstraints * cc = dynamic_cast<ColaModel*>(model.solver())->get_conic_constraints();
	// we only care about conic constraints since the rest is
	// assured by CLP and ALPS
	//model.solver()->report_feasibility();
	int num_cones = cc->num_cones();
	double * par_sol;
	double lhs = 0.0;
	double lhs_real = 0.0;
	// std::cout << "Conic Constraints feasibility report" << std::endl;
	// std::cout << std::setw(5) << std::left << "Cone"
	// 	  << std::setw(15) << std::left << "lhs"
	// 	  << std::setw(15) << std::left << "lhs_real"
	// 	  << std::endl;
	// for (int i=0; i<num_cones; ++i) {
	//   int cone_size = cc->cone_size(i);
	//   const int * members = cc->cone_members(i);
	//   par_sol = new double[cone_size];
	//   for (int j=0; j<cone_size; ++j) {
	//     par_sol[j] = full_sol[members[j]];
	//   }
	//   if (cc->type(i)==QUAD) {
	//     lhs = par_sol[0]*par_sol[0]
	//       - std::inner_product(par_sol+1, par_sol+cone_size, par_sol+1, 0.0);
	//     lhs_real = par_sol[0]
	//       -sqrt(std::inner_product(par_sol+1, par_sol+cone_size, par_sol+1, 0.0));
	//   }
	//   else if (cc->type(i)==RQUAD) {
	//     lhs = 2.0*par_sol[0]*par_sol[1]
	//       - std::inner_product(par_sol+2, par_sol+cone_size, par_sol+2, 0.0);
	//     lhs_real = lhs;
	//   }
	//   std::cout << std::setw(5) << std::left << i
	// 	    << std::setw(15) << std::left << lhs
	// 	    << std::setw(15) << std::left << lhs_real
	// 	    << std::endl;
	//   delete[] par_sol;
	// }
    }
    catch(CoinError& er) {
	std::cerr << "ERROR:" << er.message() << std::endl
		  << " from function " << er.methodName() << std::endl
		  << " from class " << er.className() << std::endl;
    }
    catch(...) {
	std::cerr << "Something went wrong!" << std::endl;
    }
    return 0;
}
