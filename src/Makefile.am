lib_LTLIBRARIES=libDietcola.la

libDietcola_la_SOURCES=\
	DcBranchActual.hpp\
	DcBranchActual.cpp\
	DcBranchBase.hpp\
	DcBranchBase.cpp\
	DcCutGenerator.hpp\
	DcCutGenerator.cpp\
	DcHeuristic.hpp\
	DcHeuristic.cpp\
	DcMessage.hpp\
	DcMessage.cpp\
	DcModel.hpp\
	DcModel.cpp\
	DcParams.hpp\
	DcParams.cpp\
	DcSolution.hpp\
	DcSolution.cpp\
	DcTreeNode.hpp\
	DcTreeNode.cpp\
	DcNodeDesc.hpp

libDietcola_la_LIBADD=\
	$(COLA_LIBS)\
	$(ALPS_LIBS)\
	$(CGL_LIBS)\
	$(LIBS)

bin_PROGRAMS=dietcola

dietcola_SOURCES=DcMain.cpp

dietcola_LDADD=libDietcola.la\
	$(COLA_LIBS)\
	$(OSICONIC_LIBS)\
	$(ALPS_LIBS)\
	$(CGL_LIBS)\
	$(LIBS)

AM_CPPFLAGS=\
	-I../../include\
	$(COLA_CFLAGS)\
	$(OSICONIC_CFLAGS)\
	$(ALPS_CFLAGS)\
	$(CGL_CFLAGS)

# header files to be installed to $(includedir)
include_HEADERS=DcMessage.hpp\
	DcModel.hpp\
	DcBranchActual.hpp\
	DcBranchBase.hpp\
	DcCutGenerator.hpp\
	DcHeuristic.hpp\
	DcNodeDesc.hpp\
	DcParams.hpp\
	DcSolution.hpp\
	DcTreeNode.hpp

# this are in the disctribution but not in include
#noinst_HEADERS=
