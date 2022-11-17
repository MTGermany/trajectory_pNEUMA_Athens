#
# makefile of general c++ projects
# source: ~/versionedProjects/lib/templates/makefile

# Achtung! Auch ohne .h File muss man bei $OBJECTS immer auch das 
# File mit der Main-Methode dazunehmen!
# 

#compiler 
CC=g++ -Wall -O3

#output:
OUTNAME=extractTraj_pNEUMA 
OUTNAME2=extractTraj2


#include path for headers
# !!! ACHTUNG keine Comments hinter Befehle/Defs => obskurer Fehler!!

LIBDIR=~/versionedProjects/lib/trunk/
BINDIR=~/bin
#STATDIR=~/trafficSim/sources/statistics  # bis jan09, nun in lib
#LIBDIR=~/trafficSim/sources/lib # bis 2011


#objects/other programs for linking

OBJECTS_MAIN=extractTraj_pNEUMA.o
OBJECTS_MAIN2=extractTraj2.o
LIBOBJECTS=general.o InOut.o RandomUtils.o Statistics.o


###########################################################
# suffix regel: mache aus *.cpp ein *.o, und zwar fuer eingabedatei $<
.cpp.o:
	${CC} -I ${LIBDIR}  -c $<

############################################################

#wenn man nur h-file veraendert, sollte make alles neu compilieren
#irgendwie muss das gehen!


#o-files werden in "." als target verzeichnis geschrieben!
#dann wird alles von hier gelinkt. sehr gut so 
#(weil InOut.o ja zu diesem build gehoert)
general.o: ${LIBDIR}/general.cpp
	${CC} -c ${LIBDIR}/general.cpp -o general.o

InOut.o: ${LIBDIR}/InOut.cpp
	${CC} -c ${LIBDIR}/InOut.cpp -o InOut.o

Statistics.o: ${LIBDIR}/Statistics.cpp
	${CC} -c ${LIBDIR}/Statistics.cpp -o Statistics.o

Math.o: ${LIBDIR}/Math.cpp
	${CC} -c ${LIBDIR}/Math.cpp -o Math.o

RandomUtils.o: ${LIBDIR}/RandomUtils.cpp
	${CC} -c ${LIBDIR}/RandomUtils.cpp -o RandomUtils.o

##########################################################


#target micNew: abhaengigkeiten in Macro $(OBJECTS_MAIN) zusammengefasst 
#wichtig ist dass naechste befehlszeile mit TAB anfaengt:
extractTraj_pNEUMA: $(OBJECTS_MAIN)  $(LIBOBJECTS)
	${CC} -o ${BINDIR}/${OUTNAME} -O3 $(OBJECTS_MAIN) $(LIBOBJECTS) -lm

extractTraj2: $(OBJECTS_MAIN2)  $(LIBOBJECTS)
	${CC} -o ${BINDIR}/${OUTNAME2} -O3 $(OBJECTS_MAIN2) $(LIBOBJECTS) -lm


#
# Misc.
#
clean:
	rm -v *.o


#-I directory
#Mit der -I-Compileroption kann man beim Übersetzen eines Programmes die Liste der Verzeichnisse erweitern, in denen nach einer Datei gesucht wird.
# gcc -Iinc hello.c
#sucht nach stdio.h zuerst als inc/stdio.h, und erst dann als /usr/include/stdio.h.


#The -I option or the INCLUDE_PATH variable described below should always be used to set the include path for header files.
