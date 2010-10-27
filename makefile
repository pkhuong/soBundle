##############################################################################
################################# makefile ###################################
##############################################################################
#									     #
#   makefile of bundle, constrained version				     #
#									     #
#   'make clean' cleans up						     #
#   'make' or 'make bundle' builds the module				     #
#									     #
#                                VERSION 3.10				     #
#                	        18 - 10 - 2004				     #
#									     #
# 		               Implementation by:			     #
#									     #
#			       Antonio Frangioni			     #
#									     #
#   			   Operations Research Group			     #
#			  Dipartimento di Informatica			     #
#   			     Universita' di Pisa			     #
#									     #
##############################################################################

# module name
NAME = bundle

# debug switches
SW = -g -ggdb -O2 -fPIC -W -Wall

# production switches
#SW = -O -w

# libreries
LIB = -lpthread -lm

# compiler
CC = g++
AR = ar

# default target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

default: $(NAME)

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean:
	rm -f *.o *~ lib$(NAME).so lib$(NAME).a


# main module (linking phase) - - - - - - - - - - - - - - - - - - - - - - - -

OBJ =   Bundle.o \
	MinQuad.o \
	BMinQuad.o \
	QPBundle.o \
	bundle_solver.o \

$(NAME): $(OBJ)
	$(CC) -shared -o lib$(NAME).so $(OBJ) $(LIB) $(SW)
	$(AR) -r lib$(NAME).a $(OBJ)

# dependencies: every .o from its .C + every (recursively) included .h- - - -

Bundle.o: Bundle.C Bundle.h OPTtypes.h OPTvect.h
	$(CC) -c $*.C -o $@ $(SW)

MinQuad.o: MinQuad.C MinQuad.h OPTtypes.h OPTvect.h
	$(CC) -c $*.C -o $@ $(SW)

BMinQuad.o: BMinQuad.C BMinQuad.h MinQuad.h OPTvect.h OPTtypes.h
	$(CC) -c $*.C -o $@ $(SW)

QPBundle.o: QPBundle.C QPBundle.h Bundle.h BMinQuad.h MinQuad.h OPTtypes.h\
	OPTvect.h
	$(CC) -c $*.C -o $@ $(SW)

bundle_solver.o: bundle_solver.C QPBundle.h Bundle.h\
	 BMinQuad.h MinQuad.h OPTtypes.h OPTvect.h
	$(CC) -c $*.C -o $@ $(SW)

############################# End of makefile ################################
