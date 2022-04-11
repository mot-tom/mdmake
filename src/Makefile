CC		= gcc
CXX		= g++
CFLAGS		= -O4 -Wall -I/usr/local/include
DEST		= ~/bin
LDFLAGS	= -L/usr/local/lib
OBJS		= main.o diamond.o diamond_layer.o oxi_layer.o oxi_amo.o diamond_oxide.o diamond_comp.o wurtz_layer.o tetra_layer.o tetra_comp.o orthodia.o orthodia2.o ortho_layer.o ortho2_layer.o zincblende.o wurtzite.o quartz.o tridymite.o rutile.o rocksalt.o corundum.o tetragonaldia.o confinementdia.o dichalcogenide.o fcc.o diamond-lmp.o super_lattice.o convert.o stack.o convert-org.o swap.o cluster.o rdf.o
PROGRAM	= mdmake

all:		$(PROGRAM)

$(PROGRAM):	$(OBJS)
		$(CXX) $(OBJS) $(LDFLAGS) -o $(PROGRAM)

clean:;	rm -f *.o *~ $(PROGRAM)

install:	$(PROGRAM)
		install -s $(PROGRAM) $(DEST)

