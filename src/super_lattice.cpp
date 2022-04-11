#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <math.h>
#include "opt.h"
using namespace std;

int super_lattice(int argc, char **argv, optpara opts){
        int nx = 1, ny = 1, nz = 1;

        if(opts.cellnum.length() > 0) {
          istringstream cell_iss(opts.cellnum.c_str());
          cell_iss >> nx >> ny >> nz;
          if(nx == 0 || ny == 0 || nz == 0) {
            cerr << "Cell number is wrong" << endl;
            return 1;
           }
         }
        cerr << "cell size " << nx << "x" << ny << "x" << nz << endl;
        cerr << "cell shape " << opts.shape << endl;
        cerr << "rand " << opts.seed << endl;

  string empty = "";    // initialize stringstream
  ostringstream os;

  srand(opts.seed);

  if(opts.atomcomp == "slr-SiGe")
  {
    double lz = 5.431, ly = lz, lx = lz/4.0;
    int numSi = 0, numGe = 0;
    os << "str.txt";
    ofstream fout_str(os.str().c_str());
    os.str(empty);
        fout_str << "Gear(5) RAND=10 " << opts.potential << " T=300.0 STEP=(0,1000100)" << endl;
        fout_str << "BLOCK FIXANGLE W=INF Q=INF LOGSTEP=10 PRINTVELOCITY dt=0.1fs scratch(100)" << endl << endl;
        fout_str << "Si2" << endl << endl;
        fout_str << lx << " 0 0" << opts.pbc_x << endl;
        fout_str << "0 " << ly << " 0" << opts.pbc_y << endl;
        fout_str << "0 0 " << lz << opts.pbc_z << endl << endl;
        fout_str << "Si 0.00 0.00 0.00" << endl << "Si 0.00 0.50 0.50" << endl;
        fout_str.close();
        numSi += 2;

    for(int i=1; i<nx; i++)
    {
        os << "layer.txt";
        ofstream fout_layer(os.str().c_str());
        os.str(empty);
        int n = rand()%2;
        string atm;
        if(!n || numSi<numGe) {atm = "Si"; numSi += 2;} else {atm = "Ge"; numGe += 2;}
        fout_layer << "Gear(5) RAND=10 " << opts.potential << " T=300.0 STEP=(0,1000100)" << endl;
        fout_layer << "BLOCK FIXANGLE W=INF Q=INF LOGSTEP=10 PRINTVELOCITY dt=0.1fs scratch(100)" << endl << endl;
        fout_layer << atm << "2" << endl << endl;
        fout_layer << lx << " 0 0" << opts.pbc_x << endl;
        fout_layer << "0 " << ly << " 0" << opts.pbc_y << endl;
        fout_layer << "0 0 " << lz << opts.pbc_z << endl << endl;
        if(i%4==0)      fout_layer << atm << " 0.00 0.00 0.00" << endl << atm << " 0.00 0.50 0.50" << endl;
        else if(i%4==1) fout_layer << atm << " 0.00 0.25 0.25" << endl << atm << " 0.00 0.75 0.75" << endl;
        else if(i%4==2) fout_layer << atm << " 0.00 0.50 0.00" << endl << atm << " 0.00 0.00 0.50" << endl;
        else if(i%4==3) fout_layer << atm << " 0.00 0.25 0.75" << endl << atm << " 0.00 0.75 0.25" << endl;
        fout_layer.close();

        os << "mdmake --stack str.txt -c layer.txt >temp.txt 2>/dev/null";
        system(os.str().c_str());
        os.str(empty);
        os << "cat temp.txt > str.txt";
        system(os.str().c_str());
        os.str(empty);
        os << "rm layer.txt temp.txt";
        system(os.str().c_str());
        os.str(empty);
    }

    numSi*=ny*nz;
    numGe*=ny*nz;
    int numall=numSi+numGe;
    double concGe=numGe/(double)numall;

    cerr << "atomnum " << numall << endl;
    cerr << "Sinum   " << numSi << endl;
    cerr << "Genum   " << numGe << endl;
    cerr << "Geconc  " << concGe << endl;

    os << "mdmake --convert mdl -c str.txt -n \"1 " << ny << " " << nz << "\" >str.mdl 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);
        os << "rm str.txt";
        system(os.str().c_str());
        os.str(empty);
    os << "mdmake --convert lmp -c str.mdl >str.lmp 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);
    os << "in.stable.SiGe";
    ofstream fout_in_stable(os.str().c_str());
    os.str(empty);
        fout_in_stable << "package      omp " << opts.omp << endl;
        fout_in_stable << "units        metal" << endl;
        fout_in_stable << "boundary     p p p" << endl;
        fout_in_stable << "atom_style   atomic" << endl;
        fout_in_stable << "atom_modify  sort 10000 1.0" << endl;
        fout_in_stable << "read_data    str.lmp" << endl;
        fout_in_stable << "pair_style   sw" << endl;
        fout_in_stable << "pair_coeff   * * SiGeSn.tomita16.eswge Si Ge" << endl;
        fout_in_stable << "neighbor     4.0 bin" << endl;
        fout_in_stable << "neigh_modify every 1 delay 0 check yes" << endl;
        fout_in_stable << "fix          1 all nve" << endl;
        fout_in_stable << "fix          2 all box/relax aniso 0.0 fixedpoint 0.0 0.0 0.0" << endl;
        fout_in_stable << "min_style    cg" << endl;
        fout_in_stable << "minimize     1e-25 1e-25 50000 100000" << endl;
        fout_in_stable << "dump         1 all custom 1 stable.final id type xs ys zs" << endl;
        fout_in_stable << "dump_modify  1 sort id" << endl;
        fout_in_stable << "run          0" << endl;
        fout_in_stable << "undump       1" << endl << endl;
        fout_in_stable.close();
    os << "mpirun -np " << opts.omp << " lmp_custom <in.stable.SiGe >/dev/null 2>&1";
    system(os.str().c_str());
    os.str(empty);
    os << "mdmake --convert mdl --lammps SiGe -c stable.final > stable.mdl 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);
  }
  else if(opts.atomcomp == "sl1-SiGe")
  {
    double lz = 5.431, ly = lz, lx = lz/4.0;
    int numSi = 0, numGe = 0;
    os << "str.txt";
    ofstream fout_str(os.str().c_str());
    os.str(empty);
        fout_str << "Gear(5) RAND=10 " << opts.potential << " T=300.0 STEP=(0,1000100)" << endl;
        fout_str << "BLOCK FIXANGLE W=INF Q=INF LOGSTEP=10 PRINTVELOCITY dt=0.1fs scratch(100)" << endl << endl;
        fout_str << "Si2" << endl << endl;
        fout_str << lx << " 0 0" << opts.pbc_x << endl;
        fout_str << "0 " << ly << " 0" << opts.pbc_y << endl;
        fout_str << "0 0 " << lz << opts.pbc_z << endl << endl;
        fout_str << "Si 0.00 0.00 0.00" << endl << "Si 0.00 0.50 0.50" << endl;
        fout_str.close();
        numSi += 2;

    for(int i=1; i<nx; i++)
    {
        os << "layer.txt";
        ofstream fout_layer(os.str().c_str());
        os.str(empty);
        string atm;
        if(i%2) {atm = "Ge"; numGe += 2;} else {atm = "Si"; numSi += 2;}
        fout_layer << "Gear(5) RAND=10 " << opts.potential << " T=300.0 STEP=(0,1000100)" << endl;
        fout_layer << "BLOCK FIXANGLE W=INF Q=INF LOGSTEP=10 PRINTVELOCITY dt=0.1fs scratch(100)" << endl << endl;
        fout_layer << atm << "2" << endl << endl;
        fout_layer << lx << " 0 0" << opts.pbc_x << endl;
        fout_layer << "0 " << ly << " 0" << opts.pbc_y << endl;
        fout_layer << "0 0 " << lz << opts.pbc_z << endl << endl;
        if(i%4==0)      fout_layer << atm << " 0.00 0.00 0.00" << endl << atm << " 0.00 0.50 0.50" << endl;
        else if(i%4==1) fout_layer << atm << " 0.00 0.25 0.25" << endl << atm << " 0.00 0.75 0.75" << endl;
        else if(i%4==2) fout_layer << atm << " 0.00 0.50 0.00" << endl << atm << " 0.00 0.00 0.50" << endl;
        else if(i%4==3) fout_layer << atm << " 0.00 0.25 0.75" << endl << atm << " 0.00 0.75 0.25" << endl;
        fout_layer.close();

        os << "mdmake --stack str.txt -c layer.txt >temp.txt 2>/dev/null";
        system(os.str().c_str());
        os.str(empty);
        os << "cat temp.txt > str.txt";
        system(os.str().c_str());
        os.str(empty);
        os << "rm layer.txt temp.txt";
        system(os.str().c_str());
        os.str(empty);
    }

    numSi*=ny*nz;
    numGe*=ny*nz;
    int numall=numSi+numGe;
    double concGe=numGe/(double)numall;

    cerr << "atomnum " << numall << endl;
    cerr << "Sinum   " << numSi << endl;
    cerr << "Genum   " << numGe << endl;
    cerr << "Geconc  " << concGe << endl;

    os << "mdmake --convert mdl -c str.txt -n \"1 " << ny << " " << nz << "\" >str.mdl 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);
        os << "rm str.txt";
        system(os.str().c_str());
        os.str(empty);
    os << "mdmake --convert lmp -c str.mdl >str.lmp 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);
    os << "in.stable.SiGe";
    ofstream fout_in_stable(os.str().c_str());
    os.str(empty);
        fout_in_stable << "package      omp " << opts.omp << endl;
        fout_in_stable << "units        metal" << endl;
        fout_in_stable << "boundary     p p p" << endl;
        fout_in_stable << "atom_style   atomic" << endl;
        fout_in_stable << "atom_modify  sort 10000 1.0" << endl;
        fout_in_stable << "read_data    str.lmp" << endl;
        fout_in_stable << "pair_style   sw" << endl;
        fout_in_stable << "pair_coeff   * * SiGeSn.tomita16.eswge Si Ge" << endl;
        fout_in_stable << "neighbor     4.0 bin" << endl;
        fout_in_stable << "neigh_modify every 1 delay 0 check yes" << endl;
        fout_in_stable << "fix          1 all nve" << endl;
        fout_in_stable << "fix          2 all box/relax aniso 0.0 fixedpoint 0.0 0.0 0.0" << endl;
        fout_in_stable << "min_style    cg" << endl;
        fout_in_stable << "minimize     1e-25 1e-25 50000 100000" << endl;
        fout_in_stable << "dump         1 all custom 1 stable.final id type xs ys zs" << endl;
        fout_in_stable << "dump_modify  1 sort id" << endl;
        fout_in_stable << "run          0" << endl;
        fout_in_stable << "undump       1" << endl << endl;
        fout_in_stable.close();
    os << "mpirun -np " << opts.omp << " lmp_custom <in.stable.SiGe >/dev/null 2>&1";
    system(os.str().c_str());
    os.str(empty);
    os << "mdmake --convert mdl --lammps SiGe -c stable.final > stable.mdl 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);
  }
  else if(opts.atomcomp == "sl2-SiGe")
  {
    double lz = 5.431, ly = lz, lx = lz/4.0;
    int numSi = 0, numGe = 0;
    os << "str.txt";
    ofstream fout_str(os.str().c_str());
    os.str(empty);
        fout_str << "Gear(5) RAND=10 " << opts.potential << " T=300.0 STEP=(0,1000100)" << endl;
        fout_str << "BLOCK FIXANGLE W=INF Q=INF LOGSTEP=10 PRINTVELOCITY dt=0.1fs scratch(100)" << endl << endl;
        fout_str << "Si2" << endl << endl;
        fout_str << lx << " 0 0" << opts.pbc_x << endl;
        fout_str << "0 " << ly << " 0" << opts.pbc_y << endl;
        fout_str << "0 0 " << lz << opts.pbc_z << endl << endl;
        fout_str << "Si 0.00 0.00 0.00" << endl << "Si 0.00 0.50 0.50" << endl;
        fout_str.close();
        numSi += 2;

    for(int i=1; i<nx; i++)
    {
        os << "layer.txt";
        ofstream fout_layer(os.str().c_str());
        os.str(empty);
        string atm;
        if(i%4==2 || i%4==3) {atm = "Ge"; numGe += 2;} else {atm = "Si"; numSi += 2;}
        fout_layer << "Gear(5) RAND=10 " << opts.potential << " T=300.0 STEP=(0,1000100)" << endl;
        fout_layer << "BLOCK FIXANGLE W=INF Q=INF LOGSTEP=10 PRINTVELOCITY dt=0.1fs scratch(100)" << endl << endl;
        fout_layer << atm << "2" << endl << endl;
        fout_layer << lx << " 0 0" << opts.pbc_x << endl;
        fout_layer << "0 " << ly << " 0" << opts.pbc_y << endl;
        fout_layer << "0 0 " << lz << opts.pbc_z << endl << endl;
        if(i%4==0)      fout_layer << atm << " 0.00 0.00 0.00" << endl << atm << " 0.00 0.50 0.50" << endl;
        else if(i%4==1) fout_layer << atm << " 0.00 0.25 0.25" << endl << atm << " 0.00 0.75 0.75" << endl;
        else if(i%4==2) fout_layer << atm << " 0.00 0.50 0.00" << endl << atm << " 0.00 0.00 0.50" << endl;
        else if(i%4==3) fout_layer << atm << " 0.00 0.25 0.75" << endl << atm << " 0.00 0.75 0.25" << endl;
        fout_layer.close();

        os << "mdmake --stack str.txt -c layer.txt >temp.txt 2>/dev/null";
        system(os.str().c_str());
        os.str(empty);
        os << "cat temp.txt > str.txt";
        system(os.str().c_str());
        os.str(empty);
        os << "rm layer.txt temp.txt";
        system(os.str().c_str());
        os.str(empty);
    }

    numSi*=ny*nz;
    numGe*=ny*nz;
    int numall=numSi+numGe;
    double concGe=numGe/(double)numall;

    cerr << "atomnum " << numall << endl;
    cerr << "Sinum   " << numSi << endl;
    cerr << "Genum   " << numGe << endl;
    cerr << "Geconc  " << concGe << endl;

    os << "mdmake --convert mdl -c str.txt -n \"1 " << ny << " " << nz << "\" >str.mdl 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);
        os << "rm str.txt";
        system(os.str().c_str());
        os.str(empty);
    os << "mdmake --convert lmp -c str.mdl >str.lmp 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);
    os << "in.stable.SiGe";
    ofstream fout_in_stable(os.str().c_str());
    os.str(empty);
        fout_in_stable << "package      omp " << opts.omp << endl;
        fout_in_stable << "units        metal" << endl;
        fout_in_stable << "boundary     p p p" << endl;
        fout_in_stable << "atom_style   atomic" << endl;
        fout_in_stable << "atom_modify  sort 10000 1.0" << endl;
        fout_in_stable << "read_data    str.lmp" << endl;
        fout_in_stable << "pair_style   sw" << endl;
        fout_in_stable << "pair_coeff   * * SiGeSn.tomita16.eswge Si Ge" << endl;
        fout_in_stable << "neighbor     4.0 bin" << endl;
        fout_in_stable << "neigh_modify every 1 delay 0 check yes" << endl;
        fout_in_stable << "fix          1 all nve" << endl;
        fout_in_stable << "fix          2 all box/relax aniso 0.0 fixedpoint 0.0 0.0 0.0" << endl;
        fout_in_stable << "min_style    cg" << endl;
        fout_in_stable << "minimize     1e-25 1e-25 50000 100000" << endl;
        fout_in_stable << "dump         1 all custom 1 stable.final id type xs ys zs" << endl;
        fout_in_stable << "dump_modify  1 sort id" << endl;
        fout_in_stable << "run          0" << endl;
        fout_in_stable << "undump       1" << endl << endl;
        fout_in_stable.close();
    os << "mpirun -np " << opts.omp << " lmp_custom <in.stable.SiGe >/dev/null 2>&1";
    system(os.str().c_str());
    os.str(empty);
    os << "mdmake --convert mdl --lammps SiGe -c stable.final > stable.mdl 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);
  }

        return 0;
}

