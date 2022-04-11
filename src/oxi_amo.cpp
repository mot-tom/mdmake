#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <math.h>
#include "opt.h"
using namespace std;

int oxi_amo(int argc, char **argv, optpara opts){
        int nx = 4, ny = 4, nz = 4;

        if(opts.cellnum.length() > 0) {
          istringstream cell_iss(opts.cellnum.c_str());
          cell_iss >> nx >> ny >> nz;
          if(nx == 0 || ny == 0 || nz == 0) {
            cerr << "Cell number is wrong" << endl;
            return 1;
           }}

  string empty = "";    // initialize stringstream
  ostringstream os;

  if(opts.atomcomp == "am-Si")
  {
    os << "mdmake -c Si -n \"" << nx << " " << ny << " " << nz << "\" >str.txt 2>/dev/null";
    system(os.str().c_str());
    os.str(empty); 

    os << "mdmake --convert lmp -c str.txt >aneal.lmp 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);

    os << "in.aneal.Si";
    ofstream fout_in_stable(os.str().c_str());
    os.str(empty);

        fout_in_stable << "package      omp " << opts.omp << endl;
        fout_in_stable << "units        metal" << endl;
        fout_in_stable << "boundary     p p p" << endl;
        fout_in_stable << "atom_style   atomic" << endl;
        fout_in_stable << "atom_modify  sort 10000 1.0" << endl;
        fout_in_stable << "read_data    aneal.lmp" << endl;
        fout_in_stable << "pair_style   sw/omp" << endl;
        fout_in_stable << "pair_coeff   * * SiGeSn.tomita16.eswge Si" << endl;
        fout_in_stable << "neighbor     4.0 bin" << endl;
        fout_in_stable << "neigh_modify every 1 delay 0 check yes" << endl;
        fout_in_stable << "velocity     all create 9000 458127641 mom yes rot yes dist gaussian" << endl;
        fout_in_stable << "timestep     0.0001" << endl;//1fs
        fout_in_stable << "fix          1 all nvt temp 9000 1000 0.05" << endl;
        fout_in_stable << "thermo_style custom step temp ke pe etotal press vol density" << endl;
        fout_in_stable << "thermo       1000" << endl;
        fout_in_stable << "run          100000" << endl;
        fout_in_stable << "unfix        1" << endl;
        fout_in_stable << "velocity     all scale 1000" << endl;
        fout_in_stable << "fix          1 all npt temp 1000 300 0.05 x 0.0 0.0 0.5 y 0.0 0.0 0.5 z 0.0 0.0 0.5" << endl;
        fout_in_stable << "run          100000" << endl;
        fout_in_stable << "unfix        1" << endl;
        fout_in_stable << "fix          1 all nve" << endl;
        fout_in_stable << "fix          2 all box/relax aniso 0.0 fixedpoint 0.0 0.0 0.0" << endl;
        fout_in_stable << "min_style    cg" << endl;
        fout_in_stable << "minimize     1e-25 1e-25 50000 100000" << endl;
        fout_in_stable << "dump         1 all custom 1 aneal.final id type xs ys zs" << endl;
        fout_in_stable << "dump_modify  1 sort id" << endl;
        fout_in_stable << "run          0" << endl;
        fout_in_stable << "undump       1" << endl << endl;
        fout_in_stable.close();

    os << "mpirun -np " << opts.omp << " lmp_custom <in.aneal.Si";
    system(os.str().c_str());
    os.str(empty);

    os << "mdmake --convert mdl --lammps Si -c aneal.final >stable.txt 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);
  }
  else if(opts.atomcomp == "am-Ge")
  {
    os << "mdmake -c Ge -n \"" << nx << " " << ny << " " << nz << "\" >str.txt 2>/dev/null";
    system(os.str().c_str());
    os.str(empty); 

    os << "mdmake --convert lmp -c str.txt >aneal.lmp 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);

    os << "in.aneal.Ge";
    ofstream fout_in_stable(os.str().c_str());
    os.str(empty);

        fout_in_stable << "package      omp " << opts.omp << endl;
        fout_in_stable << "units        metal" << endl;
        fout_in_stable << "boundary     p p p" << endl;
        fout_in_stable << "atom_style   atomic" << endl;
        fout_in_stable << "atom_modify  sort 10000 1.0" << endl;
        fout_in_stable << "read_data    aneal.lmp" << endl;
        fout_in_stable << "pair_style   sw/omp" << endl;
        fout_in_stable << "pair_coeff   * * SiGeSn.tomita16.eswge Ge" << endl;
        fout_in_stable << "neighbor     4.0 bin" << endl;
        fout_in_stable << "neigh_modify every 1 delay 0 check yes" << endl;
        fout_in_stable << "velocity     all create 9000 458127641 mom yes rot yes dist gaussian" << endl;
        fout_in_stable << "timestep     0.0001" << endl;//1fs
        fout_in_stable << "fix          1 all nvt temp 9000 1000 0.05" << endl;
        fout_in_stable << "thermo_style custom step temp ke pe etotal press vol density" << endl;
        fout_in_stable << "thermo       1000" << endl;
        fout_in_stable << "run          100000" << endl;
        fout_in_stable << "unfix        1" << endl;
        fout_in_stable << "velocity     all scale 1000" << endl;
        fout_in_stable << "fix          1 all npt temp 1000 300 0.05 x 0.0 0.0 0.5 y 0.0 0.0 0.5 z 0.0 0.0 0.5" << endl;
        fout_in_stable << "run          100000" << endl;
        fout_in_stable << "unfix        1" << endl;
        fout_in_stable << "fix          1 all nve" << endl;
        fout_in_stable << "fix          2 all box/relax aniso 0.0 fixedpoint 0.0 0.0 0.0" << endl;
        fout_in_stable << "min_style    cg" << endl;
        fout_in_stable << "minimize     1e-25 1e-25 50000 100000" << endl;
        fout_in_stable << "dump         1 all custom 1 aneal.final id type xs ys zs" << endl;
        fout_in_stable << "dump_modify  1 sort id" << endl;
        fout_in_stable << "run          0" << endl;
        fout_in_stable << "undump       1" << endl << endl;
        fout_in_stable.close();

    os << "mpirun -np " << opts.omp << " lmp_custom <in.aneal.Ge";
    system(os.str().c_str());
    os.str(empty);

    os << "mdmake --convert mdl --lammps Ge -c aneal.final >stable.txt 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);
  }
  else if(opts.atomcomp == "am-ru-TiO2")
  {
    os << "mdmake -c ru-TiO2 -n \"" << nx << " " << ny << " " << nz*2 << "\" >str.txt 2>/dev/null";
    system(os.str().c_str());
    os.str(empty); 

    os << "mdmake --convert lmpCIM -c str.txt >aneal.lmp 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);

    os << "in.aneal.TiO";
    ofstream fout_in_stable(os.str().c_str());
    os.str(empty);

        fout_in_stable << "package      omp " << opts.omp << endl;
        fout_in_stable << "units        metal" << endl;
        fout_in_stable << "boundary     p p p" << endl;
        fout_in_stable << "atom_style   charge" << endl;
        fout_in_stable << "atom_modify  sort 10000 1.0" << endl;
        fout_in_stable << "read_data    aneal.lmp" << endl;
        fout_in_stable << "pair_style   born/coul/long/omp 10.0" << endl;
        fout_in_stable << "kspace_style ewald 1e-5" << endl;
        fout_in_stable << "pair_coeff   1 1 0.006935 0.160 2.470 0.000 0 10.0" << endl;
        fout_in_stable << "pair_coeff   1 2 0.007152 0.165 2.861 0.000 0 10.0" << endl;
        fout_in_stable << "pair_coeff   2 2 0.007369 0.170 3.252 4.144 0 10.0" << endl;
        fout_in_stable << "neighbor     4.0 bin" << endl;
        fout_in_stable << "neigh_modify every 1 delay 0 check yes" << endl;
        fout_in_stable << "velocity     all create 3000 458127641 mom yes rot yes dist gaussian" << endl;
        fout_in_stable << "timestep     0.0001" << endl;//1fs
        fout_in_stable << "fix          1 all nvt temp 30000 1000 0.05" << endl;
        fout_in_stable << "variable     Odensity equal atoms*2.0/3.0/vol" << endl;
        fout_in_stable << "thermo_style custom step temp ke pe etotal press vol density v_Odensity" << endl;
        fout_in_stable << "thermo       1000" << endl;
        fout_in_stable << "run          100000" << endl;
        fout_in_stable << "unfix        1" << endl;
        fout_in_stable << "velocity     all scale 1000" << endl;
        fout_in_stable << "fix          1 all npt temp 1000 300 0.05 x 0.0 0.0 0.5 y 0.0 0.0 0.5 z 0.0 0.0 0.5" << endl;
        fout_in_stable << "run          100000" << endl;
        fout_in_stable << "unfix        1" << endl;
        fout_in_stable << "fix          1 all nve" << endl;
        fout_in_stable << "fix          2 all box/relax aniso 0.0 fixedpoint 0.0 0.0 0.0" << endl;
        fout_in_stable << "min_style    cg" << endl;
        fout_in_stable << "minimize     1e-25 1e-25 50000 100000" << endl;
        fout_in_stable << "dump         1 all custom 1 aneal.final id type xs ys zs" << endl;
        fout_in_stable << "dump_modify  1 sort id" << endl;
        fout_in_stable << "run          0" << endl;
        fout_in_stable << "undump       1" << endl << endl;
        fout_in_stable.close();

    os << "mpirun -np " << opts.omp << " lmp_custom <in.aneal.TiO";
    system(os.str().c_str());
    os.str(empty);

    os << "mdmake --convert mdl --lammps TiO2 -c aneal.final >stable.txt 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);
  }
  else if(opts.atomcomp == "am-a-Al2O3")
  {
    os << "mdmake -c a-Al2O3 -n \"" << nx << " " << ny << " " << nz << "\" >str.txt 2>/dev/null";
    system(os.str().c_str());
    os.str(empty); 

    os << "mdmake --convert lmpCIM -c str.txt >aneal.lmp 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);

    os << "in.aneal.AlO";
    ofstream fout_in_stable(os.str().c_str());
    os.str(empty);

        fout_in_stable << "package      omp " << opts.omp << endl;
        fout_in_stable << "units        metal" << endl;
        fout_in_stable << "boundary     p p p" << endl;
        fout_in_stable << "atom_style   charge" << endl;
        fout_in_stable << "atom_modify  sort 10000 1.0" << endl;
        fout_in_stable << "read_data    aneal.lmp" << endl;
        fout_in_stable << "pair_style   born/coul/long/omp 10.0" << endl;
        fout_in_stable << "kspace_style ewald 1e-5" << endl;
        fout_in_stable << "pair_coeff   1 1 0.006935 0.160 2.128 0.000 0 10.0" << endl;
        fout_in_stable << "pair_coeff   1 2 0.007152 0.165 2.690 0.000 0 10.0" << endl;
        fout_in_stable << "pair_coeff   2 2 0.007369 0.170 3.252 4.144 0 10.0" << endl;
        fout_in_stable << "neighbor     4.0 bin" << endl;
        fout_in_stable << "neigh_modify every 1 delay 0 check yes" << endl;
        fout_in_stable << "velocity     all create 3000 458127641 mom yes rot yes dist gaussian" << endl;
        fout_in_stable << "timestep     0.0001" << endl;//1fs
        fout_in_stable << "fix          1 all nvt temp 24000 1000 0.05" << endl;
        fout_in_stable << "variable     Odensity equal atoms*2.0/3.0/vol" << endl;
        fout_in_stable << "thermo_style custom step temp ke pe etotal press vol density v_Odensity" << endl;
        fout_in_stable << "thermo       1000" << endl;
        fout_in_stable << "run          100000" << endl;
        fout_in_stable << "unfix        1" << endl;
        fout_in_stable << "velocity     all scale 1000" << endl;
        fout_in_stable << "fix          1 all npt temp 1000 300 0.05 x 0.0 0.0 0.5 y 0.0 0.0 0.5 z 0.0 0.0 0.5" << endl;
        fout_in_stable << "run          100000" << endl;
        fout_in_stable << "unfix        1" << endl;
        fout_in_stable << "fix          1 all nve" << endl;
        fout_in_stable << "fix          2 all box/relax aniso 0.0 fixedpoint 0.0 0.0 0.0" << endl;
        fout_in_stable << "min_style    cg" << endl;
        fout_in_stable << "minimize     1e-25 1e-25 50000 100000" << endl;
        fout_in_stable << "dump         1 all custom 1 aneal.final id type xs ys zs" << endl;
        fout_in_stable << "dump_modify  1 sort id" << endl;
        fout_in_stable << "run          0" << endl;
        fout_in_stable << "undump       1" << endl << endl;
        fout_in_stable.close();

    os << "mpirun -np " << opts.omp << " lmp_custom <in.aneal.AlO";
    system(os.str().c_str());
    os.str(empty);

    os << "mdmake --convert mdl --lammps Al2O3 -c aneal.final >stable.txt 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);
  }
  else if(opts.atomcomp == "am-MgO")
  {
    os << "mdmake -c MgO -n \"" << nx << " " << ny << " " << nz << "\" >str.txt 2>/dev/null";
    system(os.str().c_str());
    os.str(empty); 

    os << "mdmake --convert lmpCIM -c str.txt >aneal.lmp 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);

    os << "in.aneal.MgO";
    ofstream fout_in_stable(os.str().c_str());
    os.str(empty);

        fout_in_stable << "package      omp " << opts.omp << endl;
        fout_in_stable << "units        metal" << endl;
        fout_in_stable << "boundary     p p p" << endl;
        fout_in_stable << "atom_style   charge" << endl;
        fout_in_stable << "atom_modify  sort 10000 1.0" << endl;
        fout_in_stable << "read_data    aneal.lmp" << endl;
        fout_in_stable << "pair_style   born/coul/long/omp 10.0" << endl;
        fout_in_stable << "kspace_style ewald 1e-5" << endl;
        fout_in_stable << "pair_coeff   1 1 0.006935 0.160 2.322 0.041 0 10.0" << endl;
        fout_in_stable << "pair_coeff   1 2 0.007152 0.165 2.787 0.414 0 10.0" << endl;
        fout_in_stable << "pair_coeff   2 2 0.007369 0.170 3.252 4.144 0 10.0" << endl;
        fout_in_stable << "neighbor     4.0 bin" << endl;
        fout_in_stable << "neigh_modify every 1 delay 0 check yes" << endl;
        fout_in_stable << "velocity     all create 3000 458127641 mom yes rot yes dist gaussian" << endl;
        fout_in_stable << "timestep     0.0001" << endl;//1fs
        fout_in_stable << "fix          1 all nvt temp 31000 1000 0.05" << endl;
        fout_in_stable << "variable     Odensity equal atoms*2.0/3.0/vol" << endl;
        fout_in_stable << "thermo_style custom step temp ke pe etotal press vol density v_Odensity" << endl;
        fout_in_stable << "thermo       1000" << endl;
        fout_in_stable << "run          100000" << endl;
        fout_in_stable << "unfix        1" << endl;
        fout_in_stable << "velocity     all scale 1000" << endl;
        fout_in_stable << "fix          1 all npt temp 1000 300 0.05 x 0.0 0.0 0.5 y 0.0 0.0 0.5 z 0.0 0.0 0.5" << endl;
        fout_in_stable << "run          100000" << endl;
        fout_in_stable << "unfix        1" << endl;
        fout_in_stable << "fix          1 all nve" << endl;
        fout_in_stable << "fix          2 all box/relax aniso 0.0 fixedpoint 0.0 0.0 0.0" << endl;
        fout_in_stable << "min_style    cg" << endl;
        fout_in_stable << "minimize     1e-25 1e-25 50000 100000" << endl;
        fout_in_stable << "dump         1 all custom 1 aneal.final id type xs ys zs" << endl;
        fout_in_stable << "dump_modify  1 sort id" << endl;
        fout_in_stable << "run          0" << endl;
        fout_in_stable << "undump       1" << endl << endl;
        fout_in_stable.close();

    os << "mpirun -np " << opts.omp << " lmp_custom <in.aneal.MgO";
    system(os.str().c_str());
    os.str(empty);

    os << "mdmake --convert mdl --lammps MgO -c aneal.final >stable.txt 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);
  }
  else if(opts.atomcomp == "am-SrO")
  {
    os << "mdmake -c SrO -n \"" << nx << " " << ny << " " << nz << "\" >str.txt 2>/dev/null";
    system(os.str().c_str());
    os.str(empty); 

    os << "mdmake --convert lmpCIM -c str.txt >aneal.lmp 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);

    os << "in.aneal.SrO";
    ofstream fout_in_stable(os.str().c_str());
    os.str(empty);

        fout_in_stable << "package      omp " << opts.omp << endl;
        fout_in_stable << "units        metal" << endl;
        fout_in_stable << "boundary     p p p" << endl;
        fout_in_stable << "atom_style   charge" << endl;
        fout_in_stable << "atom_modify  sort 10000 1.0" << endl;
        fout_in_stable << "read_data    aneal.lmp" << endl;
        fout_in_stable << "pair_style   born/coul/long/omp 10.0" << endl;
        fout_in_stable << "kspace_style ewald 1e-5" << endl;
        fout_in_stable << "pair_coeff   1 1 0.006935 0.160 3.264 2.331 0 10.0" << endl;
        fout_in_stable << "pair_coeff   1 2 0.007152 0.165 3.258 3.108 0 10.0" << endl;
        fout_in_stable << "pair_coeff   2 2 0.007369 0.170 3.252 4.144 0 10.0" << endl;
        fout_in_stable << "neighbor     4.0 bin" << endl;
        fout_in_stable << "neigh_modify every 1 delay 0 check yes" << endl;
        fout_in_stable << "velocity     all create 3000 458127641 mom yes rot yes dist gaussian" << endl;
        fout_in_stable << "timestep     0.0001" << endl;//1fs
        fout_in_stable << "fix          1 all nvt temp 34000 1000 0.05" << endl;
        fout_in_stable << "variable     Odensity equal atoms*2.0/3.0/vol" << endl;
        fout_in_stable << "thermo_style custom step temp ke pe etotal press vol density v_Odensity" << endl;
        fout_in_stable << "thermo       1000" << endl;
        fout_in_stable << "run          100000" << endl;
        fout_in_stable << "unfix        1" << endl;
        fout_in_stable << "velocity     all scale 1000" << endl;
        fout_in_stable << "fix          1 all npt temp 1000 300 0.05 x 0.0 0.0 0.5 y 0.0 0.0 0.5 z 0.0 0.0 0.5" << endl;
        fout_in_stable << "run          100000" << endl;
        fout_in_stable << "unfix        1" << endl;
        fout_in_stable << "fix          1 all nve" << endl;
        fout_in_stable << "fix          2 all box/relax aniso 0.0 fixedpoint 0.0 0.0 0.0" << endl;
        fout_in_stable << "min_style    cg" << endl;
        fout_in_stable << "minimize     1e-25 1e-25 50000 100000" << endl;
        fout_in_stable << "dump         1 all custom 1 aneal.final id type xs ys zs" << endl;
        fout_in_stable << "dump_modify  1 sort id" << endl;
        fout_in_stable << "run          0" << endl;
        fout_in_stable << "undump       1" << endl << endl;
        fout_in_stable.close();

    os << "mpirun -np " << opts.omp << " lmp_custom <in.aneal.SrO";
    system(os.str().c_str());
    os.str(empty);

    os << "mdmake --convert mdl --lammps SrO -c aneal.final >stable.txt 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);
  }
  else if(opts.atomcomp == "am-aq-SiO2")
  {
    os << "mdmake -c aq-SiO2 -n \"" << nx << " " << ny << " " << nz << "\" >str.txt 2>/dev/null";
    system(os.str().c_str());
    os.str(empty); 

    os << "mdmake --convert lmp -c str.txt >aneal.lmp 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);

    os << "in.aneal.SiO";
    ofstream fout_in_stable(os.str().c_str());
    os.str(empty);

        fout_in_stable << "package      omp " << opts.omp << endl;
        fout_in_stable << "units        metal" << endl;
        fout_in_stable << "boundary     p p p" << endl;
        fout_in_stable << "atom_style   atomic" << endl;
        fout_in_stable << "atom_modify  sort 10000 1.0" << endl;
        fout_in_stable << "read_data    aneal.lmp" << endl;
        fout_in_stable << "pair_style   sw/wfnho" << endl;
        fout_in_stable << "pair_coeff   * * SiO.sw_wfnho Si O" << endl;
        fout_in_stable << "neighbor     4.0 bin" << endl;
        fout_in_stable << "neigh_modify every 1 delay 0 check yes" << endl;
        fout_in_stable << "velocity     all create 12000 458127641 mom yes rot yes dist gaussian" << endl;
        fout_in_stable << "timestep     0.0001" << endl;//1fs
        fout_in_stable << "fix          1 all nvt temp 12000 1000 0.05" << endl;
        fout_in_stable << "variable     Odensity equal atoms*2.0/3.0/vol" << endl;
        fout_in_stable << "thermo_style custom step temp ke pe etotal press vol density v_Odensity" << endl;
        fout_in_stable << "thermo       1000" << endl;
        fout_in_stable << "run          100000" << endl;
        fout_in_stable << "unfix        1" << endl;
        fout_in_stable << "velocity     all scale 1000" << endl;
        fout_in_stable << "fix          1 all npt temp 1000 300 0.05 x 0.0 0.0 0.5 y 0.0 0.0 0.5 z 0.0 0.0 0.5" << endl;
        fout_in_stable << "run          100000" << endl;
        fout_in_stable << "unfix        1" << endl;
        fout_in_stable << "fix          1 all nve" << endl;
        fout_in_stable << "fix          2 all box/relax aniso 0.0 fixedpoint 0.0 0.0 0.0" << endl;
        fout_in_stable << "min_style    cg" << endl;
        fout_in_stable << "minimize     1e-25 1e-25 50000 100000" << endl;
        fout_in_stable << "dump         1 all custom 1 aneal.final id type xs ys zs" << endl;
        fout_in_stable << "dump_modify  1 sort id" << endl;
        fout_in_stable << "run          0" << endl;
        fout_in_stable << "undump       1" << endl << endl;
        fout_in_stable.close();

    os << "mpirun -np " << opts.omp << " lmp_custom <in.aneal.SiO";
    system(os.str().c_str());
    os.str(empty);

    os << "mdmake --convert mdl --lammps SiO2 -c aneal.final >stable.txt 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);
  }
  else if(opts.atomcomp == "am-bq-SiO2")
  {
    os << "mdmake -c bq-SiO2 -n \"" << nx << " " << ny << " " << nz << "\" >str.txt 2>/dev/null";
    system(os.str().c_str());
    os.str(empty); 

    os << "mdmake --convert lmp -c str.txt >aneal.lmp 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);

    os << "in.aneal.SiO";
    ofstream fout_in_stable(os.str().c_str());
    os.str(empty);

        fout_in_stable << "package      omp " << opts.omp << endl;
        fout_in_stable << "units        metal" << endl;
        fout_in_stable << "boundary     p p p" << endl;
        fout_in_stable << "atom_style   atomic" << endl;
        fout_in_stable << "atom_modify  sort 10000 1.0" << endl;
        fout_in_stable << "read_data    aneal.lmp" << endl;
        fout_in_stable << "pair_style   sw/wfnho" << endl;
        fout_in_stable << "pair_coeff   * * SiO.sw_wfnho Si O" << endl;
        fout_in_stable << "neighbor     4.0 bin" << endl;
        fout_in_stable << "neigh_modify every 1 delay 0 check yes" << endl;
        fout_in_stable << "velocity     all create 12000 458127641 mom yes rot yes dist gaussian" << endl;
        fout_in_stable << "timestep     0.0001" << endl;//1fs
        fout_in_stable << "fix          1 all nvt temp 12000 1000 0.05" << endl;
        fout_in_stable << "variable     Odensity equal atoms*2.0/3.0/vol" << endl;
        fout_in_stable << "thermo_style custom step temp ke pe etotal press vol density v_Odensity" << endl;
        fout_in_stable << "thermo       1000" << endl;
        fout_in_stable << "run          100000" << endl;
        fout_in_stable << "unfix        1" << endl;
        fout_in_stable << "velocity     all scale 1000" << endl;
        fout_in_stable << "fix          1 all npt temp 1000 300 0.05 x 0.0 0.0 0.5 y 0.0 0.0 0.5 z 0.0 0.0 0.5" << endl;
        fout_in_stable << "run          100000" << endl;
        fout_in_stable << "unfix        1" << endl;
        fout_in_stable << "fix          1 all nve" << endl;
        fout_in_stable << "fix          2 all box/relax aniso 0.0 fixedpoint 0.0 0.0 0.0" << endl;
        fout_in_stable << "min_style    cg" << endl;
        fout_in_stable << "minimize     1e-25 1e-25 50000 100000" << endl;
        fout_in_stable << "dump         1 all custom 1 aneal.final id type xs ys zs" << endl;
        fout_in_stable << "dump_modify  1 sort id" << endl;
        fout_in_stable << "run          0" << endl;
        fout_in_stable << "undump       1" << endl << endl;
        fout_in_stable.close();

    os << "mpirun -np " << opts.omp << " lmp_custom <in.aneal.SiO";
    system(os.str().c_str());
    os.str(empty);

    os << "mdmake --convert mdl --lammps SiO2 -c aneal.final >stable.txt 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);
  }
  else if(opts.atomcomp == "am-at-SiO2")
  {
    os << "mdmake -c at-SiO2 -n \"" << nx << " " << ny << " " << nz << "\" >str.txt 2>/dev/null";
    system(os.str().c_str());
    os.str(empty); 

    os << "mdmake --convert lmp -c str.txt >aneal.lmp 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);

    os << "in.aneal.SiO";
    ofstream fout_in_stable(os.str().c_str());
    os.str(empty);

        fout_in_stable << "package      omp " << opts.omp << endl;
        fout_in_stable << "units        metal" << endl;
        fout_in_stable << "boundary     p p p" << endl;
        fout_in_stable << "atom_style   atomic" << endl;
        fout_in_stable << "atom_modify  sort 10000 1.0" << endl;
        fout_in_stable << "read_data    aneal.lmp" << endl;
        fout_in_stable << "pair_style   sw/wfnho" << endl;
        fout_in_stable << "pair_coeff   * * SiO.sw_wfnho Si O" << endl;
        fout_in_stable << "neighbor     4.0 bin" << endl;
        fout_in_stable << "neigh_modify every 1 delay 0 check yes" << endl;
        fout_in_stable << "velocity     all create 12000 458127641 mom yes rot yes dist gaussian" << endl;
        fout_in_stable << "timestep     0.0001" << endl;//1fs
        fout_in_stable << "fix          1 all nvt temp 12000 1000 0.05" << endl;
        fout_in_stable << "variable     Odensity equal atoms*2.0/3.0/vol" << endl;
        fout_in_stable << "thermo_style custom step temp ke pe etotal press vol density v_Odensity" << endl;
        fout_in_stable << "thermo       1000" << endl;
        fout_in_stable << "run          100000" << endl;
        fout_in_stable << "unfix        1" << endl;
        fout_in_stable << "velocity     all scale 1000" << endl;
        fout_in_stable << "fix          1 all npt temp 1000 300 0.05 x 0.0 0.0 0.5 y 0.0 0.0 0.5 z 0.0 0.0 0.5" << endl;
        fout_in_stable << "run          100000" << endl;
        fout_in_stable << "unfix        1" << endl;
        fout_in_stable << "fix          1 all nve" << endl;
        fout_in_stable << "fix          2 all box/relax aniso 0.0 fixedpoint 0.0 0.0 0.0" << endl;
        fout_in_stable << "min_style    cg" << endl;
        fout_in_stable << "minimize     1e-25 1e-25 50000 100000" << endl;
        fout_in_stable << "dump         1 all custom 1 aneal.final id type xs ys zs" << endl;
        fout_in_stable << "dump_modify  1 sort id" << endl;
        fout_in_stable << "run          0" << endl;
        fout_in_stable << "undump       1" << endl << endl;
        fout_in_stable.close();

    os << "mpirun -np " << opts.omp << " lmp_custom <in.aneal.SiO";
    system(os.str().c_str());
    os.str(empty);

    os << "mdmake --convert mdl --lammps SiO2 -c aneal.final >stable.txt 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);
  }
  else if(opts.atomcomp == "am-bt-SiO2")
  {
    os << "mdmake -c bt-SiO2 -n \"" << nx << " " << ny << " " << nz << "\" >str.txt 2>/dev/null";
    system(os.str().c_str());
    os.str(empty); 

    os << "mdmake --convert lmp -c str.txt >aneal.lmp 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);

    os << "in.aneal.SiO";
    ofstream fout_in_stable(os.str().c_str());
    os.str(empty);

        fout_in_stable << "package      omp " << opts.omp << endl;
        fout_in_stable << "units        metal" << endl;
        fout_in_stable << "boundary     p p p" << endl;
        fout_in_stable << "atom_style   atomic" << endl;
        fout_in_stable << "atom_modify  sort 10000 1.0" << endl;
        fout_in_stable << "read_data    aneal.lmp" << endl;
        fout_in_stable << "pair_style   sw/wfnho" << endl;
        fout_in_stable << "pair_coeff   * * SiO.sw_wfnho Si O" << endl;
        fout_in_stable << "neighbor     4.0 bin" << endl;
        fout_in_stable << "neigh_modify every 1 delay 0 check yes" << endl;
        fout_in_stable << "velocity     all create 12000 458127641 mom yes rot yes dist gaussian" << endl;
        fout_in_stable << "timestep     0.0001" << endl;//1fs
        fout_in_stable << "fix          1 all nvt temp 12000 1000 0.05" << endl;
        fout_in_stable << "variable     Odensity equal atoms*2.0/3.0/vol" << endl;
        fout_in_stable << "thermo_style custom step temp ke pe etotal press vol density v_Odensity" << endl;
        fout_in_stable << "thermo       1000" << endl;
        fout_in_stable << "run          100000" << endl;
        fout_in_stable << "unfix        1" << endl;
        fout_in_stable << "velocity     all scale 1000" << endl;
        fout_in_stable << "fix          1 all npt temp 1000 300 0.05 x 0.0 0.0 0.5 y 0.0 0.0 0.5 z 0.0 0.0 0.5" << endl;
        fout_in_stable << "run          100000" << endl;
        fout_in_stable << "unfix        1" << endl;
        fout_in_stable << "fix          1 all nve" << endl;
        fout_in_stable << "fix          2 all box/relax aniso 0.0 fixedpoint 0.0 0.0 0.0" << endl;
        fout_in_stable << "min_style    cg" << endl;
        fout_in_stable << "minimize     1e-25 1e-25 50000 100000" << endl;
        fout_in_stable << "dump         1 all custom 1 aneal.final id type xs ys zs" << endl;
        fout_in_stable << "dump_modify  1 sort id" << endl;
        fout_in_stable << "run          0" << endl;
        fout_in_stable << "undump       1" << endl << endl;
        fout_in_stable.close();

    os << "mpirun -np " << opts.omp << " lmp_custom <in.aneal.SiO";
    system(os.str().c_str());
    os.str(empty);

    os << "mdmake --convert mdl --lammps SiO2 -c aneal.final >stable.txt 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);
  }
  else if(opts.atomcomp == "am-st-SiO2")
  {
    os << "mdmake -c st-SiO2 -n \"" << nx << " " << ny << " " << nz << "\" >str.txt 2>/dev/null";
    system(os.str().c_str());
    os.str(empty); 

    os << "mdmake --convert lmp -c str.txt >aneal.lmp 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);

    os << "in.aneal.SiO";
    ofstream fout_in_stable(os.str().c_str());
    os.str(empty);

        fout_in_stable << "package      omp " << opts.omp << endl;
        fout_in_stable << "units        metal" << endl;
        fout_in_stable << "boundary     p p p" << endl;
        fout_in_stable << "atom_style   atomic" << endl;
        fout_in_stable << "atom_modify  sort 10000 1.0" << endl;
        fout_in_stable << "read_data    aneal.lmp" << endl;
        fout_in_stable << "pair_style   sw/wfnho" << endl;
        fout_in_stable << "pair_coeff   * * SiO.sw_wfnho Si O" << endl;
        fout_in_stable << "neighbor     4.0 bin" << endl;
        fout_in_stable << "neigh_modify every 1 delay 0 check yes" << endl;
        fout_in_stable << "velocity     all create 12000 458127641 mom yes rot yes dist gaussian" << endl;
        fout_in_stable << "timestep     0.0001" << endl;//1fs
        fout_in_stable << "fix          1 all nvt temp 12000 1000 0.05" << endl;
        fout_in_stable << "variable     Odensity equal atoms*2.0/3.0/vol" << endl;
        fout_in_stable << "thermo_style custom step temp ke pe etotal press vol density v_Odensity" << endl;
        fout_in_stable << "thermo       1000" << endl;
        fout_in_stable << "run          100000" << endl;
        fout_in_stable << "unfix        1" << endl;
        fout_in_stable << "velocity     all scale 1000" << endl;
        fout_in_stable << "fix          1 all npt temp 1000 300 0.05 x 0.0 0.0 0.5 y 0.0 0.0 0.5 z 0.0 0.0 0.5" << endl;
        fout_in_stable << "run          100000" << endl;
        fout_in_stable << "unfix        1" << endl;
        fout_in_stable << "fix          1 all nve" << endl;
        fout_in_stable << "fix          2 all box/relax aniso 0.0 fixedpoint 0.0 0.0 0.0" << endl;
        fout_in_stable << "min_style    cg" << endl;
        fout_in_stable << "minimize     1e-25 1e-25 50000 100000" << endl;
        fout_in_stable << "dump         1 all custom 1 aneal.final id type xs ys zs" << endl;
        fout_in_stable << "dump_modify  1 sort id" << endl;
        fout_in_stable << "run          0" << endl;
        fout_in_stable << "undump       1" << endl << endl;
        fout_in_stable.close();

    os << "mpirun -np " << opts.omp << " lmp_custom <in.aneal.SiO";
    system(os.str().c_str());
    os.str(empty);

    os << "mdmake --convert mdl --lammps SiO2 -c aneal.final >stable.txt 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);
  }
  else if(opts.atomcomp == "am-ac-SiO2")
  {
    os << "mdmake -c ac-SiO2 -n \"" << nx << " " << ny << " " << nz << "\" >str.txt 2>/dev/null";
    system(os.str().c_str());
    os.str(empty); 

    os << "mdmake --convert lmp -c str.txt >aneal.lmp 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);

    os << "in.aneal.SiO";
    ofstream fout_in_stable(os.str().c_str());
    os.str(empty);

        fout_in_stable << "package      omp " << opts.omp << endl;
        fout_in_stable << "units        metal" << endl;
        fout_in_stable << "boundary     p p p" << endl;
        fout_in_stable << "atom_style   atomic" << endl;
        fout_in_stable << "atom_modify  sort 10000 1.0" << endl;
        fout_in_stable << "read_data    aneal.lmp" << endl;
        fout_in_stable << "pair_style   sw/wfnho" << endl;
        fout_in_stable << "pair_coeff   * * SiO.sw_wfnho Si O" << endl;
        fout_in_stable << "neighbor     4.0 bin" << endl;
        fout_in_stable << "neigh_modify every 1 delay 0 check yes" << endl;
        fout_in_stable << "velocity     all create 12000 458127641 mom yes rot yes dist gaussian" << endl;
        fout_in_stable << "timestep     0.0001" << endl;//1fs
        fout_in_stable << "fix          1 all nvt temp 12000 1000 0.05" << endl;
        fout_in_stable << "variable     Odensity equal atoms*2.0/3.0/vol" << endl;
        fout_in_stable << "thermo_style custom step temp ke pe etotal press vol density v_Odensity" << endl;
        fout_in_stable << "thermo       1000" << endl;
        fout_in_stable << "run          100000" << endl;
        fout_in_stable << "unfix        1" << endl;
        fout_in_stable << "velocity     all scale 1000" << endl;
        fout_in_stable << "fix          1 all npt temp 1000 300 0.05 x 0.0 0.0 0.5 y 0.0 0.0 0.5 z 0.0 0.0 0.5" << endl;
        fout_in_stable << "run          100000" << endl;
        fout_in_stable << "unfix        1" << endl;
        fout_in_stable << "fix          1 all nve" << endl;
        fout_in_stable << "fix          2 all box/relax aniso 0.0 fixedpoint 0.0 0.0 0.0" << endl;
        fout_in_stable << "min_style    cg" << endl;
        fout_in_stable << "minimize     1e-25 1e-25 50000 100000" << endl;
        fout_in_stable << "dump         1 all custom 1 aneal.final id type xs ys zs" << endl;
        fout_in_stable << "dump_modify  1 sort id" << endl;
        fout_in_stable << "run          0" << endl;
        fout_in_stable << "undump       1" << endl << endl;
        fout_in_stable.close();

    os << "mpirun -np " << opts.omp << " lmp_custom <in.aneal.SiO";
    system(os.str().c_str());
    os.str(empty);

    os << "mdmake --convert mdl --lammps SiO2 -c aneal.final >stable.txt 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);
  }
  else if(opts.atomcomp == "am-bc-SiO2")
  {
    os << "mdmake -c bc-SiO2 -n \"" << nx << " " << ny << " " << nz << "\" >str.txt 2>/dev/null";
    system(os.str().c_str());
    os.str(empty); 

    os << "mdmake --convert lmp -c str.txt >aneal.lmp 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);

    os << "in.aneal.SiO";
    ofstream fout_in_stable(os.str().c_str());
    os.str(empty);

        fout_in_stable << "package      omp " << opts.omp << endl;
        fout_in_stable << "units        metal" << endl;
        fout_in_stable << "boundary     p p p" << endl;
        fout_in_stable << "atom_style   atomic" << endl;
        fout_in_stable << "atom_modify  sort 10000 1.0" << endl;
        fout_in_stable << "read_data    aneal.lmp" << endl;
        fout_in_stable << "pair_style   sw/wfnho" << endl;
        fout_in_stable << "pair_coeff   * * SiO.sw_wfnho Si O" << endl;
        fout_in_stable << "neighbor     4.0 bin" << endl;
        fout_in_stable << "neigh_modify every 1 delay 0 check yes" << endl;
        fout_in_stable << "velocity     all create 12000 458127641 mom yes rot yes dist gaussian" << endl;
        fout_in_stable << "timestep     0.0001" << endl;//1fs
        fout_in_stable << "fix          1 all nvt temp 12000 1000 0.05" << endl;
        fout_in_stable << "variable     Odensity equal atoms*2.0/3.0/vol" << endl;
        fout_in_stable << "thermo_style custom step temp ke pe etotal press vol density v_Odensity" << endl;
        fout_in_stable << "thermo       1000" << endl;
        fout_in_stable << "run          100000" << endl;
        fout_in_stable << "unfix        1" << endl;
        fout_in_stable << "velocity     all scale 1000" << endl;
        fout_in_stable << "fix          1 all npt temp 1000 300 0.05 x 0.0 0.0 0.5 y 0.0 0.0 0.5 z 0.0 0.0 0.5" << endl;
        fout_in_stable << "run          100000" << endl;
        fout_in_stable << "unfix        1" << endl;
        fout_in_stable << "fix          1 all nve" << endl;
        fout_in_stable << "fix          2 all box/relax aniso 0.0 fixedpoint 0.0 0.0 0.0" << endl;
        fout_in_stable << "min_style    cg" << endl;
        fout_in_stable << "minimize     1e-25 1e-25 50000 100000" << endl;
        fout_in_stable << "dump         1 all custom 1 aneal.final id type xs ys zs" << endl;
        fout_in_stable << "dump_modify  1 sort id" << endl;
        fout_in_stable << "run          0" << endl;
        fout_in_stable << "undump       1" << endl << endl;
        fout_in_stable.close();

    os << "mpirun -np " << opts.omp << " lmp_custom <in.aneal.SiO";
    system(os.str().c_str());
    os.str(empty);

    os << "mdmake --convert mdl --lammps SiO2 -c aneal.final >stable.txt 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);
  }
  else if(opts.atomcomp == "am-bc-GeO2")
  {
    os << "mdmake -c bc-GeO2 -l -n \"" << nx << " " << ny << " " << nz << "\" >str.txt 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);   

    os << "head.aneal1";
    ofstream fout_in_aneal1(os.str().c_str());
    os.str(empty);

    fout_in_aneal1 << "Gear(5) PRECISE(12) RAND=10 ESW T=7500.0 STEP=(0,20000)" << endl;
    fout_in_aneal1 << "BLOCK W=INF FIXANGLE Q=1.0 dt=0.1fs" << endl;
    fout_in_aneal1.close();

    os << "head.aneal2";
    ofstream fout_in_aneal2(os.str().c_str());
    os.str(empty);

    fout_in_aneal2 << "Gear(5) PRECISE(12) RAND=10 ESW T=1000.0 STEP=(0,20000)" << endl;
    fout_in_aneal2 << "BLOCK W=INF FIXANGLE Q=1.0 dt=0.1fs" << endl;
    fout_in_aneal2.close();

    os << "head.stable";
    ofstream fout_in_stable(os.str().c_str());
    os.str(empty);

    fout_in_stable << "OPT(P) RAND=10 ESWGE T=300.0 STEP=(0,20000)" << endl;
    fout_in_stable << "BLOCK W=INF FIXANGLE Q=1.0" << endl;
    fout_in_stable.close();

    os << "tail -n +3 str.txt | cat head.aneal1 - >aneal1.mdl";
    system(os.str().c_str());
    os.str(empty);

    os << "mdlabo < aneal1.mdl > aneal1.final";
    system(os.str().c_str());
    os.str(empty);

    os << "mdmake --convert mdl -c aneal1.final -l >aneal2.txt 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);

    os << "tail -n +3 aneal2.txt | cat head.aneal2 - >aneal2.mdl";
    system(os.str().c_str());
    os.str(empty);

    os << "mdlabo < aneal2.mdl > aneal2.final";
    system(os.str().c_str());
    os.str(empty);

    os << "mdmake --convert mdl -c aneal2.final -l >opt.txt 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);

    os << "tail -n +3 opt.txt | cat head.aneal2 - >opt.mdl";
    system(os.str().c_str());
    os.str(empty);

    os << "mdlabo < opt.mdl > stable.final";
    system(os.str().c_str());
    os.str(empty);

    os << "mdmake --convert mdl -c stable.final -l >stable.txt 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);
  }

  if(0){//test
    os << "mdmake --convert xyz -c stable.txt >test.xyz 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);
  }

  return 0;
}

