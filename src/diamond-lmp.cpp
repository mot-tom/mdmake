#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <math.h>
#include "opt.h"
using namespace std;
int UNIT_ATOM_NUM_DIA_LMP = 8;

/* from Mdlabo */
typedef enum {Si, Ge, C, Sn, Dummy} AtomIonType;

 typedef enum {ESWGE_OPTION_133DEG, ESWGE_OPTION_144DEG,
               ESWGE_OPTION_B3LYP, ESWGE_OPTION_MP2,
               ESWGE_OPTION_HF, ESWGE_OPTION_B3LYP_HF, ESWGE_OPTION_MP2_HF,
               ESWGE_OPTION_TEST} ESWGE_option;
/* from Mdlabo */

extern int shape(string shape, double x, double y, double z);
 
int diamond_lmp(int argc, char **argv, optpara opts, string lmpname){
        string pbc_x, pbc_y, pbc_z;
        string atm, atm2, atm3, atm4;
        string atmname, atmname2, atmname3, atmname4;
        double atmtype, neigh;
        double mass = 0, mass2 = 0, mass3 = 0, mass4 = 0;
        double conc2 = 0, conc3 = 0, conc4 = 0;
        double lattice = 5.431;//lattice_costant[A]
        double c11, c12, c44;
        int nx = 4, ny = 4, nz = 4;

        ostringstream osfile_input;
        ostringstream osfile_struct;
        osfile_input << "input." << lmpname;
        osfile_struct << "struct." << lmpname;
        ofstream fout_in(osfile_input.str().c_str());
        ofstream fout_st(osfile_struct.str().c_str());
        if(!fout_in){
          cerr << "Can't open " << osfile_input.str() << "!" << endl;
          return 1;
         }
        else cout << "\"" << osfile_input.str() << "\" is made." << endl;
        if(!fout_st){
          cerr << "Can't open " << osfile_struct.str() << "!" << endl;
          return 1;
         }
        else cout << "\"" << osfile_struct.str() << "\" is made." << endl;

        ESWGE_option eswge_option;
        if(opts.potential == "ESWGE(hf)")
          eswge_option = ESWGE_OPTION_HF;
        if(opts.potential == "ESWGE(b3lyp)")
          eswge_option = ESWGE_OPTION_B3LYP;
        if(opts.potential == "ESWGE(b3lyp+hf)")
          eswge_option = ESWGE_OPTION_B3LYP_HF;
        if(opts.potential == "ESWGE(mp2)")
          eswge_option = ESWGE_OPTION_MP2;
        if(opts.potential == "ESWGE(mp2+hf)")
          eswge_option = ESWGE_OPTION_MP2_HF;

        if(opts.pbc_x==" N") pbc_x = "f";
        else pbc_x = "p";
        if(opts.pbc_y==" N") pbc_y = "f";
        else pbc_y = "p";
        if(opts.pbc_z==" N") pbc_z = "f";
        else pbc_z = "p";

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

        int atmall = 8*nx*ny*nz;
        int atmnum2 = 0, atmnum3 = 0, atmnum4 = 0;
        double reconc2 = 0, reconc3 = 0, reconc4 = 0;

        double sxx = 0, syy = 0, szz = 0, syz = 0, szx = 0, sxy = 0;//stress[GPa]
        double exx = 0, eyy = 0, ezz = 0, eyz = 0, ezx = 0, exy = 0;//strain

        int i, j, k, l, m;
        double lx, ly, lz, txy, txz, tyx, tyz, tzx, tzy;

        string hess;

        cerr << "atomcomp " << opts.atomcomp << endl;
        cerr << "atomnum " << atmall << " atoms" << endl;

        if(opts.atomcomp == "SiGe") {
          atm = "1";//Si
          atm2 = "2";//Ge
          atmtype = 2;
          if(opts.atomratio.length() > 0) {
            istringstream conc_iss(opts.atomratio.c_str());
            conc_iss >> conc2;
           }
          else conc2 = 0.5;
          atmnum2 = (int)(atmall*conc2);
	   reconc2 = (double)atmnum2/atmall;
          cerr << "atm2     " << atm2 << endl;
          cerr << "atm2num  " << atmnum2 << endl;
          cerr << "atm2conc " << reconc2 << endl;
          cerr << "rand     " << opts.seed << endl;
          mass = 28.0855;
          mass2 = 72.63;
          lattice = 5.431 + 0.2*reconc2 + 0.027*reconc2*reconc2;
          c11 = 0.00767769330941004 - 0.00211568198181548*reconc2;
          c12 = -0.00213584937950066 + 0.000539603259448393*reconc2;
          c44 = 0.0125628140703518 - 0.00240724580988776*reconc2;
        }
        else if(opts.atomcomp == "SiC") {
          atm = "1";//Si
          atm2 = "2";//C
          atmtype = 2;
          if(opts.atomratio.length() > 0) {
            istringstream conc_iss(opts.atomratio.c_str());
            conc_iss >> conc2;
           }
          else conc2 = 0.5;
          atmnum2 = (int)(atmall*conc2);
	   reconc2 = (double)atmnum2/atmall;
          cerr << "atm2     " << atm2 << endl;
          cerr << "atm2num  " << atmnum2 << endl;
          cerr << "atm2conc " << reconc2 << endl;
          cerr << "rand     " << opts.seed << endl;
          mass = 28.0855;
          mass2 = 12.0107;
          lattice = 5.431 - 2.4239*reconc2 + 0.5705*reconc2*reconc2;
        }
        else if(opts.atomcomp == "GeSiSn") {
          atm = "1";
          atm2 = "2";
          atm3 = "3";
          atmtype = 3;
          if(opts.atomratio.length() > 0) {
            istringstream conc_iss(opts.atomratio.c_str());
            conc_iss >> conc2 >> conc3;
           }
          else {conc2 = 0.333; conc3 = 0.333;}
          atmnum2 = (int)(atmall*conc2);
	   reconc2 = (double)atmnum2/atmall;
          atmnum3 = (int)(atmall*conc3);
	   reconc3 = (double)atmnum3/atmall;
          cerr << "atm2     " << atm2 << endl;
          cerr << "atm2num  " << atmnum2 << endl;
          cerr << "atm2conc " << reconc2 << endl;
          cerr << "atm3     " << atm3 << endl;
          cerr << "atm3num  " << atmnum3 << endl;
          cerr << "atm3conc " << reconc3 << endl;
          cerr << "rand     " << opts.seed << endl;
          mass = 72.63;
          mass2 = 28.0855;
          mass3 = 118.710;
          lattice = 5.658 - 0.227*reconc2 + 0.831*reconc3;
        }
        else if(opts.atomcomp == "Si") {
          atm = "1";//Si
          atmname = "Si";
          atmtype = 1;
          neigh = 4.0;
          mass = 28.0855;
          lattice = 5.431;
          c11 = 0.00767769330941004;
          c12 = -0.00213584937950066;
          c44 = 0.0125628140703518;
        }
        else if(opts.atomcomp == "Ge") {
          atm = "1";//Ge
          atmname = "Ge";
          atmtype = 1;
          neigh = 4.0;
          mass = 72.63;
          lattice = 5.658;
          c11 = 0.00979337529122551;
          c12 = -0.00267545263894905;
          c44 = 0.0149700598802395;
        }
        else if(opts.atomcomp == "Sn") {
          atm = "1";//Sn
          atmname = "Sn";
          atmtype = 1;
          neigh = 4.0;
          mass = 118.710;
          lattice = 6.489;
          c11 = 0.0194049414495866;
          c12 = -0.00578397542698767;
          c44 = 0.0276243093922652;
        }
        else if(opts.atomcomp == "C") {
          atm = "1";//C
          atmname = "C";
          atmtype = 1;
          neigh = 4.0;
          mass = 12.0107;
          lattice = 3.567;
          c11 = 0.000949273446777954;
          c12 = -0.0000978469720702131;
          c44 = 0.00173010380622837;
        }

        if(opts.stress.length() > 0){
                istringstream s_iss(opts.stress.c_str());
                s_iss >> sxx >> syy >> szz >> syz >> szx >> sxy;//[GPa]
                cerr << "stress Sxx=" << sxx;
                cerr << "GPa Syy=" << syy << "GPa Szz=" << szz;
                cerr << "GPa Syz=" << syz << "GPa Szx=" << szx << "GPa Sxy=" << sxy << "GPa" << endl;
                exx = c11*sxx + c12*syy + c12*szz;
                eyy = c12*sxx + c11*syy + c12*szz;
                ezz = c12*sxx + c12*syy + c11*szz;
                eyz = c44*syz/2; ezx = c44*szx/2; exy = c44*sxy/2;
	}
        if(opts.strain.length() > 0){
                istringstream e_iss(opts.strain.c_str());
                e_iss >> exx >> eyy >> ezz >> eyz >> ezx >> exy;//[%]
                cerr << "strain Exx=" << exx;
                cerr << "% Eyy=" << eyy << "% Ezz=" << ezz;
                cerr << "% Eyz=" << eyz << "% Ezx=" << ezx << "% Exy=" << exy << "%" << endl;
                exx = 0.01*exx; eyy = 0.01*eyy; ezz = 0.01*ezz;
                eyz = 0.01*eyz; ezx = 0.01*ezx; exy = 0.01*exy;
	}

        lx = nx*lattice*(1+exx);
        ly = ny*lattice*(1+eyy);
        lz = nz*lattice*(1+ezz);
        txy = nx*lattice*exy/2;
        txz = nx*lattice*ezx/2;
        tyx = ny*lattice*exy/2;
        tyz = ny*lattice*eyz/2;
        tzx = nz*lattice*ezx/2;
        tzy = nz*lattice*eyz/2;

        Coordinate atom[8] = {
                {0, 0, 0}, {0.25, 0.25, 0.25}, {0, 0.5, 0.5}, {0.25, 0.75, 0.75},
                {0.5, 0, 0.5}, {0.75, 0.25, 0.75}, {0.5, 0.5, 0}, {0.75, 0.75, 0.25}
        };

        fout_in << "package      omp 1" << endl;
        fout_in << "units        metal" << endl;
        fout_in << "boundary     " << pbc_x << " " << pbc_y << " " << pbc_z << endl;
        fout_in << "atom_style   atomic" << endl;
        fout_in << "atom_modify  sort 10000 1.0" << endl;
        fout_in << "read_data    " << osfile_struct.str() << endl;

        fout_in << "pair_style   sw" << endl;

        fout_in << "neighbor     " << neigh << " bin" << endl;
        fout_in << "neigh_modify every 1 delay 0 check yes" << endl;
        fout_in << "fix          1 all nve" << endl;
        fout_in << "fix          2 all box/relax aniso 0.0 fixedpoint 0.0 0.0 0.0" << endl;
        fout_in << "thermo       10" << endl;
        fout_in << "min_style    cg" << endl;
        fout_in << "minimize     1e-25 1e-25 50000 100000" << endl;
        fout_in << "dump         1 all custom 1 final." << lmpname << " id type xu yu zu" << endl;
        fout_in << "run          0" << endl;
        fout_in << "undump       1" << endl << endl;

        if(opts.stress.length() > 0 || opts.strain.length() > 0) fout_st << "straind ";
        fout_st << opts.atomcomp;
        if(atmnum2) fout_st << " " << atm2 << "=" << reconc2;
        if(atmnum3) fout_st << " " << atm3 << "=" << reconc3;
        if(atmnum2+atmnum3) fout_st << " Rand=" << opts.seed;
        if(opts.stress.length() > 0){
                fout_st << " Sxx=" << sxx;
                fout_st << "GPa Syy=" << syy << "GPa Szz=" << szz;
                fout_st << "GPa Syz=" << syz << "GPa Szx=" << szx << "GPa Sxy=" << sxy << "GPa";
	}
        if(opts.strain.length() > 0){
                fout_st << " Exx=" << exx;
                fout_st << " Eyy=" << eyy << " Ezz=" << ezz;
                fout_st << " Eyz=" << eyz << " Ezx=" << ezx << " Exy=" << exy;
	}
        fout_st << " structure" << endl << endl;

        fout_st << " " << atmall << " atoms" << endl;
        fout_st << " " << atmtype << " atom types" << endl;
        fout_st << " 0.00 " << lx << " xlo xhi" << endl;
        fout_st << " 0.00 " << ly << " ylo yhi" << endl;
        fout_st << " 0.00 " << lz << " zlo zhi" << endl;
        if(txy || txz || tyx || tyz || tzx || tzy)
          fout_st << " " << txy + tyx << " " << txz + tzx << " " << tyz + tzy << " xy xz yz" << endl;
        fout_st << endl;

        fout_st << " Masses" << endl << endl;
        fout_st << "1 " << mass << endl;
        if(mass2) fout_st << "2 " << mass2 << endl;
        if(mass3) fout_st << "3 " << mass3 << endl;
        if(mass4) fout_st << "4 " << mass4 << endl;
        fout_st << endl;

        fout_st << " Atoms" << endl << endl;

	 vector<string> atmlist;
	 for(m=0; m<atmall; m++){
		if(m<atmnum4) atmlist.push_back(atm4);
		else if(m<atmnum4+atmnum3) atmlist.push_back(atm3);
		else if(m<atmnum4+atmnum3+atmnum2) atmlist.push_back(atm2);
		else atmlist.push_back(atm);
	 }
	 srand(opts.seed);
	 for(m=0; m<atmlist.size(); m++){
		int n = rand()%(m+1);
		string t = atmlist[m];
		atmlist[m] = atmlist[n];
		atmlist[n] = t;
	 }

        int delatm = 0;
	 m = 0;
        for(i=0; i<nx; i++){
          for(j=0; j<ny; j++){
            for(k=0; k<nz; k++){
                for(l=0; l<UNIT_ATOM_NUM_DIA_LMP; l++){
                      double ax = (atom[l].x+i)/nx, ay = (atom[l].y+j)/ny, az = (atom[l].z+k)/nz;
                      if(shape(opts.shape, ax, ay, az))
                        fout_st << m+1 << " " << atmlist[m] << " " << ax*lx << " " << ay*ly << " " << az*lz << " " << endl;
                      else delatm++;
                      m++;
                }
            }
          }
        }
        if(opts.shape != "cuboid"){
          cerr << "deleted " << delatm << " atoms" << endl;
          cerr << "remains " << atmall - delatm << " atoms" << endl;
         }

        fout_in.close();
        fout_st.close();

        return 0;
}

