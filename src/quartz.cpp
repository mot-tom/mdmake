#include <iostream>
#include <sstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <math.h>
#include "opt.h"
using namespace std;
int UNIT_ATOM_NUM_QUARTZ_Si = 3;
int UNIT_ATOM_NUM_QUARTZ_O = 6;

extern int shape(string shape, double x, double y, double z);
 
int quartz(int argc, char **argv, optpara opts){
        string atm, atm2, atm3, atm4;
        string crystal = "bq";
        double conc2 = 0, conc3 = 0, conc4 = 0;
        double lattice, lattice_a = 4.9133, lattice_c = 5.4053;//lattice_costant[A]
        double c11, c12, c44;
        int nx = 4, ny = 4, nz = 4;

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

        int atmall = 9*nx*ny*nz;
        int atmnumSi = 3*nx*ny*nz;
        int atmnum2 = 0, atmnum3 = 0, atmnum4 = 0;
        double reconc2 = 0, reconc3 = 0, reconc4 = 0;

        double sxx = 0, syy = 0, szz = 0, syz = 0, szx = 0, sxy = 0;//stress[GPa]
        double exx = 0, eyy = 0, ezz = 0, eyz = 0, ezx = 0, exy = 0;//strain

        int i, j, k, l, m, n;
        double lx, ly, lz, txy, txz, tyx, tyz, tzx, tzy;

        string hess;

        cerr << "mdpotential " << opts.potential << endl;
        cerr << "atomcomp " << opts.atomcomp << endl;

        if(opts.atomcomp == "bq-SiO2") {
          crystal = "bq";
          atm = "Si";
          lattice_a = 4.9133;
          lattice_c = 5.4053;
        }
        else if(opts.atomcomp == "bq-GeO2") {
          crystal = "bq";
          atm = "Ge";
          lattice_a = 4.9858;
          lattice_c = 5.6473;
        }
        else if(opts.atomcomp == "aq-SiO2") {
          crystal = "aq";
          atm = "Si";
          lattice_a = 4.9133;
          lattice_c = 5.4053;
        }
        else if(opts.atomcomp == "aq-GeO2") {
          crystal = "aq";
          atm = "Ge";
          lattice_a = 4.9858;
          lattice_c = 5.6473;
        }
        cerr << "atomnum " << atmall << " atoms" << endl;

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
                cerr << "% Eyz=" << eyz << "% Ezx=" << ezx << "% Exy=" << exy << endl;
                exx = 0.01*exx; eyy = 0.01*eyy; ezz = 0.01*ezz;
                eyz = 0.01*eyz; ezx = 0.01*ezx; exy = 0.01*exy;
	}

        lx = nx*lattice_a*(1+exx);
        ly = ny*lattice_a*cos(-30 * M_PI/180)*(1+eyy);
        lz = nz*lattice_c*(1+ezz);
        txy = nx*lattice_a*exy/2;
        txz = nx*lattice_a*ezx/2;
        tyx = ny*lattice_a*cos(-30 * M_PI/180)*exy/2;
        tyz = ny*lattice_a*cos(-30 * M_PI/180)*eyz/2;
        tzx = nz*lattice_c*ezx/2;
        tzy = nz*lattice_c*eyz/2;

        Coordinate atom_si_bq[3] = {
                {0, 0, 0}, {0.5, 0, 0.666667}, {0, 0.5, 0.333333}
        };
        Coordinate atom_o_bq[6] = {
                {0.296667, 0.086667, 0.833333}, {0.086667, 0.296667, 0.166667}, {0.716667, 0.3, 0.5}, {0.3, 0.716667, 0.5},
                {0.703333, 0.913333, 0.833333}, {0.913333, 0.703333, 0.166667}
        };

        Coordinate atom_si_aq[3] = {
                {0, 0, 0}, {0.466667, 0.933333, 0.666667}, {0.933333, 0.466667, 0.333333}
        };
        Coordinate atom_o_aq[6] = {
                {0.336667, 0.066667, 0.876667}, {0.066667, 0.336667, 0.123333}, {0.623333, 0.206667, 0.456667}, {0.206667, 0.623333, 0.543333},
                {0.733333, 0.873333, 0.79}, {0.873333, 0.733333, 0.21}
        };

	 cout << "Gear(5) RAND=10 " << opts.potential << " T=300.0 STEP=(0,1000100)" << endl;
        cout << "BLOCK FIXANGLE W=INF Q=INF LOGSTEP=10 PRINTVELOCITY dt=0.1fs scratch(100)" << endl << endl;

        if(opts.stress.length() > 0 || opts.strain.length() > 0) cout << "straind ";
        cout << opts.atomcomp;
        if(atmnum4) cout << " " << atm4 << "=" << reconc4 << " Rand=" << opts.seed;
        else if(atmnum3) cout << " " << atm3 << "=" << reconc3 << " Rand=" << opts.seed;
        if(opts.stress.length() > 0){
                cout << " Sxx=" << sxx;
                cout << "GPa Syy=" << syy << "GPa Szz=" << szz;
                cout << "GPa Syz=" << syz << "GPa Szx=" << szx << "GPa Sxy=" << sxy << "GPa";
	}
        if(opts.strain.length() > 0){
                cout << " Exx=" << exx;
                cout << " Eyy=" << eyy << " Ezz=" << ezz;
                cout << " Eyz=" << eyz << " Ezx=" << ezx << " Exy=" << exy << "%";
	}
        cout << " structure" << endl << endl;

        cout << lx << " " << txy << " " << txz << opts.pbc_x << endl;
        cout << tyx << " " << ly << " " << tyz << opts.pbc_y << endl;
        cout << tzx << " " << tzy << " " << lz << opts.pbc_z << endl << endl;

	 vector<string> atmlist;
	 for(m=0; m<atmnumSi; m++){
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
                for(l=0; l<UNIT_ATOM_NUM_QUARTZ_Si; l++){
                      double ax, ay, az;
                      if(crystal == "bq")
                        ax = ((atom_si_bq[l].x+i) + (atom_si_bq[l].y+j)*sin(-30 * M_PI/180))/nx, ay = (atom_si_bq[l].y+j)/ny, az = (atom_si_bq[l].z+k)/nz;
                      else if(crystal == "aq")
                        ax = ((atom_si_aq[l].x+i) + (atom_si_aq[l].y+j)*sin(-30 * M_PI/180))/nx, ay = (atom_si_aq[l].y+j)/ny, az = (atom_si_aq[l].z+k)/nz;
                      if(ax < -1.0e-07) ax += 1; if(ax > -1.0e-07 && ax < 1.0e-07) ax = 0;
                      if(shape(opts.shape, ax, ay, az))
                        cout << atmlist[m] << " " << ax << " " << ay << " " << az << " " << endl;
                      else delatm++;
                      m++;
                }
            }
          }
        }
        for(i=0; i<nx; i++){
          for(j=0; j<ny; j++){
            for(k=0; k<nz; k++){
                for(l=0; l<UNIT_ATOM_NUM_QUARTZ_O; l++){
                      double ax, ay, az;
                      if(crystal == "bq")
                        ax = ((atom_o_bq[l].x+i) + (atom_o_bq[l].y+j)*sin(-30 * M_PI/180))/nx, ay = (atom_o_bq[l].y+j)/ny, az = (atom_o_bq[l].z+k)/nz;
                      else if(crystal == "aq")
                        ax = ((atom_o_aq[l].x+i) + (atom_o_aq[l].y+j)*sin(-30 * M_PI/180))/nx, ay = (atom_o_aq[l].y+j)/ny, az = (atom_o_aq[l].z+k)/nz;
                      if(ax < -1.0e-07) ax += 1; if(ax > -1.0e-07 && ax < 1.0e-07) ax = 0;
                      if(shape(opts.shape, ax, ay, az))
                        cout << "O" << " " << ax << " " << ay << " " << az << " " << endl;
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

        return 0;
}

