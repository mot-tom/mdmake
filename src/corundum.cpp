#include <iostream>
#include <sstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <math.h>
#include "opt.h"
using namespace std;
int UNIT_ATOM_NUM_CORUNDUM_Al = 12;
int UNIT_ATOM_NUM_CORUNDUM_O = 18;

extern int shape(string shape, double x, double y, double z);
 
int corundum(int argc, char **argv, optpara opts){
        string atm, atm2, atm3, atm4;
        double conc2 = 0, conc3 = 0, conc4 = 0;
        double lattice, lattice_a = 4.9133, lattice_c = 5.4053;//lattice_costant[A]
        double c11, c12, c44;
        int nx = 4, ny = 4, nz = 1;

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

        int atmall = 30*nx*ny*nz;
        int atmnumAl = 12*nx*ny*nz;
        int atmnum2 = 0, atmnum3 = 0, atmnum4 = 0;
        double reconc2 = 0, reconc3 = 0, reconc4 = 0;

        double sxx = 0, syy = 0, szz = 0, syz = 0, szx = 0, sxy = 0;//stress[GPa]
        double exx = 0, eyy = 0, ezz = 0, eyz = 0, ezx = 0, exy = 0;//strain

        int i, j, k, l, m, n;
        double lx, ly, lz, txy, txz, tyx, tyz, tzx, tzy;

        string hess;

        cerr << "mdpotential " << opts.potential << endl;
        cerr << "atomcomp " << opts.atomcomp << endl;

        if(opts.atomcomp == "a-Al2O3") {
          atm = "Al";
          lattice_a = 4.759;
          lattice_c = 12.991;
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

        Coordinate atom_al[12] = {
                {0, 0, 0.148}, {0.333333333, 0.666666667, 0.018667}, {0.666666667, 0.333333333, 0.185333},
                {0, 0, 0.352}, {0.333333333, 0.666666667, 0.314667}, {0.666666667, 0.333333333, 0.481333},
                {0, 0, 0.648}, {0.333333333, 0.666666667, 0.518667}, {0.666666667, 0.333333333, 0.685333},
                {0, 0, 0.852}, {0.333333333, 0.666666667, 0.814667}, {0.666666667, 0.333333333, 0.981333}
        };
        Coordinate atom_o[18] = {
                {0.360667, 0.333333, 0.083333}, {0.666667, 0.027333, 0.083333}, {0.972667, 0.639333, 0.083333},
                {0.000000, 0.306000, 0.250000}, {0.306000, 0.000000, 0.250000}, {0.694000, 0.694000, 0.250000},
                {0.027333, 0.666667, 0.416667}, {0.333333, 0.360667, 0.416667}, {0.639333, 0.972667, 0.416667},
                {0.360667, 0.027333, 0.583333}, {0.666667, 0.639333, 0.583333}, {0.972667, 0.333333, 0.583333},
                {0.000000, 0.694000, 0.750000}, {0.306000, 0.306000, 0.750000}, {0.694000, 0.000000, 0.750000},
                {0.333333, 0.972667, 0.916667}, {0.639333, 0.666667, 0.916667}, {0.027333, 0.360667, 0.916667}
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
	 for(m=0; m<atmnumAl; m++){
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
                for(l=0; l<UNIT_ATOM_NUM_CORUNDUM_Al; l++){
                      double ax, ay, az;
                      ax = ((atom_al[l].x+i) + (atom_al[l].y+j)*sin(-30 * M_PI/180))/nx, ay = (atom_al[l].y+j)/ny, az = (atom_al[l].z+k-0.148)/nz;
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
                for(l=0; l<UNIT_ATOM_NUM_CORUNDUM_O; l++){
                      double ax, ay, az;
                      ax = ((atom_o[l].x+i) + (atom_o[l].y+j)*sin(-30 * M_PI/180))/nx, ay = (atom_o[l].y+j)/ny, az = (atom_o[l].z+k-0.148)/nz;
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

