#include <iostream>
#include <sstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <math.h>
#include "opt.h"
using namespace std;
int UNIT_ATOM1_NUM_DC = 2;
int UNIT_ATOM2_NUM_DC = 1;

extern int shape(string shape, double x, double y, double z);
 
int dichalcogenide(int argc, char **argv, optpara opts){
        string atmf, atms, atm2, atm3, atm4;
        double conc2 = 0, conc3 = 0, conc4 = 0;
        double lattice_a = 3.127, lattice_c = 5.8;//lattice_costant[A]
        double c11, c12, c44;
        int nx = 8, ny = 8, nz = 4;

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

        int atmall = 3*nx*ny*nz;
        int atmnum2 = 0, atmnum3 = 0, atmnum4 = 0;
        double reconc2 = 0, reconc3 = 0, reconc4 = 0;

        double sxx = 0, syy = 0, szz = 0, syz = 0, szx = 0, sxy = 0;//stress[GPa]
        double exx = 0, eyy = 0, ezz = 0, eyz = 0, ezx = 0, exy = 0;//strain

        int i, j, k, l, m, n, o;
        double lx, ly, lz, txy, txz, tyx, tyz, tzx, tzy;

        string hess;

        cerr << "mdpotential " << opts.potential << endl;
        cerr << "atomcomp " << opts.atomcomp << endl;
        cerr << "atomnum " << atmall << " atoms" << endl;

        if(opts.atomcomp == "MoS2") {
          atmf = "S1";
          atms = "S2";
          atm2 = "Mo";
          lattice_a = 3.127;
          lattice_c = 12.066 / 2.0;
          c11 = 0.01;
          c12 = 0;
          c44 = 0.02;
        }
        else if(opts.atomcomp == "MoSX") {
          atmf = "S1";
          atms = "S2";
          atm2 = "Mo";
          atm3 = "X ";
          if(opts.atomratio.length() > 0) {
            istringstream conc_iss(opts.atomratio.c_str());
            conc_iss >> conc3;
           }
          else conc3 = 0.01;
          if(conc3 > 0.666666) {
            cerr << "Wrong atom composition" << endl;
            return 1;
           }
          atmnum3 = ((int)(atmall*conc3));
	   reconc3 = (double)atmnum3/atmall;
          cerr << "atm3     " << atm3 << endl;
          cerr << "atm3num  " << atmnum3 << endl;
          cerr << "atm3conc " << reconc3 << endl;
          cerr << "rand     " << opts.seed << endl;
          lattice_a = 3.127;
          lattice_c = 12.066 / 2.0;
          c11 = 0.01;
          c12 = 0;
          c44 = 0.02;
        }
        else if(opts.atomcomp == "MoXS2") {
          atmf = "S1";
          atms = "S2";
          atm2 = "Mo";
          atm4 = "X ";
          if(opts.atomratio.length() > 0) {
            istringstream conc_iss(opts.atomratio.c_str());
            conc_iss >> conc4;
           }
          else conc4 = 0.01;
          if(conc4 > 0.333333) {
            cerr << "Wrong atom composition" << endl;
            return 1;
           }
          atmnum4 = ((int)(atmall*conc4));
	   reconc4 = (double)atmnum4/atmall;
          cerr << "atm3     " << atm4 << endl;
          cerr << "atm3num  " << atmnum4 << endl;
          cerr << "atm3conc " << reconc4 << endl;
          cerr << "rand     " << opts.seed << endl;
          lattice_a = 3.127;
          lattice_c = 12.066 / 2.0;
          c11 = 0.01;
          c12 = 0;
          c44 = 0.02;
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

        lx = nx*lattice_a*(1+exx);
        ly = ny*lattice_a*cos(-30 * M_PI/180)*(1+eyy);
        lz = nz*lattice_c*(1+ezz);
        txy = nx*lattice_a*exy/2;
        txz = nx*lattice_a*ezx/2;
        tyx = ny*lattice_a*cos(-30 * M_PI/180)*exy/2;
        tyz = ny*lattice_a*cos(-30 * M_PI/180)*eyz/2;
        tzx = nz*lattice_c*ezx/2;
        tzy = nz*lattice_c*eyz/2;

        Coordinate atom1_a[2] = {
                {0, 0, 0.25}, {0, 0, 0.75}
         };
        Coordinate atom1_b[2] = {
                {0.333333, 0.666667, 0.25}, {0.333333, 0.666667, 0.75}
         };
        Coordinate atom2_a    = {0.333333, 0.666667, 0.5};
        Coordinate atom2_b    = {0, 0, 0.5};

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
                cout << " Eyz=" << eyz << " Ezx=" << ezx << " Exy=" << exy;
	}
        cout << " structure" << endl << endl;

        cout << lx << " " << txy << " " << txz << opts.pbc_x << endl;
        cout << tyx << " " << ly << " " << tyz << opts.pbc_y << endl;
        cout << tzx << " " << tzy << " " << lz << opts.pbc_z << endl << endl;

	 vector<string> atmlist1, atmlist2, atmlist3;
	 for(n=0; n<nx*ny; n++){
	   for(m=0; m<3*nz; m++){
              if(m%6 == 3 || m%6 == 5){
                  if(m+n*(3*nz)<atmnum3) atmlist3.push_back(atm3);
                  else atmlist3.push_back(atms);
                }
              else if(m%6 == 0 || m%6 == 2){
                if(m+n*(3*nz)<atmnum3) atmlist1.push_back(atm3);
                else atmlist1.push_back(atmf);
                }
		else {
                if(m+n*(3*nz)<atmnum4) atmlist2.push_back(atm4);
                else atmlist2.push_back(atm2);
                }
	   }
	 }

	 srand(opts.seed);
	 for(m=0; m<atmlist1.size(); m++){
		int n = rand()%(m+1);
		string t = atmlist1[m];
		atmlist1[m] = atmlist1[n];
		atmlist1[n] = t;
	 }
	 for(m=0; m<atmlist2.size(); m++){
		int n = rand()%(m+1);
		string t = atmlist2[m];
		atmlist2[m] = atmlist2[n];
		atmlist2[n] = t;
	 }
	 for(m=0; m<atmlist3.size(); m++){
		int n = rand()%(m+1);
		string t = atmlist3[m];
		atmlist3[m] = atmlist3[n];
		atmlist3[n] = t;
	 }

	 m = 0; n = 0; o = 0;
        for(i=0; i<nx; i++){
          for(j=0; j<ny; j++){
            for(k=0; k<nz; k++){
              if(k%2){
                      double ax = ((atom1_b[0].x+i) + (atom1_b[0].y+j)*sin(-30 * M_PI/180))/nx, ay = (atom1_b[0].y+j)/ny, az = (atom1_b[0].z+k)/nz;
                      if(ax < -1.0e-07) ax += 1; if(ax > -1.0e-07 && ax < 1.0e-07) ax = 0;
                      if(shape(opts.shape, ax, ay, az))
                        cout << atmlist3[o] << " " << ax << " " << ay << " " << az << " " << endl;
                      o++;
                      ax = ((atom2_b.x+i) + (atom2_b.y+j)*sin(-30 * M_PI/180))/nx, ay = (atom2_b.y+j)/ny, az = (atom2_b.z+k)/nz;
                      if(ax < -1.0e-07) ax += 1; if(ax > -1.0e-07 && ax < 1.0e-07) ax = 0;
                      if(shape(opts.shape, ax, ay, az))
                        cout << atmlist2[n] << " " << ax << " " << ay << " " << az << " " << endl;
                      n++;
                      ax = ((atom1_b[1].x+i) + (atom1_b[1].y+j)*sin(-30 * M_PI/180))/nx; ay = (atom1_b[1].y+j)/ny; az = (atom1_b[1].z+k)/nz;
                      if(ax < -1.0e-07) ax += 1; if(ax > -1.0e-07 && ax < 1.0e-07) ax = 0;
                      if(shape(opts.shape, ax, ay, az))
                        cout << atmlist3[o] << " " << ax << " " << ay << " " << az << " " << endl;
                      o++;
              }
             else{
                      double ax = ((atom1_a[0].x+i) + (atom1_a[0].y+j)*sin(-30 * M_PI/180))/nx, ay = (atom1_a[0].y+j)/ny, az = (atom1_a[0].z+k)/nz;
                      if(ax < -1.0e-07) ax += 1; if(ax > -1.0e-07 && ax < 1.0e-07) ax = 0;
                      if(shape(opts.shape, ax, ay, az))
                        cout << atmlist1[m] << " " << ax << " " << ay << " " << az << " " << endl;
                      m++;
                      ax = ((atom2_a.x+i) + (atom2_a.y+j)*sin(-30 * M_PI/180))/nx; ay = (atom2_a.y+j)/ny; az = (atom2_a.z+k)/nz;
                      if(ax < -1.0e-07) ax += 1; if(ax > -1.0e-07 && ax < 1.0e-07) ax = 0;
                      if(shape(opts.shape, ax, ay, az))
                        cout << atmlist2[n] << " " << ax << " " << ay << " " << az << " " << endl;
                      n++;
                      ax = ((atom1_a[1].x+i) + (atom1_a[1].y+j)*sin(-30 * M_PI/180))/nx; ay = (atom1_a[1].y+j)/ny; az = (atom1_a[1].z+k)/nz;
                      if(ax < -1.0e-07) ax += 1; if(ax > -1.0e-07 && ax < 1.0e-07) ax = 0;
                      if(shape(opts.shape, ax, ay, az))
                        cout << atmlist1[m] << " " << ax << " " << ay << " " << az << " " << endl;
                      m++;
              }
            }
          }
        }

        return 0;
}

