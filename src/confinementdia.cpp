#include <iostream>
#include <sstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <math.h>
#include "opt.h"
using namespace std;
int UNIT_ATOM_NUM_CONFDIA = 8;

extern int shape(string shape, double x, double y, double z);
int conf_confdia(string shape, double x, double y, double z, double conc);
 
int confinementdia(int argc, char **argv, optpara opts){
        string atm, atm2, atm3, atm4;
        double conc2 = 0, conc3 = 0, conc4 = 0;
        double lattice = 5.431;//lattice_costant[A]
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

        int atmall = 8*nx*ny*nz;
        int atmnum2 = 0, atmnum3 = 0, atmnum4 = 0;
        double reconc2 = 0, reconc3 = 0, reconc4 = 0;

        double sxx = 0, syy = 0, szz = 0, syz = 0, szx = 0, sxy = 0;//stress[GPa]
        double exx = 0, eyy = 0, ezz = 0, eyz = 0, ezx = 0, exy = 0;//strain

        int i, j, k, l, m;
        double lx, ly, lz, txy, txz, tyx, tyz, tzx, tzy;

        string hess;

        cerr << "mdpotential " << opts.potential << endl;
        cerr << "atomcomp " << opts.atomcomp << endl;
        cerr << "atomnum " << atmall << " atoms" << endl;

        if(opts.atomcomp == "Si-in-Ge") {
          atm = "Si";
          atm2 = "Ge";
          if(opts.atomratio.length() > 0) {
            istringstream conc_iss(opts.atomratio.c_str());
            conc_iss >> conc2;
           }
          else conc2 = 0.5;
          cerr << "atm2     " << atm2 << endl;
          cerr << "atm2conc " << conc2 << endl;
          cerr << "rand     " << opts.seed << endl;
          lattice = 5.431 + 0.2*conc2 + 0.027*conc2*conc2;
          c11 = 0.00767769330941004 - 0.00211568198181548*conc2;
          c12 = -0.00213584937950066 + 0.000539603259448393*conc2;
          c44 = 0.0125628140703518 - 0.00240724580988776*conc2;
        }
        else if(opts.atomcomp == "Ge-in-Si") {
          atm = "Ge";
          atm2 = "Si";
          if(opts.atomratio.length() > 0) {
            istringstream conc_iss(opts.atomratio.c_str());
            conc_iss >> conc2;
           }
          else conc2 = 0.5;
          cerr << "atm2     " << atm2 << endl;
          cerr << "atm2conc " << conc2 << endl;
          cerr << "rand     " << opts.seed << endl;
          lattice = 5.431 + 0.2*conc2 + 0.027*conc2*conc2;
          c11 = 0.00767769330941004 - 0.00211568198181548*conc2;
          c12 = -0.00213584937950066 + 0.000539603259448393*conc2;
          c44 = 0.0125628140703518 - 0.00240724580988776*conc2;
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

	 cout << "Gear(5) RAND=10 " << opts.potential << " T=300.0 STEP=(0,1000100)" << endl;
        cout << "BLOCK FIXANGLE W=INF Q=INF LOGSTEP=10 PRINTVELOCITY dt=0.1fs scratch(100)" << endl << endl;

        if(opts.stress.length() > 0 || opts.strain.length() > 0) cout << "straind ";
        cout << opts.atomcomp;
        if(atmnum2) cout << " " << atm2 << "=" << reconc2 << " Rand=" << opts.seed;
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

	 vector<string> atmlist1, atmlist2;
	 for(m=0; m<atmall*2; m++){
		if(m<atmnum4) {
                if(m%2) atmlist1.push_back(atm4);
                else atmlist2.push_back(atm4);
                }
		else if(m<atmnum4+atmnum3) {
                if(m%2) atmlist1.push_back(atm3);
                else atmlist2.push_back(atm3);
                }
		else {
                if(m%2) atmlist1.push_back(atm);
                else atmlist2.push_back(atm2);
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

        int delatm = 0;
        atmnum2 = 0;
	 m = 0;
        for(i=0; i<nx; i++){
          for(j=0; j<ny; j++){
            for(k=0; k<nz; k++){
                for(l=0; l<UNIT_ATOM_NUM_CONFDIA; l++){
                      double ax = (atom[l].x+i)/nx, ay = (atom[l].y+j)/ny, az = (atom[l].z+k)/nz;
                      if(conf_confdia(opts.shape, ax, ay, az, 1 - conc2))
                        cout << atmlist1[m] << " " << ax << " " << ay << " " << az << " " << endl;
                      else if(shape(opts.shape, ax, ay, az)){
                        cout << atmlist2[m] << " " << ax << " " << ay << " " << az << " " << endl;
                        atmnum2++;
                         }
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
        cerr << "replaced " << atmnum2 << " atoms " << atm << "->" << atm2 << endl;
        cerr << atm  << " num " << atmall - delatm - atmnum2 << " atoms" << endl;
        cerr << atm2 << " num " << atmnum2 << " atoms " << endl;

        return 0;
}

int conf_confdia(string shape, double x, double y, double z, double conc){
  if(shape == "sphere" && ((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) + (z-0.5)*(z-0.5)) > 0.25*conc) return 0;
  if(shape == "cylinder" && ((y-0.5)*(y-0.5) + (z-0.5)*(z-0.5)) > 0.25*conc) return 0;
  if(shape == "cylinder_y" && ((x-0.5)*(x-0.5) + (z-0.5)*(z-0.5)) > 0.25*conc) return 0;
  if(shape == "cylinder_z" && ((y-0.5)*(y-0.5) + (x-0.5)*(x-0.5)) > 0.25*conc) return 0;
  if(shape == "cuboid" && ((x-0.5)*(x-0.5) > 0.25*sqrt(conc) || (y-0.5)*(y-0.5) > 0.25*sqrt(conc) || (z-0.5)*(z-0.5) > 0.25*sqrt(conc))) return 0;

  //normalized wire
  if(shape == "cuboid_xn" && ((y-0.5)*(y-0.5) > 0.25/4.0*conc || (z-0.5)*(z-0.5) > 0.25/4.0*conc)) return 0;
  if(shape == "cuboid_yn" && ((x-0.5)*(x-0.5) > 0.25/4.0*conc || (z-0.5)*(z-0.5) > 0.25/4.0*conc)) return 0;
  if(shape == "cuboid_zn" && ((x-0.5)*(x-0.5) > 0.25/4.0*conc || (y-0.5)*(y-0.5) > 0.25/4.0*conc)) return 0;
  if(shape == "cylinder_xn" && ((y-0.5)*(y-0.5) + (z-0.5)*(z-0.5)) > 0.25/M_PI*conc) return 0;
  if(shape == "cylinder_yn" && ((x-0.5)*(x-0.5) + (z-0.5)*(z-0.5)) > 0.25/M_PI*conc) return 0;
  if(shape == "cylinder_zn" && ((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5)) > 0.25/M_PI*conc) return 0;
  if(shape == "prism_xn" && (z-0.3+((0.5-0.166667)*(1-sqrt(conc))) > y
  || z-0.3+((0.5-0.166667)*(1-sqrt(conc))) > (1-y) || z < 0.3+(0.166667*(1-sqrt(conc))))) return 0;
  if(shape == "prism_yn" && (z-0.3+((0.5-0.166667)*(1-sqrt(conc))) > x
  || z-0.3+((0.5-0.166667)*(1-sqrt(conc))) > (1-x) || z < 0.3+(0.166667*(1-sqrt(conc))))) return 0;
  if(shape == "prism_zn" && (x-0.3+((0.5-0.166667)*(1-sqrt(conc))) > y
  || x-0.3+((0.5-0.166667)*(1-sqrt(conc))) > (1-y) || x < 0.3+(0.166667*(1-sqrt(conc))))) return 0;
  return 1;
}
