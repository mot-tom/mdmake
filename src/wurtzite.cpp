#include <iostream>
#include <sstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <math.h>
#include "opt.h"
using namespace std;

extern int shape(string shape, double x, double y, double z);
 
int wurtzite(int argc, char **argv, optpara opts){
        string atm, atm2, atm3, atm4;
        string stkptn;
        int stkfreq;
        double conc2 = 0, conc3 = 0, conc4 = 0;
        double lattice, lattice_a = 3.0730, lattice_c = 2.5133;//lattice_costant[A]
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

        int atmall = 2*nx*ny*nz;
        int atmnum2 = 0, atmnum3 = 0, atmnum4 = 0;
        double reconc2 = 0, reconc3 = 0, reconc4 = 0;

        double sxx = 0, syy = 0, szz = 0, syz = 0, szx = 0, sxy = 0;//stress[GPa]
        double exx = 0, eyy = 0, ezz = 0, eyz = 0, ezx = 0, exy = 0;//strain

        int i, j, k, l, m, n;
        double lx, ly, lz, txy, txz, tyx, tyz, tzx, tzy;

        string hess;

        cerr << "mdpotential " << opts.potential << endl;
        cerr << "atomcomp " << opts.atomcomp << endl;

        if(opts.atomcomp == "2H-SiC") {
          stkfreq = 2;
          nz *= stkfreq; atmall *= stkfreq;
          stkptn = "AB";
          atm = "Si";
          atm2 = "C ";
          lattice_a = 3.0730;
          lattice_c = 5.04 / 2.0;
          c11 = 0.01;
          c12 = 0;
          c44 = 0.02;
        }
        else if(opts.atomcomp == "2H-Si") {
          stkfreq = 2;
          nz *= stkfreq; atmall *= stkfreq;
          stkptn = "AB";
          atm = "Si";
          atm2 = "Si";
          lattice_a = 3.84030;
          lattice_c = 9.40676 / 3.0;
          c11 = 0.01;
          c12 = 0;
          c44 = 0.02;
        }
        else if(opts.atomcomp == "2H-SiCX") {
          stkfreq = 2;
          nz *= stkfreq; atmall *= stkfreq;
          stkptn = "AB";
          atm = "Si";
          atm2 = "C ";
          atm3 = "X ";
          if(opts.atomratio.length() > 0) {
            istringstream conc_iss(opts.atomratio.c_str());
            conc_iss >> conc3;
           }
          else conc3 = 0.01;
          atmnum3 = ((int)(atmall*conc3));
	   reconc3 = (double)atmnum3/atmall;
          cerr << "atm3     " << atm3 << endl;
          cerr << "atm3num  " << atmnum3 << endl;
          cerr << "atm3conc " << reconc3 << endl;
          cerr << "rand     " << opts.seed << endl;
          lattice_a = 3.044;
          lattice_c = 4.995 / 2.0;
          c11 = 0.01;
          c12 = 0;
          c44 = 0.02;
        }
        else if(opts.atomcomp == "3H-SiGe-alloy") {
          stkfreq = 3;
          nz *= stkfreq; atmall *= stkfreq;
          stkptn = "ABC";
          atm = "Si";
          atm2 = "Si";
          atm3 = "Ge";
          if(opts.atomratio.length() > 0) {
            istringstream conc_iss(opts.atomratio.c_str());
            conc_iss >> conc3;
           }
          else conc3 = 0.5;
          atmnum3 = (int)(atmall*conc3);
	   reconc3 = (double)atmnum3/atmall;
          cerr << "atm2     " << atm3 << endl;
          cerr << "atm2num  " << atmnum3 << endl;
          cerr << "atm2conc " << reconc3 << endl;
          cerr << "rand     " << opts.seed << endl;
          lattice = 5.431 + 0.2*reconc3 + 0.027*reconc3*reconc3;
          lattice_a = lattice * sqrt(2) / 2.0;
          lattice_c = lattice * sqrt(3) / 3.0;
          c11 = 0.00767769330941004 - 0.00211568198181548*reconc3;
          c12 = -0.00213584937950066 + 0.000539603259448393*reconc3;
          c44 = 0.0125628140703518 - 0.00240724580988776*reconc3;
        }
        else if(opts.atomcomp == "3H-GeSiSn-alloy") {
          stkfreq = 3;
          nz *= stkfreq; atmall *= stkfreq;
          stkptn = "ABC";
          atm = "Ge";
          atm2 = "Ge";
          atm3 = "Si";
          atm4 = "Sn";
          if(opts.atomratio.length() > 0) {
            istringstream conc_iss(opts.atomratio.c_str());
            conc_iss >> conc3 >> conc4;
           }
          else {conc3 = 0.333; conc4 = 0.333;}
          atmnum3 = (int)(atmall*conc3);
	   reconc3 = (double)atmnum3/atmall;
          atmnum4 = (int)(atmall*conc4);
	   reconc4 = (double)atmnum4/atmall;
          cerr << "atm2     " << atm3 << endl;
          cerr << "atm2num  " << atmnum3 << endl;
          cerr << "atm2conc " << reconc3 << endl;
          cerr << "atm3     " << atm4 << endl;
          cerr << "atm3num  " << atmnum4 << endl;
          cerr << "atm3conc " << reconc4 << endl;
          cerr << "rand     " << opts.seed << endl;
          lattice = 5.658 - 0.227*reconc3 + 0.831*reconc4;
          lattice_a = lattice * sqrt(2) / 2.0;
          lattice_c = lattice * sqrt(3) / 3.0;
          c11 = 0.01;
          c12 = 0;
          c44 = 0.02;
        }
        else if(opts.atomcomp == "3H-SiC") {
          stkfreq = 3;
          nz *= stkfreq; atmall *= stkfreq;
          stkptn = "ABC";
          atm = "Si";
          atm2 = "C ";
          lattice_a = 3.04834;
          lattice_c = 7.46687 / 3.0;
          c11 = 0.01;
          c12 = 0;
          c44 = 0.02;
        }
        else if(opts.atomcomp == "3H-Si") {
          stkfreq = 3;
          nz *= stkfreq; atmall *= stkfreq;
          stkptn = "ABC";
          atm = "Si";
          atm2 = "Si";
          lattice_a = 3.84030;
          lattice_c = 9.40676 / 3.0;
          c11 = 0.01;
          c12 = 0;
          c44 = 0.02;
        }
        else if(opts.atomcomp == "4H-SiC") {
          stkfreq = 4;
          nz *= stkfreq; atmall *= stkfreq;
          stkptn = "ABCB";
          atm = "Si";
          atm2 = "C ";
          lattice_a = 3.045;
          lattice_c = 9.971 / 4.0;
          c11 = 0.01;
          c12 = 0;
          c44 = 0.02;
        }
        else if(opts.atomcomp == "4H-Si") {
          stkfreq = 4;
          nz *= stkfreq; atmall *= stkfreq;
          stkptn = "ABCB";
          atm = "Si";
          atm2 = "Si";
          lattice_a = 3.84030;
          lattice_c = 9.40676 / 3.0;
          c11 = 0.01;
          c12 = 0;
          c44 = 0.02;
        }
        else if(opts.atomcomp == "6H-SiC") {
          stkfreq = 6;
          nz *= stkfreq; atmall *= stkfreq;
          stkptn = "ABCACB";
          atm = "Si";
          atm2 = "C ";
          lattice_a = 3.045;
          lattice_c = 14.966 / 6.0;
          c11 = 0.01;
          c12 = 0;
          c44 = 0.02;
        }
        else if(opts.atomcomp == "8H-SiC") {
          stkfreq = 8;
          nz *= stkfreq; atmall *= stkfreq;
          stkptn = "ABCABACB";
          atm = "Si";
          atm2 = "C ";
          lattice_a = 3.047;
          lattice_c = 19.926 / 8.0;
          c11 = 0.01;
          c12 = 0;
          c44 = 0.02;
        }
        else if(opts.atomcomp == "10H-SiC") {
          stkfreq = 10;
          nz *= stkfreq; atmall *= stkfreq;
          stkptn = "ABCABCBACB";
          atm = "Si";
          atm2 = "C ";
          lattice_a = 3.047;
          lattice_c = (24.913 + 0.002) / 10.0;
          c11 = 0.01;
          c12 = 0;
          c44 = 0.02;
        }
        else if(opts.atomcomp == "12H-SiC") {
          stkfreq = 12;
          nz *= stkfreq; atmall *= stkfreq;
          stkptn = "ABCABCACBACB";
          atm = "Si";
          atm2 = "C ";
          lattice_a = 3.047;
          lattice_c = 29.883 / 12.0;
          c11 = 0.01;
          c12 = 0;
          c44 = 0.02;
        }
        else if(opts.atomcomp == "15R-SiC") {
          stkfreq = 15;
          nz *= stkfreq; atmall *= stkfreq;
          stkptn = "ABCACBCABACABCB";
          atm = "Si";
          atm2 = "C ";
          lattice_a = 3.047;
          lattice_c = 37.7 / 15.0;
          c11 = 0.01;
          c12 = 0;
          c44 = 0.02;
        }
        else if(opts.atomcomp == "18H-SiC") {
          stkfreq = 18;
          nz *= stkfreq; atmall *= stkfreq;
          stkptn = "ABCACBABCBABCBABCB";
          atm = "Si";
          atm2 = "C ";
          lattice_a = 3.046;
          lattice_c = 44.860 / 18.0;
          c11 = 0.01;
          c12 = 0;
          c44 = 0.02;
        }
        else if(opts.atomcomp == "2H-AlN") {
          stkfreq = 2;
          nz *= stkfreq; atmall *= stkfreq;
          stkptn = "AB";
          atm = "Al";
          atm2 = "N ";
          lattice_a = 3.111;
          lattice_c = 4.978 / 2.0;
          c11 = 0.01;
          c12 = 0;
          c44 = 0.02;
        }
        else if(opts.atomcomp == "2H-BN") {
          stkfreq = 2;
          nz *= stkfreq; atmall *= stkfreq;
          stkptn = "AB";
          atm = "B ";
          atm2 = "N ";
          lattice_a = 2.55;
          lattice_c = 4.17 / 2.0;
          c11 = 0.01;
          c12 = 0;
          c44 = 0.02;
        }
        else if(opts.atomcomp == "2H-BeO") {
          stkfreq = 2;
          nz *= stkfreq; atmall *= stkfreq;
          stkptn = "AB";
          atm = "Be";
          atm2 = "O ";
          lattice_a = 2.66;
          lattice_c = 4.37 / 2.0;
          c11 = 0.01;
          c12 = 0;
          c44 = 0.02;
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

        Coordinate atom1_a = {0, 0, 0};
        Coordinate atom2_a = {0, 0, 0.75};
        Coordinate atom1_b = {0.666667, 0.333333, 0};
        Coordinate atom2_b = {0.666667, 0.333333, 0.75};
        Coordinate atom1_c = {0.333333, 0.666667, 0};
        Coordinate atom2_c = {0.333333, 0.666667, 0.75};

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

	 vector<string> atmlist1, atmlist2;
	 for(m=0; m<atmall; m++){
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
	 m = 0; n = 0;
        for(i=0; i<nx; i++){
          for(j=0; j<ny; j++){
            for(k=0; k<nz; k++){
              if(stkptn[k%stkfreq] == 'C'){
                      double ax = ((atom1_c.x+i) + (atom1_c.y+j)*sin(-30 * M_PI/180))/nx, ay = (atom1_c.y+j)/ny, az = (atom1_c.z+k)/nz;
                      if(ax < -1.0e-07) ax += 1; if(ax > -1.0e-07 && ax < 1.0e-07) ax = 0;
                      if(shape(opts.shape, ax, ay, az))
                        cout << atmlist1[m] << " " << ax << " " << ay << " " << az << " " << endl;
                      else delatm++;
                      m++;
                      ax = ((atom2_c.x+i) + (atom2_c.y+j)*sin(-30 * M_PI/180))/nx; ay = (atom2_c.y+j)/ny; az = (atom2_c.z+k)/nz;
                      if(ax < -1.0e-07) ax += 1; if(ax > -1.0e-07 && ax < 1.0e-07) ax = 0;
                      if(shape(opts.shape, ax, ay, az))
                        cout << atmlist2[n] << " " << ax << " " << ay << " " << az << " " << endl;
                      else delatm++;
                      n++;
              }
              else if(stkptn[k%stkfreq] == 'B'){
                      double ax = ((atom1_b.x+i) + (atom1_b.y+j)*sin(-30 * M_PI/180))/nx, ay = (atom1_b.y+j)/ny, az = (atom1_b.z+k)/nz;
                      if(ax < -1.0e-07) ax += 1; if(ax > -1.0e-07 && ax < 1.0e-07) ax = 0;
                      if(shape(opts.shape, ax, ay, az))
                        cout << atmlist1[m] << " " << ax << " " << ay << " " << az << " " << endl;
                      else delatm++;
                      m++;
                      ax = ((atom2_b.x+i) + (atom2_b.y+j)*sin(-30 * M_PI/180))/nx; ay = (atom2_b.y+j)/ny; az = (atom2_b.z+k)/nz;
                      if(ax < -1.0e-07) ax += 1; if(ax > -1.0e-07 && ax < 1.0e-07) ax = 0;
                      if(shape(opts.shape, ax, ay, az))
                        cout << atmlist2[n] << " " << ax << " " << ay << " " << az << " " << endl;
                      else delatm++;
                      n++;
              }
              else{
                      double ax = ((atom1_a.x+i) + (atom1_a.y+j)*sin(-30 * M_PI/180))/nx, ay = (atom1_a.y+j)/ny, az = (atom1_a.z+k)/nz;
                      if(ax < -1.0e-07) ax += 1; if(ax > -1.0e-07 && ax < 1.0e-07) ax = 0;
                      if(shape(opts.shape, ax, ay, az))
                        cout << atmlist1[m] << " " << ax << " " << ay << " " << az << " " << endl;
                      else delatm++;
                      m++;
                      ax = ((atom2_a.x+i) + (atom2_a.y+j)*sin(-30 * M_PI/180))/nx; ay = (atom2_a.y+j)/ny; az = (atom2_a.z+k)/nz;
                      if(ax < -1.0e-07) ax += 1; if(ax > -1.0e-07 && ax < 1.0e-07) ax = 0;
                      if(shape(opts.shape, ax, ay, az))
                        cout << atmlist2[n] << " " << ax << " " << ay << " " << az << " " << endl;
                      else delatm++;
                      n++;
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

