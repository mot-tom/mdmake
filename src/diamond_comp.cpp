#include <iostream>
#include <sstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <math.h>
#include "opt.h"
using namespace std;
int UNIT_ATOM_NUM_DIA_COMP_Si = 4;

extern int shape(string shape, double x, double y, double z);

int diamond_comp(int argc, char **argv, optpara opts){
        string atm, atm2, atm3, atm4;
        string stkptn;
        int stknum, stknumSi;
        Coordinate stkpos[30];
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

        int atmall = 4*nx*ny*nz;
        int atmnumSi = 4*nx*ny*nz;
        int atmnum2 = 0, atmnum3 = 0, atmnum4 = 0;
        double reconc2 = 0, reconc3 = 0, reconc4 = 0;
        double sxx = 0, syy = 0, szz = 0, syz = 0, szx = 0, sxy = 0;//stress[GPa]
        double exx = 0, eyy = 0, ezz = 0, eyz = 0, ezx = 0, exy = 0;//strain

        int i, j, k, l, m;
        double lx, ly, lz, txy, txz, tyx, tyz, tzx, tzy;

        string hess;

        if(opts.atomcomp == "g-Si3N4") {
          stknum = 14;
          stknumSi = 6;
          atmall *= stknum;
          atmnumSi *= stknumSi;
          stkptn = "FFFFFFNNNNNNNN";
          stkpos[0].x = 0.000, stkpos[0].y = 0.000, stkpos[0].z =0.000;
          stkpos[1].x = 0.250, stkpos[1].y = 0.250, stkpos[1].z =0.250;
          stkpos[2].x = 0.625, stkpos[2].y = 0.625, stkpos[2].z =0.625;
          stkpos[3].x = 0.875, stkpos[3].y = 0.875, stkpos[3].z =0.625;
          stkpos[4].x = 0.875, stkpos[4].y = 0.625, stkpos[4].z =0.875;
          stkpos[5].x = 0.625, stkpos[5].y = 0.875, stkpos[5].z =0.875;

          stkpos[6].x = 0.379, stkpos[6].y = 0.379, stkpos[6].z =0.379;
          stkpos[7].x = 0.129, stkpos[7].y = 0.129, stkpos[7].z =0.371;
          stkpos[8].x = 0.129, stkpos[8].y = 0.371, stkpos[8].z =0.129;
          stkpos[9].x = 0.371, stkpos[9].y = 0.129, stkpos[9].z =0.129;
          stkpos[10].x = 0.629, stkpos[10].y = 0.629, stkpos[10].z =0.871;
          stkpos[11].x = 0.629, stkpos[11].y = 0.871, stkpos[11].z =0.629;
          stkpos[12].x = 0.871, stkpos[12].y = 0.629, stkpos[12].z =0.629;
          stkpos[13].x = 0.871, stkpos[13].y = 0.871, stkpos[13].z =0.871;
          atm = "Si";
          lattice = 7.166;
        }
        else if(opts.atomcomp == "Y2O3") {
          stknum = 24;
          stknumSi = 8;
          atmall *= stknum;
          atmnumSi *= stknumSi;
          stkptn = "FFFFFFFFoooooooooooooooo";
          stkpos[0].x = 0.250, stkpos[0].y = 0.000, stkpos[0].z =0.000;
          stkpos[1].x = 0.000, stkpos[1].y = 0.250, stkpos[1].z =0.000;
          stkpos[2].x = 0.000, stkpos[2].y = 0.000, stkpos[2].z =0.250;
          stkpos[3].x = 0.250, stkpos[3].y = 0.250, stkpos[3].z =0.250;
          stkpos[4].x = 0.500, stkpos[4].y = 0.500, stkpos[4].z =0.750;
          stkpos[5].x = 0.500, stkpos[5].y = 0.750, stkpos[5].z =0.500;
          stkpos[6].x = 0.750, stkpos[6].y = 0.500, stkpos[6].z =0.500;
          stkpos[7].x = 0.750, stkpos[7].y = 0.750, stkpos[7].z =0.750;

          stkpos[8].x = 0.125, stkpos[8].y = 0.125, stkpos[8].z =0.125;
          stkpos[9].x = 0.375, stkpos[9].y = 0.375, stkpos[9].z =0.375;
          stkpos[10].x = 0.625, stkpos[10].y = 0.625, stkpos[10].z =0.625;
          stkpos[11].x = 0.875, stkpos[11].y = 0.875, stkpos[11].z =0.875;
          stkpos[12].x = 0.375, stkpos[12].y = 0.625, stkpos[12].z =0.125;
          stkpos[13].x = 0.625, stkpos[13].y = 0.875, stkpos[13].z =0.125;
          stkpos[14].x = 0.875, stkpos[14].y = 0.375, stkpos[14].z =0.125;
          stkpos[15].x = 0.125, stkpos[15].y = 0.875, stkpos[15].z =0.375;
          stkpos[16].x = 0.625, stkpos[16].y = 0.125, stkpos[16].z =0.375;
          stkpos[17].x = 0.875, stkpos[17].y = 0.625, stkpos[17].z =0.375;
          stkpos[18].x = 0.125, stkpos[18].y = 0.375, stkpos[18].z =0.625;
          stkpos[19].x = 0.375, stkpos[19].y = 0.875, stkpos[19].z =0.625;
          stkpos[20].x = 0.875, stkpos[20].y = 0.125, stkpos[20].z =0.625;
          stkpos[21].x = 0.125, stkpos[21].y = 0.625, stkpos[21].z =0.875;
          stkpos[22].x = 0.375, stkpos[22].y = 0.125, stkpos[22].z =0.875;
          stkpos[23].x = 0.625, stkpos[23].y = 0.375, stkpos[23].z =0.875;
          atm = "Y ";
          lattice = 10.5981;
        }
        else if(opts.atomcomp == "C-ZrO2") {
          stknum = 3;
          stknumSi = 1;
          atmall *= stknum;
          atmnumSi *= stknumSi;
          stkptn = "FOO";
          stkpos[0].x = 0.000, stkpos[0].y = 0.000, stkpos[0].z =0.000;
          stkpos[1].x = 0.250, stkpos[1].y = 0.250, stkpos[1].z =0.250;
          stkpos[2].x = 0.750, stkpos[2].y = 0.750, stkpos[2].z =0.750;
          atm = "Zr";
          lattice = 5.1289;
        }
        else if(opts.atomcomp == "C-SnO2") {
          stknum = 3;
          stknumSi = 1;
          atmall *= stknum;
          atmnumSi *= stknumSi;
          stkptn = "FOO";
          stkpos[0].x = 0.000, stkpos[0].y = 0.000, stkpos[0].z =0.000;
          stkpos[1].x = 0.250, stkpos[1].y = 0.250, stkpos[1].z =0.250;
          stkpos[2].x = 0.750, stkpos[2].y = 0.750, stkpos[2].z =0.750;
          atm = "Sn";
          lattice = 4.925;
        }
        else if(opts.atomcomp == "C-SiO2") {
          stknum = 3;
          stknumSi = 1;
          atmall *= stknum;
          atmnumSi *= stknumSi;
          stkptn = "FOO";
          stkpos[0].x = 0.000, stkpos[0].y = 0.000, stkpos[0].z =0.000;
          stkpos[1].x = 0.250, stkpos[1].y = 0.250, stkpos[1].z =0.250;
          stkpos[2].x = 0.750, stkpos[2].y = 0.750, stkpos[2].z =0.750;
          atm = "Si";
          lattice = 4.925;
        }
        else if(opts.atomcomp == "bc-SiO2-test") {
          stknum = 6;
          stknumSi = 2;
          atmall *= stknum;
          atmnumSi *= stknumSi;
          stkptn = "FFOOOO";
          stkpos[0].x = 0.000, stkpos[0].y = 0.000, stkpos[0].z =0.000;
          stkpos[1].x = 0.250, stkpos[1].y = 0.250, stkpos[1].z =0.250;
          stkpos[2].x = 0.125, stkpos[2].y = 0.125, stkpos[2].z =0.125;
          stkpos[3].x = 0.125, stkpos[3].y = 0.375, stkpos[3].z =0.375;
          stkpos[4].x = 0.375, stkpos[4].y = 0.125, stkpos[4].z =0.375;
          stkpos[5].x = 0.375, stkpos[5].y = 0.375, stkpos[5].z =0.125;
          atm = "Si";
          lattice = 7.166;
        }
        else if(opts.atomcomp == "Si-test") {
          stknum = 2;
          stknumSi = 2;
          atmall *= stknum;
          atmnumSi *= stknumSi;
          stkptn = "FF";
          stkpos[0].x = 0.000, stkpos[0].y = 0.000, stkpos[0].z =0.000;
          stkpos[1].x = 0.250, stkpos[1].y = 0.250, stkpos[1].z =0.250;
          atm = "Si";
          lattice = 5.431;
        }

        cerr << "mdpotential " << opts.potential << endl;
        cerr << "atomcomp " << opts.atomcomp << endl;
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

        Coordinate atom[4] = {
                {0, 0, 0}, {0, 0.5, 0.5}, 
                {0.5, 0, 0.5}, {0.5, 0.5, 0}
        };

	 cout << "Gear(5) RAND=10 " << opts.potential << " T=300.0 STEP=(0,1000100)" << endl;
        cout << "BLOCK FIXANGLE W=INF Q=INF LOGSTEP=10 PRINTVELOCITY dt=0.1fs scratch(100)" << endl << endl;

        if(opts.stress.length() > 0 || opts.strain.length() > 0) cout << "straind ";
        cout << opts.atomcomp;
        if(atmnum2) cout << " " << atm2 << "=" << reconc2;
        if(atmnum3) cout << " " << atm3 << "=" << reconc3;
        if(atmnum2+atmnum3) cout << " Rand=" << opts.seed;
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
       for(int n=0; n<stknum; n++){
        for(i=0; i<nx; i++){
          for(j=0; j<ny; j++){
            for(k=0; k<nz; k++){
                for(l=0; l<UNIT_ATOM_NUM_DIA_COMP_Si; l++){
                    double tx, ty, tz, ax, ay, az;
                      tx = atom[l].x+stkpos[n].x, ty = atom[l].y+stkpos[n].y, tz = atom[l].z+stkpos[n].z;
                      if(tx>=1.0) tx -= 1.0;
                      if(ty>=1.0) ty -= 1.0;
                      if(tz>=1.0) tz -= 1.0;
                      ax = (tx+i)/nx, ay = (ty+j)/ny, az = (tz+k)/nz;
                      if(shape(opts.shape, ax, ay, az)){
                        if(stkptn[n] == 'O') cout << "O " << " " << ax << " " << ay << " " << az << " " << endl;
                        else if(stkptn[n] == 'o' && l) cout << "O " << " " << ax << " " << ay << " " << az << " " << endl;
                        else if(stkptn[n] == 'N') cout << "N " << " " << ax << " " << ay << " " << az << " " << endl;
                        else if(stkptn[n] == 'F'){cout << atmlist[m] << " " << ax << " " << ay << " " << az << " " << endl; m++;}
                        else delatm++;
                            }
                      else delatm++;
                    }
                }
            }
          }
        }
        if(delatm){
          cerr << "deleted " << delatm << " atoms" << endl;
          cerr << "remains " << atmall - delatm << " atoms" << endl;
         }

        return 0;
}

