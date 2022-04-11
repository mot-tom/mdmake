#include <iostream>
#include <sstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <math.h>
#include "opt.h"
using namespace std;
int UNIT_ATOM_NUM_DIA_OXIDE_Si = 8;
int UNIT_ATOM_NUM_DIA_OXIDE_O = 16;

extern int shape(string shape, double x, double y, double z);

int diamond_oxide(int argc, char **argv, optpara opts){
        string atm, atm2, atm3, atm4;
        string crystal = "bc";
        double conc2 = 0, conc3 = 0, conc4 = 0;
        double lattice = 7.166;//lattice_costant[A]
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

        int atmall = 24*nx*ny*nz;
        int atmnumSi = 8*nx*ny*nz;
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

        if(opts.atomcomp == "bc-SiO2") {
          atm = "Si";
          crystal = "bc";
          lattice = 7.166;
        }
        else if(opts.atomcomp == "bc-GeO2") {
          atm = "Ge";
          crystal = "bc";
          lattice = 7.859;//7.166*1.70/1.55
        }
        else if(opts.atomcomp == "bc2-SiO2") {
          atm = "Si";
          crystal = "bc2";
          lattice = 7.166;
        }
        else if(opts.atomcomp == "bc2-GeO2") {
          atm = "Ge";
          crystal = "bc2";
          lattice = 7.859;//7.166*1.70/1.55
        }
        else if(opts.atomcomp == "bc3-SiO2") {
          atm = "Si";
          crystal = "bc3";
          lattice = 7.166;
        }
        else if(opts.atomcomp == "bc3-GeO2") {
          atm = "Ge";
          crystal = "bc3";
          lattice = 7.859;//7.166*1.70/1.55
        }
        else if(opts.atomcomp == "ac-SiO2") {
          atm = "Si";
          crystal = "ac";
          lattice = 6.92563;
        }
        else if(opts.atomcomp == "ac-GeO2") {
          atm = "Ge";
          crystal = "ac";
          lattice = 6.92563*1.70/1.55;
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

        Coordinate atom_si_bc[8] = {
                {0, 0, 0}, {0.25, 0.25, 0.25}, {0, 0.5, 0.5}, {0.25, 0.75, 0.75},
                {0.5, 0, 0.5}, {0.75, 0.25, 0.75}, {0.5, 0.5, 0}, {0.75, 0.75, 0.25}
        };
        Coordinate atom_o_bc[16] = {
                {0.125, 0.125, 0.125}, {0.125, 0.375, 0.375}, {0.125, 0.625, 0.625}, {0.125, 0.875, 0.875},
                {0.375, 0.625, 0.875}, {0.375, 0.875, 0.625}, {0.375, 0.125, 0.375}, {0.375, 0.375, 0.125},
                {0.625, 0.625, 0.125}, {0.625, 0.875, 0.375}, {0.625, 0.125, 0.625}, {0.625, 0.375, 0.875},
                {0.875, 0.125, 0.875}, {0.875, 0.375, 0.625}, {0.875, 0.625, 0.375}, {0.875, 0.875, 0.125}
        };

        Coordinate atom_si_bc2[8] = {
                {0, 0, 0}, {0.25, 0.25, 0.25}, {0, 0.5, 0.5}, {0.25, 0.75, 0.75},
                {0.5, 0, 0.5}, {0.75, 0.25, 0.75}, {0.5, 0.5, 0}, {0.75, 0.75, 0.25}
        };
        Coordinate atom_o_bc2[16] = {
                {0.125, 0.125, 0.125}, {0.125, 0.423333, 0.326667}, {0.076667, 0.673333, 0.625}, {0.173333, 0.875, 0.923333},
                {0.375, 0.576667, 0.826667}, {0.375, 0.875, 0.625}, {0.326667, 0.125, 0.423333}, {0.423333, 0.326667, 0.125},
                {0.673333, 0.625, 0.076667}, {0.576667, 0.826667, 0.375}, {0.625, 0.076667, 0.673333}, {0.625, 0.375, 0.875},
                {0.923333, 0.173333, 0.875}, {0.826667, 0.375, 0.576667}, {0.875, 0.625, 0.375}, {0.875, 0.923333, 0.173333}
        };

        Coordinate atom_si_bc3[8] = {
                {0, 0, 0}, {0.25, 0.25, 0.25}, {0, 0.5, 0.5}, {0.25, 0.75, 0.75},
                {0.5, 0, 0.5}, {0.75, 0.25, 0.75}, {0.5, 0.5, 0}, {0.75, 0.75, 0.25}
        };
        Coordinate atom_o_bc3[16] = {
                {0.076667, 0.173333, 0.125}, {0.173333, 0.423333, 0.375}, {0.076667, 0.673333, 0.625}, {0.173333, 0.923333, 0.875},
                {0.326667, 0.576667, 0.875}, {0.423333, 0.826667, 0.625}, {0.326667, 0.076667, 0.375}, {0.423333, 0.326667, 0.125},
                {0.576667, 0.673333, 0.125}, {0.673333, 0.923333, 0.375}, {0.576667, 0.173333, 0.625}, {0.673333, 0.423333, 0.875},
                {0.826667, 0.076667, 0.875}, {0.923333, 0.326667, 0.625}, {0.826667, 0.576667, 0.375}, {0.923333, 0.826667, 0.125}
        };

        Coordinate atom_si_ac[8] = {
                {0, 0, 0}, {0.2, 0.3, 0.25}, {0.9, 0.5, 0.5}, {0.2, 0.7, 0.75},
                {0.4, 0, 0.5}, {0.7, 0.2, 0.75}, {0.5, 0.5, 0}, {0.7, 0.8, 0.25}
        };
        Coordinate atom_o_ac[16] = {
                {0.133333, 0.175, 0.075}, {0.025, 0.433333, 0.325}, {0.025, 0.566666, 0.675}, {0.133333, 0.825, 0.925},
                {0.375, 0.566666, 0.825}, {0.266666, 0.825, 0.566666}, {0.266666, 0.175, 0.433333}, {0.375, 0.433333, 0.175},
                {0.633333, 0.666666, 0.075}, {0.525, 0.925, 0.325}, {0.525, 0.075, 0.675}, {0.633333, 0.333333, 0.925},
                {0.875, 0.075, 0.825}, {0.766666, 0.333333, 0.566666}, {0.766666, 0.666666, 0.433333}, {0.875, 0.925, 0.175}
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
        for(i=0; i<nx; i++){
          for(j=0; j<ny; j++){
            for(k=0; k<nz; k++){
                for(l=0; l<UNIT_ATOM_NUM_DIA_OXIDE_Si; l++){
                      double ax, ay, az;
                      if(crystal == "bc")
                        ax = (atom_si_bc[l].x+i)/nx, ay = (atom_si_bc[l].y+j)/ny, az = (atom_si_bc[l].z+k)/nz;
                      else if(crystal == "bc2")
                        ax = (atom_si_bc2[l].x+i)/nx, ay = (atom_si_bc2[l].y+j)/ny, az = (atom_si_bc2[l].z+k)/nz;
                      else if(crystal == "bc3")
                        ax = (atom_si_bc3[l].x+i)/nx, ay = (atom_si_bc3[l].y+j)/ny, az = (atom_si_bc3[l].z+k)/nz;
                      else if(crystal == "ac")
                        ax = (atom_si_ac[l].x+i)/nx, ay = (atom_si_ac[l].y+j)/ny, az = (atom_si_ac[l].z+k)/nz;
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
                for(l=0; l<UNIT_ATOM_NUM_DIA_OXIDE_O; l++){
                      double ax, ay, az;
                      if(crystal == "bc")
                        ax = (atom_o_bc[l].x+i)/nx, ay = (atom_o_bc[l].y+j)/ny, az = (atom_o_bc[l].z+k)/nz;
                      else if(crystal == "bc2")
                        ax = (atom_o_bc2[l].x+i)/nx, ay = (atom_o_bc2[l].y+j)/ny, az = (atom_o_bc2[l].z+k)/nz;
                      else if(crystal == "bc3")
                        ax = (atom_o_bc3[l].x+i)/nx, ay = (atom_o_bc3[l].y+j)/ny, az = (atom_o_bc3[l].z+k)/nz;
                      else if(crystal == "ac")
                        ax = (atom_o_ac[l].x+i)/nx, ay = (atom_o_ac[l].y+j)/ny, az = (atom_o_ac[l].z+k)/nz;
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

