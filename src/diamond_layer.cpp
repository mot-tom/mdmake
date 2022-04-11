#include <iostream>
#include <sstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <math.h>
#include "opt.h"
using namespace std;
int UNIT_ATOM_NUM_DIA_LAYER = 2;

extern int shape(string shape, double x, double y, double z);

int diamond_layer(int argc, char **argv, optpara opts){
        string atm, atm2, atm3, atm4;
        double conc2 = 0, conc3 = 0, conc4 = 0;
        double lattice = 5.431;//lattice_costant[A]
        double c11, c12, c44;
        int nx = 4, ny = 4, nz = 2;
        int one_layer = 0;

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

        int atmall = 2*nx*ny*nz*2;
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

        if(opts.atomcomp == "zl-Si") {
          atm = "Si";
          lattice = 5.431;
          c11 = 0.00767769330941004;
          c12 = -0.00213584937950066;
          c44 = 0.0125628140703518;
        }
        else if(opts.atomcomp == "zl1-Si") {
          atm = "Si";
          one_layer = 1;
          lattice = 5.431;
          c11 = 0.00767769330941004;
          c12 = -0.00213584937950066;
          c44 = 0.0125628140703518;
        }
        else if(opts.atomcomp == "zl-Ge") {
          atm = "Ge";
          lattice = 5.658;
          c11 = 0.00979337529122551;
          c12 = -0.00267545263894905;
          c44 = 0.0149700598802395;
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
        lz = nz*2*lattice/4*(1+ezz)*1;
        txy = nx*lattice*exy/2;
        txz = nx*lattice*ezx/2;
        tyx = ny*lattice*exy/2;
        tyz = ny*lattice*eyz/2;
        tzx = nz*2*lattice/4*ezx/2;
        tzy = nz*2*lattice/4*eyz/2;

        Coordinate atom_zl1[2] = {
                {0, 0, 0}, {0.5, 0.5, 0}
        };
        Coordinate atom_zl2[2] = {
                {0.25, 0.25, 0}, {0.75, 0.75, 0}
        };
        Coordinate atom_zl3[2] = {
                {0, 0.5, 0}, {0.5, 0, 0}
        };
        Coordinate atom_zl4[2] = {
                {0.25, 0.75, 0}, {0.75, 0.25, 0}
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

        int count=0;
        m = 0;
        for(k=0; k<nz; k++){
          for(j=0; j<ny; j++){
            for(i=0; i<nx; i++){
                for(l=0; l<UNIT_ATOM_NUM_DIA_LAYER; l++){
                      double ax, ay, az;
                      if(count%4==0)      ax = (atom_zl4[l].x+i)/nx, ay = (atom_zl4[l].y+j)/ny, az = 0.5-(atom_zl4[l].z-k)/(nz*2);
                      else if(count%4==1) ax = (atom_zl3[l].x+i)/nx, ay = (atom_zl3[l].y+j)/ny, az = 0.5-(atom_zl3[l].z-k)/(nz*2);
                      else if(count%4==2) ax = (atom_zl2[l].x+i)/nx, ay = (atom_zl2[l].y+j)/ny, az = 0.5-(atom_zl2[l].z-k)/(nz*2);
                      else if(count%4==3) ax = (atom_zl1[l].x+i)/nx, ay = (atom_zl1[l].y+j)/ny, az = 0.5-(atom_zl1[l].z-k)/(nz*2);
                      if(k==0) 
                        cout << atmlist[m] << " " << ax << " " << ay << " " << az << " F" << endl;
                      else
                        cout << atmlist[m] << " " << ax << " " << ay << " " << az << " G " << nz - k << endl;
                      m++;
                }
              if(!one_layer){
                for(l=0; l<UNIT_ATOM_NUM_DIA_LAYER; l++){
                      double ax, ay, az;
                      if(count%4==0)      ax = (atom_zl1[l].x+i)/nx, ay = (atom_zl1[l].y+j)/ny, az = 0.5-(atom_zl1[l].z+k+1)/(nz*2);
                      else if(count%4==1) ax = (atom_zl2[l].x+i)/nx, ay = (atom_zl2[l].y+j)/ny, az = 0.5-(atom_zl2[l].z+k+1)/(nz*2);
                      else if(count%4==2) ax = (atom_zl3[l].x+i)/nx, ay = (atom_zl3[l].y+j)/ny, az = 0.5-(atom_zl3[l].z+k+1)/(nz*2);
                      else if(count%4==3) ax = (atom_zl4[l].x+i)/nx, ay = (atom_zl4[l].y+j)/ny, az = 0.5-(atom_zl4[l].z+k+1)/(nz*2);
                      if(k==0) 
                        cout << atmlist[m] << " " << ax << " " << ay << " " << az << " F" << endl;
                      else
                        cout << atmlist[m] << " " << ax << " " << ay << " " << az << " G " << nz - k << endl;
                      m++;
                }}
            }
          }
        count++;
        }

        return 0;
}

