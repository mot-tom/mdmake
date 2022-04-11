#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <math.h>
#include "opt.h"
using namespace std;

#define MAXELEMENT 25
#define BINARY 27

typedef enum {Si, Ge, C, Sn, Mo, S, O, SX, Ga, As, In, Au, Ni, Ti, Al, Mg, Sr, Y, N, Zr, B, Be, Ar, Sb, Te} AtomIonType;

struct info        // information of atom
{
  string name;     // atom name
  double x, y, z;  // coordinate of atom
  int group;       // group number of atom
  string option;   // options of atom
};

struct add         // information of atom
{
  string name;     // atom name
  double x, y, z;  // coordinate of atom
  int point[BINARY];     
};

int swap(int argc, char **argv, optpara opts){

  ostringstream os;
  string empty = "";
  string path = "../";

  int nx = 7, ny = 1, nz = 10;
  if(opts.cellnum.length() > 0) {
    istringstream cell_iss(opts.cellnum.c_str());
    cell_iss >> nx >> ny >> nz;
    if(ny == 0 || nz == 0) {
      cerr << "Number is wrong" << endl;
      return 1;
     }
   }

  const int max_atom  = 300000;       // maximum number of atom
  const int max_group = 1000;         // maximum number of group
  double leng[3] = {0.0, 0.0, 0.0};    // length of each side of crystal  
  double tilt[3][3] = {{0.0, 0.0, 0.0}, 
                       {0.0, 0.0, 0.0},  // x
                       {0.0, 0.0, 0.0}}; // tilt length of each side of crystal

  int atomnum[MAXELEMENT];             // number of atom, respectively
  for(int i=0;i<MAXELEMENT;i++){
    atomnum[i] = 0;
  }


  // file check
  ifstream file_fin(opts.atomcomp.c_str());
  //file_fin.open(opts.atomcomp.c_str());
  if(!file_fin)
  {
    cerr << "Can't open " << opts.atomcomp << "!" << endl;
    return 1;
  }
  cerr << "Read \"" << opts.atomcomp << "\"." << endl;

  string temp_line, file_type;
  int line = 0;
  while(getline(file_fin, temp_line)){
    if(temp_line.find("Gear(") != string::npos && line == 0){file_type = "mdli";}
    if(temp_line.find("OPT(") != string::npos && line == 0){file_type = "mdli";}
    line++;
  }
  file_fin.close();

  if(file_type == "mdli") cerr << "\"" << opts.atomcomp << "\" is input file for Mdlabo." << endl;
  else{
    cerr << "\"" << opts.atomcomp << "\" is wrong file type!" << endl;
    return 1;
  }

  // option 1(select atoms)
  if(opts.method == "random") cerr << "swap method: random" << endl;
  else if(opts.method == "none"){

  }else{
    cerr << "Output type of \"" << opts.method << "\" is wrong!" << endl;
    return 1;
  }

  // option 2(periodic)
  int px = 0, py = 0, pz = 0;
  if(nx == 0){
    cerr << "periodic: nonperiodic"<< endl;
  }
  else if(nx == 1){
    cerr << "periodic: x"<< endl;
    px = 1;
  }
  else if(nx == 2){
    cerr << "periodic: y"<< endl;
    py = 1;
  }
  else if(nx == 3){
    cerr << "periodic: xy"<< endl;
    px = 1;
    py = 1;
  }
  else if(nx == 4){
    cerr << "periodic: z"<< endl;
    pz = 1;
  }
  else if(nx == 5){
    cerr << "periodic: xz"<< endl;
    px = 1;
    pz = 1;
  }
  else if(nx == 6){
    cerr << "periodic: yz"<< endl;
    py = 1;
    pz = 1;
  }
  else if(nx == 7){
    cerr << "periodic: xyz"<< endl;
    px = 1;
    py = 1;
    pz = 1;
  }
  else{
    cerr << "option " << nx << " :This option is unimplemented"<< endl;
    //cerr << "e.g. periodic x=1 y=2 z=4 xyz=7 non=8"<< endl;
    return 1;
  }

  // option 3(end condition)
  if(ny == 1) cerr << "end condition1 :STEP<"<< nz << ">" << endl;
  else{
    cerr << "end condition" << ny << " :This option is unimplemented"<< endl;
    return 1;
  }
  
  // variable(atom)
  string atom_name; //unused
  string atom_type[4]; //unused
  double atom_mass[4]; //unused
  int atom_num[4]; 
  double atom_CIM[4]; //unused
  double atom_CMAS[4]; //unused
  double atom_PIM[4]; //unused
  int all_type;

  int cluster_atom = 0; // type of cluster atoms.

  // read mdlfile
  file_fin.open(opts.atomcomp.c_str());
  int blank = 0, lc = 0;
  string title, first;
  info* zero_atom = new info[max_atom];
  int all_atom = 0;
  
  while(getline(file_fin, temp_line)){
    istringstream line_iss(temp_line.c_str());

    if(blank == 2){
      line_iss >> zero_atom[0].x >> zero_atom[0].y >> zero_atom[0].z;

      if(zero_atom[0].x != 0.0 && lc == 0)
        leng[0] = zero_atom[0].x;
      if(zero_atom[0].y != 0.0 && lc == 0)
        tilt[0][1] = zero_atom[0].y;
      if(zero_atom[0].z != 0.0 && lc == 0)
        tilt[0][2] = zero_atom[0].z;

      if(zero_atom[0].x != 0.0 && lc == 1)
        tilt[1][0] = zero_atom[0].x;
      if(zero_atom[0].y != 0.0 && lc == 1)
        leng[1] = zero_atom[0].y;
      if(zero_atom[0].z != 0.0 && lc == 1)
        tilt[1][2] = zero_atom[0].z;

      if(zero_atom[0].x != 0.0 && lc == 2)
        tilt[2][0] = zero_atom[0].x;
      if(zero_atom[0].y != 0.0 && lc == 2)
        tilt[2][1] = zero_atom[0].y;
      if(zero_atom[0].z != 0.0 && lc == 2)      
        leng[2] = zero_atom[0].z;
      lc++;
    }

    if(blank == 3 && temp_line != ""){
      if(all_atom+1 > max_atom){
        cout << "Dynamic memory space is not enough for atom!" << endl;
        return 1;
      }
      
      line_iss >> first;
      
      if(all_atom == 0){
        if(first != "Si" && first != "SX" && first != "Ge" && first != "Sn" && first != "C"
        && first != "Mo" && first != "S"  && first != "S1" && first != "S2" && first != "O"
        && first != "Ga" && first != "As" && first != "In" && first != "Au" && first != "Ni"
        && first != "Ti" && first != "Al" && first != "Mg" && first != "Sr" && first != "Y"
        && first != "N"  && first != "Zr" && first != "B"  && first != "Be" && first != "Ar"
	&& first != "Sb" && first != "Te") //takizawa
        {
          cout << "Setting of atom name is wrong! (" << first << ")" << endl;
          return 1;
        }
      }
      
      if(first != "V"){
        zero_atom[all_atom].name = first;
        line_iss >> zero_atom[all_atom].x >> zero_atom[all_atom].y >> zero_atom[all_atom].z;
        /*getline(line_iss, zero_atom[all_atom].option);
        if(zero_atom[all_atom].option != ""){
          string str; 
          int group;
          istringstream is(zero_atom[all_atom].option);
          is >> str >> group;
        }*/
      }
      all_atom++;

      if(first == "Si")
        atomnum[Si]++;
      if(first == "SX")
        atomnum[SX]++;
      if(first == "Ge")
        atomnum[Ge]++;
      if(first == "Sn")
        atomnum[Sn]++;
      if(first == "C")
        atomnum[C]++;
      if(first == "Mo")
        atomnum[Mo]++;
      if(first == "S")
        atomnum[S]++;
      if(first == "S1")
        atomnum[S]++;
      if(first == "S2")
        atomnum[S]++;
      if(first == "O")
        atomnum[O]++;
      if(first == "Ga")
        atomnum[Ga]++;
      if(first == "As")
        atomnum[As]++;
      if(first == "In")
        atomnum[In]++;
      if(first == "Au")
        atomnum[Au]++;
      if(first == "Ni")
        atomnum[Ni]++;
      if(first == "Ti")
        atomnum[Ti]++;
      if(first == "Al")
        atomnum[Al]++;
      if(first == "Mg")
        atomnum[Mg]++;
      if(first == "Sr")
        atomnum[Sr]++;
      if(first == "Y")
        atomnum[Y]++;
      if(first == "N")
        atomnum[N]++;
      if(first == "Zr")
        atomnum[Zr]++;
      if(first == "B")
        atomnum[B]++;
      if(first == "Be")
        atomnum[Be]++;
      if(first == "Ar")
        atomnum[Ar]++;
      if(first == "Sb") //takizawa
        atomnum[Sb]++;
      if(first == "Te") //takizawa
        atomnum[Te]++;
    }
    if(temp_line == "")
      blank++;

    if(blank == 1) title = temp_line;
  }
  file_fin.close();

  //atom check
  if(1){
    if(atomnum[Mo] && atomnum[S]){
      all_type = 3;
      atom_name = "MoS2";
      atom_type[0] = "S1"; atom_mass[0] = 32.065; atom_num[0] = atomnum[S];
      atom_type[1] = "S2"; atom_mass[1] = 32.065; atom_num[1] = 0;
      atom_type[2] = "Mo"; atom_mass[2] = 95.94;  atom_num[2] = atomnum[Mo];
    }
    else if(atomnum[Ga] && atomnum[As]){
      all_type = 2;
      atom_name = "GaAs";
      atom_type[0] = "Ga"; atom_mass[0] = 69.723;  atom_num[0] = atomnum[Ga];
      atom_type[1] = "As"; atom_mass[1] = 74.9216; atom_num[1] = atomnum[As];

      atom_CIM[0] = 3.0; atom_CMAS[0] = atom_CIM[0]*0.4725; atom_PIM[0] = atom_CIM[0]*0.6;
      atom_CIM[1] =-3.0; atom_CMAS[1] = atom_CIM[1]*0.4725; atom_PIM[1] = atom_CIM[0]*0.6;
    }
    else if(atomnum[In] && atomnum[As]){
      all_type = 2;
      atom_name = "InAs";
      atom_type[0] = "In"; atom_mass[0] = 114.818; atom_num[0] = atomnum[Ga];
      atom_type[1] = "As"; atom_mass[1] = 74.9216; atom_num[1] = atomnum[As];

      atom_CIM[0] = 3.0; atom_CMAS[0] = atom_CIM[0]*0.4725; atom_PIM[0] = atom_CIM[0]*0.6;
      atom_CIM[1] =-3.0; atom_CMAS[1] = atom_CIM[1]*0.4725; atom_PIM[1] = atom_CIM[0]*0.6;
    }
    else if(atomnum[Si] && atomnum[Ge] && atomnum[Sn]){
      all_type = 3;
      atom_name = "GeSiSn";
      atom_type[0] = "Ge"; atom_mass[0] = 72.63;   atom_num[0] = atomnum[Ge];
      atom_type[1] = "Si"; atom_mass[1] = 28.0855; atom_num[1] = atomnum[Si];
      atom_type[2] = "Sn"; atom_mass[2] = 118.710; atom_num[2] = atomnum[Sn];
    }
    else if(atomnum[Ge] && atomnum[Sb] && atomnum[Te]){ //takizawa
      all_type = 3;
      atom_name = "GeSbTe";
      atom_type[0] = "Ge"; atom_mass[0] = 72.63;   atom_num[0] = atomnum[Ge];
      atom_type[1] = "Sb"; atom_mass[1] = 121.76; atom_num[1] = atomnum[Sb];
      atom_type[2] = "Te"; atom_mass[2] = 127.60; atom_num[2] = atomnum[Te];
    }
    else if(atomnum[Si] && atomnum[SX] && atomnum[O]){
      all_type = 3;
      atom_name = "SiO2";
      atom_type[0] = "SX"; atom_mass[0] = 255;     atom_num[0] = atomnum[SX];
      atom_type[1] = "Si"; atom_mass[1] = 28.0855; atom_num[1] = atomnum[Si];
      atom_type[2] = "O";  atom_mass[2] = 15.9994; atom_num[2] = atomnum[O];
    }
    else if(atomnum[Si] && atomnum[Ge]){
      all_type = 2;
      atom_name = "SiGe";
      atom_type[0] = "Si"; atom_mass[0] = 28.0855; atom_num[0] = atomnum[Si];
      atom_type[1] = "Ge"; atom_mass[1] = 72.63;   atom_num[1] = atomnum[Ge];
    }
    else if(atomnum[Si] && atomnum[Sn]){
      all_type = 2;
      atom_name = "SiSn";
      atom_type[0] = "Si"; atom_mass[0] = 28.0855; atom_num[0] = atomnum[Si];
      atom_type[1] = "Sn"; atom_mass[1] = 118.710; atom_num[1] = atomnum[Sn];
    }
    else if(atomnum[Si] && atomnum[C]){
      all_type = 2;
      atom_name = "SiC";
      atom_type[0] = "Si"; atom_mass[0] = 28.0855; atom_num[0] = atomnum[Si];
      atom_type[1] = "C";  atom_mass[1] = 12.0107; atom_num[1] = atomnum[C];
    }
    else if(atomnum[Si] && atomnum[N]){
      all_type = 2;
      atom_name = "Si3N4";
      atom_type[0] = "Si"; atom_mass[0] = 28.0855; atom_num[0] = atomnum[Si];
      atom_type[1] = "N";  atom_mass[1] = 14.0067; atom_num[1] = atomnum[N];

      atom_CIM[0] = 4.0; atom_CMAS[0] = atom_CIM[0]*0.4725; atom_PIM[0] = atom_CIM[0]*0.6;
      atom_CIM[1] =-3.0; atom_CMAS[1] = atom_CIM[1]*0.4725; atom_PIM[1] = atom_CIM[0]*0.6;
    }
    else if(atomnum[Si] && atomnum[O] && atomnum[Ar]){
      all_type = 3;
      atom_name = "SiO2Ar";
      atom_type[0] = "Si"; atom_mass[0] = 28.0855; atom_num[0] = atomnum[Si];
      atom_type[1] = "O";  atom_mass[1] = 15.9994; atom_num[1] = atomnum[O];
      atom_type[2] = "Ar"; atom_mass[2] = 39.948;  atom_num[2] = atomnum[Ar];

      atom_CIM[0] = 4.0; atom_CMAS[0] = atom_CIM[0]*0.4725; atom_PIM[0] = atom_CIM[0]*0.6;
      atom_CIM[1] =-2.0; atom_CMAS[1] = atom_CIM[1]*0.4725; atom_PIM[1] = atom_CIM[1]*0.6;
      atom_CIM[2] = 0.0; atom_CMAS[2] = atom_CIM[2]*0.4725; atom_PIM[2] = atom_CIM[2]*0.6;
    }
    else if(atomnum[Si] && atomnum[Ar]){
      all_type = 2;
      atom_name = "SiAr";
      atom_type[0] = "Si"; atom_mass[0] = 28.0855; atom_num[0] = atomnum[Si];
      atom_type[1] = "Ar"; atom_mass[1] = 39.948;  atom_num[1] = atomnum[Ar];

      atom_CIM[0] = 4.0; atom_CMAS[0] = atom_CIM[0]*0.4725; atom_PIM[0] = atom_CIM[0]*0.6;
      atom_CIM[1] = 0.0; atom_CMAS[1] = atom_CIM[1]*0.4725; atom_PIM[1] = atom_CIM[1]*0.6;
    }
    else if(atomnum[Si] && atomnum[O]){
      all_type = 2;
      atom_name = "SiO2";
      atom_type[0] = "Si"; atom_mass[0] = 28.0855; atom_num[0] = atomnum[Si];
      atom_type[1] = "O";  atom_mass[1] = 15.9994; atom_num[1] = atomnum[O];

      atom_CIM[0] = 4.0; atom_CMAS[0] = atom_CIM[0]*0.4725; atom_PIM[0] = atom_CIM[0]*0.6;
      atom_CIM[1] =-2.0; atom_CMAS[1] = atom_CIM[1]*0.4725; atom_PIM[1] = atom_CIM[0]*0.6;
    }
    else if(atomnum[Ge] && atomnum[O]){
      all_type = 2;
      atom_name = "GeO2";
      atom_type[0] = "Ge"; atom_mass[0] = 72.63;   atom_num[0] = atomnum[Ge];
      atom_type[1] = "O";  atom_mass[1] = 15.9994; atom_num[1] = atomnum[O];

      atom_CIM[0] = 4.0; atom_CMAS[0] = atom_CIM[0]*0.4725; atom_PIM[0] = atom_CIM[0]*0.6;
      atom_CIM[1] =-2.0; atom_CMAS[1] = atom_CIM[1]*0.4725; atom_PIM[1] = atom_CIM[1]*0.6;
    }
    else if(atomnum[Sn] && atomnum[O]){
      all_type = 2;
      atom_name = "SnO2";
      atom_type[0] = "Sn"; atom_mass[0] = 118.710; atom_num[0] = atomnum[Sn];
      atom_type[1] = "O";  atom_mass[1] = 15.9994; atom_num[1] = atomnum[O];

      atom_CIM[0] = 4.0; atom_CMAS[0] = atom_CIM[0]*0.4725; atom_PIM[0] = atom_CIM[0]*0.6;
      atom_CIM[1] =-2.0; atom_CMAS[1] = atom_CIM[1]*0.4725; atom_PIM[1] = atom_CIM[1]*0.6;
    }
    else if(atomnum[Ti] && atomnum[O]){
      all_type = 2;
      atom_name = "TiO2";
      atom_type[0] = "Ti"; atom_mass[0] = 47.867;  atom_num[0] = atomnum[Ti];
      atom_type[1] = "O";  atom_mass[1] = 15.9994; atom_num[1] = atomnum[O];

      atom_CIM[0] = 4.0; atom_CMAS[0] = atom_CIM[0]*0.4725; atom_PIM[0] = atom_CIM[0]*0.6;
      atom_CIM[1] =-2.0; atom_CMAS[1] = atom_CIM[1]*0.4725; atom_PIM[1] = atom_CIM[1]*0.6;
    }
    else if(atomnum[Zr] && atomnum[O]){
      all_type = 2;
      atom_name = "ZrO2";
      atom_type[0] = "Zr"; atom_mass[0] = 91.224;  atom_num[0] = atomnum[Zr];
      atom_type[1] = "O";  atom_mass[1] = 15.9994; atom_num[1] = atomnum[O];

      atom_CIM[0] = 4.0; atom_CMAS[0] = atom_CIM[0]*0.4725; atom_PIM[0] = atom_CIM[0]*0.6;
      atom_CIM[1] =-2.0; atom_CMAS[1] = atom_CIM[1]*0.4725; atom_PIM[1] = atom_CIM[1]*0.6;
    }
    else if(atomnum[Mg] && atomnum[O]){
      all_type = 2;
      atom_name = "MgO";
      atom_type[0] = "Mg"; atom_mass[0] = 24.3050; atom_num[0] = atomnum[Mg];
      atom_type[1] = "O";  atom_mass[1] = 15.9994; atom_num[1] = atomnum[O];

      atom_CIM[0] = 2.0; atom_CMAS[0] = atom_CIM[0]*0.4725; atom_PIM[0] = atom_CIM[0]*0.6;
      atom_CIM[1] =-2.0; atom_CMAS[1] = atom_CIM[1]*0.4725; atom_PIM[1] = atom_CIM[0]*0.6;
    }
    else if(atomnum[Sr] && atomnum[O]){
      all_type = 2;
      atom_name = "SrO";
      atom_type[0] = "Sr"; atom_mass[0] = 87.62;   atom_num[0] = atomnum[Mg];
      atom_type[1] = "O";  atom_mass[1] = 15.9994; atom_num[1] = atomnum[O];

      atom_CIM[0] = 2.0; atom_CMAS[0] = atom_CIM[0]*0.4725; atom_PIM[0] = atom_CIM[0]*0.6;
      atom_CIM[1] =-2.0; atom_CMAS[1] = atom_CIM[1]*0.4725; atom_PIM[1] = atom_CIM[1]*0.6;
    }
    else if(atomnum[Be] && atomnum[O]){
      all_type = 2;
      atom_name = "BeO";
      atom_type[0] = "Be"; atom_mass[0] = 9.01218; atom_num[0] = atomnum[Be];
      atom_type[1] = "O";  atom_mass[1] = 15.9994; atom_num[1] = atomnum[O];

      atom_CIM[0] = 2.0; atom_CMAS[0] = atom_CIM[0]*0.4725; atom_PIM[0] = atom_CIM[0]*0.6;
      atom_CIM[1] =-2.0; atom_CMAS[1] = atom_CIM[1]*0.4725; atom_PIM[1] = atom_CIM[1]*0.6;
    }
    else if(atomnum[B] && atomnum[N]){
      all_type = 2;
      atom_name = "BN";
      atom_type[0] = "B";  atom_mass[0] = 10.811;  atom_num[0] = atomnum[B];
      atom_type[1] = "N";  atom_mass[1] = 14.0067; atom_num[1] = atomnum[N];

      atom_CIM[0] = 3.0; atom_CMAS[0] = atom_CIM[0]*0.4725; atom_PIM[0] = atom_CIM[0]*0.6;
      atom_CIM[1] =-3.0; atom_CMAS[1] = atom_CIM[1]*0.4725; atom_PIM[1] = atom_CIM[1]*0.6;
    }
    else if(atomnum[Al] && atomnum[N]){
      all_type = 2;
      atom_name = "AlN";
      atom_type[0] = "Al"; atom_mass[0] = 26.9815; atom_num[0] = atomnum[Al];
      atom_type[1] = "N";  atom_mass[1] = 14.0067; atom_num[1] = atomnum[N];
  
      atom_CIM[0] = 3.0; atom_CMAS[0] = atom_CIM[0]*0.4725; atom_PIM[0] = atom_CIM[0]*0.6;
      atom_CIM[1] =-3.0; atom_CMAS[1] = atom_CIM[1]*0.4725; atom_PIM[1] = atom_CIM[1]*0.6;
    }
    else if(atomnum[Al] && atomnum[O]){
      all_type = 2;
      atom_name = "Al2O3";
      atom_type[0] = "Al"; atom_mass[0] = 26.9815; atom_num[0] = atomnum[Al];
      atom_type[1] = "O";  atom_mass[1] = 15.9994; atom_num[1] = atomnum[O];

      atom_CIM[0] = 3.0; atom_CMAS[0] = atom_CIM[0]*0.4725; atom_PIM[0] = atom_CIM[0]*0.6;
      atom_CIM[1] =-2.0; atom_CMAS[1] = atom_CIM[1]*0.4725; atom_PIM[1] = atom_CIM[1]*0.6;
    }
    else if(atomnum[Y] && atomnum[O]){
      all_type = 2;
      atom_name = "Y2O3";
      atom_type[0] = "Y";  atom_mass[0] = 88.9059; atom_num[0] = atomnum[Y];
      atom_type[1] = "O";  atom_mass[1] = 15.9994; atom_num[1] = atomnum[O];

      atom_CIM[0] = 3.0; atom_CMAS[0] = atom_CIM[0]*0.4725; atom_PIM[0] = atom_CIM[0]*0.6;
      atom_CIM[1] =-2.0; atom_CMAS[1] = atom_CIM[1]*0.4725; atom_PIM[1] = atom_CIM[1]*0.6;
    }
    else if(atomnum[Si]){
      all_type = 1;
      atom_name = "Si";
      atom_type[0] = atom_name; atom_mass[0] = 28.0855; atom_num[0] = atomnum[Si];
    }
    else if(atomnum[Ge]){
      all_type = 1;
      atom_name = "Ge";
      atom_type[0] = atom_name; atom_mass[0] = 72.63; atom_num[0] = atomnum[Ge];
    }
    else if(atomnum[Sn]){
      all_type = 1;
      atom_name = "Sn";
      atom_type[0] = atom_name; atom_mass[0] = 118.710; atom_num[0] = atomnum[Sn];
    }
    else if(atomnum[C]){
      all_type = 1;
      atom_name = "C";
      atom_type[0] = atom_name; atom_mass[0] = 12.0107; atom_num[0] = atomnum[C];
    }
    else if(atomnum[Au]){
      all_type = 1;
      atom_name = "Au";
      atom_type[0] = atom_name; atom_mass[0] = 196.9666; atom_num[0] = atomnum[Au];
    }
    else if(atomnum[Ni]){
      all_type = 1;
      atom_name = "Ni";
      atom_type[0] = atom_name; atom_mass[0] = 58.6934; atom_num[0] = atomnum[Au];
    }

    if(all_type == 1){
      cerr << "This program is not support the kind of single atom!" << endl;
      return 1;
    }

    // notice
    if(atom_name != "SiGe"){
      cerr << "This program is only support SiGe now!" << endl;
      return 1;
    }
    cluster_atom = 1;
  }
  
  //sort_atom
  int kindone = 0; //Si
  int kindtwo = 1; //Ge
  add atom_one[all_atom];
  cerr << "all:" << all_atom << " swap:" << atom_num[kindtwo] << endl;
  int k_one = 0, k_two = atom_num[kindone];
  for(int i=0; i<all_atom; i++){
    if(zero_atom[i].name == atom_type[kindone]){
      if(k_one >= atom_num[kindone]){
        cerr << "program error:1-1" << endl;
        return 1;
      }
      atom_one[k_one].x    = zero_atom[i].x;
      atom_one[k_one].y    = zero_atom[i].y;
      atom_one[k_one].z    = zero_atom[i].z;
      atom_one[k_one].name = zero_atom[i].name;
      for(int j=0; j<BINARY; j++) atom_one[k_one].point[j] = 0;
      ++k_one;
    }else if(zero_atom[i].name == atom_type[kindtwo]){
      if(k_two >= atom_num[kindone] + atom_num[kindtwo]){
        cerr << "program error:1-2" << endl;
        return 1;
      }
      atom_one[k_two].x    = zero_atom[i].x;
      atom_one[k_two].y    = zero_atom[i].y;
      atom_one[k_two].z    = zero_atom[i].z;
      atom_one[k_two].name = zero_atom[i].name;
      for(int j=0; j<BINARY; j++) atom_one[k_two].point[j] = 0;    
      ++k_two;
    }
  }
  delete[] zero_atom;

  int cut = pow(3.0, 2);
  int ok;
  double dx, dy, dz, r, rx, ry, rz, rxy, ryz, rxz, rxyz;
  for(int i=0; i<all_atom; i++){
    for(int j=i+1; j<all_atom; j++){
       ok = 0;
       dx = fabs(atom_one[i].x - atom_one[j].x) * leng[0];
       dy = fabs(atom_one[i].y - atom_one[j].y) * leng[1];
       dz = fabs(atom_one[i].z - atom_one[j].z) * leng[2];

       r    = pow(dx, 2.0) + pow(dy, 2.0) + pow(dz, 2.0);
       rx   = pow(leng[0] - dx, 2.0) + pow(dy, 2.0) + pow(dz, 2.0);
       ry   = pow(dx, 2.0) + pow(leng[1] - dy, 2.0) + pow(dz, 2.0);
       rz   = pow(dx, 2.0) + pow(dy, 2.0) + pow(leng[2] - dz, 2.0);
       rxy  = pow(leng[0] - dx, 2.0) + pow(leng[1] - dy, 2.0) + pow(dz, 2.0);
       ryz  = pow(dx, 2.0) + pow(leng[1] - dy, 2.0) + pow(leng[2] - dz, 2.0);
       rxz  = pow(leng[0] - dx, 2.0) + pow(dy, 2.0) + pow(leng[2] - dz, 2.0);
       rxyz = pow(leng[0] - dx, 2.0) + pow(leng[1] - dy, 2.0) + pow(leng[2] - dz, 2.0);

       if(dx <= leng[0]/2.0 && dy <= leng[1]/2.0 && dz <= leng[2]/2.0 && r    <cut) ok = 1;
       if(dx <= leng[0]/2.0 && dy <= leng[1]/2.0 && dz >  leng[2]/2.0 && rz   <cut && pz==1) ok = 1;
       if(dx <= leng[0]/2.0 && dy >  leng[1]/2.0 && dz <= leng[2]/2.0 && ry   <cut && py==1) ok = 1;
       if(dx <= leng[0]/2.0 && dy >  leng[1]/2.0 && dz >  leng[2]/2.0 && ryz  <cut && py==1 && pz==1) ok = 1;
       if(dx >  leng[0]/2.0 && dy <= leng[1]/2.0 && dz <= leng[2]/2.0 && rx   <cut && px==1) ok = 1;
       if(dx >  leng[0]/2.0 && dy <= leng[1]/2.0 && dz >  leng[2]/2.0 && rxz  <cut && px==1 && pz==1) ok = 1;
       if(dx >  leng[0]/2.0 && dy >  leng[1]/2.0 && dz <= leng[2]/2.0 && rxy  <cut && px==1 && py==1) ok = 1;
       if(dx >  leng[0]/2.0 && dy >  leng[1]/2.0 && dz >  leng[2]/2.0 && rxyz <cut && px==1 && py==1 && pz==1) ok = 1;
       
       if(ok == 1){
         ++atom_one[i].point[0];
         ++atom_one[j].point[0];
         atom_one[i].point[atom_one[i].point[0]] = j;
         atom_one[j].point[atom_one[j].point[0]] = i;
       }

    }
  }


  cerr << endl << "STEP:" << endl;
  int STEP = 0;
  int swap1, swap2;
  int ratio1, ratio2;
  while(STEP<nz){
    //swap atoms
    if(opts.method == "random"){
      ratio1 = 0;
      ratio2 = 0;
      srand((unsigned int)time(NULL));
      swap1 =  rand() % atom_num[kindone];                       //Si
      swap2 = (rand() % atom_num[kindtwo]) + atom_num[kindone];  //Ge
      for(int i=1; i<=atom_one[swap1].point[0]; i++){
        if(atom_one[swap1].point[i] >= atom_num[kindone]) ratio1 = ratio1 + 1;
      }
      for(int i=1; i<=atom_one[swap2].point[0]; i++){
        if(atom_one[swap2].point[i] >= atom_num[kindone]) ratio2 = ratio2 + 1;
      }
      if(ratio1 > ratio2){
        swap(atom_one[swap1].x, atom_one[swap2].x);
        swap(atom_one[swap1].y, atom_one[swap2].y);
        swap(atom_one[swap1].z, atom_one[swap2].z);
        for(int i=1; i<=atom_one[swap1].point[0]; i++){
          for(int j=1; j<=atom_one[atom_one[swap1].point[i]].point[0]; j++){
            if(atom_one[atom_one[swap1].point[i]].point[j] = swap1){
              atom_one[atom_one[swap1].point[i]].point[j] = swap2;
            }
          }
        }
        for(int i=1; i<=atom_one[swap2].point[0]; i++){
          for(int j=1; j<=atom_one[atom_one[swap2].point[i]].point[0]; j++){
            if(atom_one[atom_one[swap2].point[i]].point[j] = swap2){
              atom_one[atom_one[swap2].point[i]].point[j] = swap1;
            }
          }
        }
        swap(atom_one[swap1].point, atom_one[swap2].point);
        ++STEP;
        if(STEP<10) cerr << "  " << STEP;
	if(10<=STEP) cerr << " " << STEP;
	if(STEP%10==0 && STEP!=nz) cerr << endl;

      }
    }else if(opts.method == "none"){ //takizawa

    }
  }
  cerr << " Finish!" << endl << endl;

  //output
  cout << "Gear(5) RAND=10 " << opts.potential << " T=300.0 STEP=(0,1000100)" << endl;
  cout << "BLOCK FIXANGLE W=INF Q=INF LOGSTEP=10 PRINTVELOCITY dt=0.1fs scratch(100)" << endl << endl;
  cout << title << endl << endl;
  cout << leng[0] << " " << tilt[0][1] << " " << tilt[0][2] << opts.pbc_x << endl;
  cout << tilt[1][0] << " " << leng[1] << " " << tilt[1][2] << opts.pbc_y << endl;
  cout << tilt[2][0] << " " << tilt[2][1] << " " << leng[2] << opts.pbc_z << endl << endl;
  for(int i=0;i<all_atom;i++){
    cout << atom_one[i].name << " " << atom_one[i].x << " " << atom_one[i].y << " " << atom_one[i].z << endl;
  }
  cout << endl;

  return 0;
}

