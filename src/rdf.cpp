//L652
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

#define MAXELEMENT 27 // takizawa
typedef enum {Si, Ge, C, Sn, Mo, S, O, SX, Ga, As, In, Au, Ni, Ti, Al, Mg, Sr, Y, N, Zr, B, Be, Ar, Sb, Te, Na, Cl} AtomIonType; //takizawa

struct info        // information of atom
{
  string name;     // atom name
  double x, y, z;  // coordinate of atom
  int lay;         // layer number of atom
  string option;   // options of atom
};

extern int shape(string shape, double x, double y, double z);

double sunsq(double x, double y, double z){ //takizawa
  return sqrt(pow(x, 2.0) + pow(y, 2.0) + pow(z, 2.0));
}

int rdf(int argc, char **argv, optpara opts){
  const int max_atom = 300000;         // maximum number of atom

  int atomA = 1, atomB = 1; 
  double cutoffC = 4.0;

  if(opts.cellnum.length() > 0) {
    istringstream cell_iss(opts.cellnum.c_str());
    cell_iss >> atomA >> atomB >> cutoffC;
    if(atomA == 0 || atomB == 0 || cutoffC == 0) {
      cerr << "Cell number is wrong" << endl;
      return 1;
     }
   }

  cerr << "cell shape " << opts.shape << endl;
  if(opts.shape != "cuboid"){
    cerr << "Sorry! This program only corresponds 'cuboid'" << endl;
    return 1;
  }

  double leng[3] = {0.0, 0.0, 0.0};    // length of each side of crystal  
  double tilt[3][3] = {{0.0, 0.0, 0.0}, 
                       {0.0, 0.0, 0.0},  // x
                       {0.0, 0.0, 0.0}}; // tilt length of each side of crystal

  int atomnum[MAXELEMENT];             // number of atom, respectively
  for(int i=0;i<MAXELEMENT;i++){
    atomnum[i] = 0;
  }

  ifstream file_fin(opts.atomcomp.c_str());
  
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
    cerr << "input type of \"" << opts.atomcomp << "\" is wrong for RDF!" << endl;
    return 1;  
  }

  //takizawa
  //bond type
  if(opts.atomratio == "rdf"){
    cerr << "Mode: RDF" << endl;
  }else if(opts.atomratio == "num"){
    cerr << "Mode: Coordinate Number" << endl;
  }else{
    cerr << "option of \"" << opts.atomratio << "\" is wrong!" << endl;
    return 1;
  }

  file_fin.open(opts.atomcomp.c_str());
  int blank = 0, lc = 0;
  string title, first;
  info* zero_atom = new info[max_atom];
  int all_atom = 0;
  if(file_type == "mdli"){
    while(getline(file_fin, temp_line)){
      istringstream line_iss(temp_line.c_str());

      if(blank == 2)
      {
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
          tilt[1][2] = zero_atom[0].z;  double atom_CIM[4];
  double atom_CMAS[4];
  double atom_PIM[4];

        if(zero_atom[0].x != 0.0 && lc == 2)
          tilt[2][0] = zero_atom[0].x;
        if(zero_atom[0].y != 0.0 && lc == 2)
          tilt[2][1] = zero_atom[0].y;
        if(zero_atom[0].z != 0.0 && lc == 2)      
          leng[2] = zero_atom[0].z;

        if(lc < 3){
          cerr << "|ul" << lc << "0 ul" << lc << "1 ul" << lc << "2|";
          if(lc == 1) cerr << "=";
          else cerr << " ";
          if(lc == 0)
            cerr << showpoint << "|" << setw(6) << leng[0]    << " " << tilt[0][1] << " " << tilt[0][2] << "|";
          else if(lc == 1)
            cerr << showpoint << "|" << setw(6) << tilt[1][0] << " " << leng[1]    << " " << tilt[1][2] << "|";
          else if(lc == 2)
            cerr << showpoint << "|" << setw(6) << tilt[2][0] << " " << tilt[2][1] << " " << leng[2] << "|";
          cerr << endl;
         }
        lc++;
      }

      if(blank == 3 && temp_line != "")
     {
        if(all_atom+1 > max_atom)
        {
          cerr << "Dynamic memory space is not enough for atom!" << endl;
          return 1;
        }
      
        line_iss >> first;
      
        if(all_atom == 0)
        {
          if(first != "Si" && first != "SX" && first != "Ge" && first != "Sn" && first != "C"
          && first != "Mo" && first != "S"  && first != "S1" && first != "S2" && first != "O"
          && first != "Ga" && first != "As" && first != "In" && first != "Au" && first != "Ni"
          && first != "Ti" && first != "Al" && first != "Mg" && first != "Sr" && first != "Y"
          && first != "N"  && first != "Zr" && first != "B"  && first != "Be" && first != "Ar"
	  && first != "Sb" && first != "Te" && first != "Na" && first != "Cl") // takizawa
          {
            cerr << "Setting of atom name is wrong! (" << first << ")" << endl;
            return 1;
          }
        }
      
        if(first != "V")
        {
          zero_atom[all_atom].name = first;
          line_iss >> zero_atom[all_atom].x >> zero_atom[all_atom].y >> zero_atom[all_atom].z;
          getline(line_iss, zero_atom[all_atom].option);
          if(zero_atom[all_atom].option != ""){
            string str; int lay;
            istringstream is(zero_atom[all_atom].option);
            is >> str >> lay;
            if(str=="F") zero_atom[all_atom].lay = 255;
            if(str=="G") zero_atom[all_atom].lay = lay;
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
          if(first == "Na") //takizawa
            atomnum[Na]++;
          if(first == "Cl") //takizawa
            atomnum[Cl]++; 
        }
      }

      if(temp_line == "")
        blank++;

      if(blank == 1) title = temp_line;
    }
  }else{
        return 1;  
  }
  file_fin.close();

  /* takizawa
  leng[0] *= atomA;
  leng[1] *= atomB;
  leng[2] *= cutoffC;
  tilt[0][1] *= atomA;
  tilt[0][2] *= atomA;
  tilt[1][0] *= atomB;
  tilt[1][2] *= atomB;
  tilt[2][0] *= cutoffC;
  tilt[2][1] *= cutoffC;
  */

  string atom_name;
  string atom_type[4];
  double atom_mass[4];
  int atom_num[4];
  double atom_CIM[4];
  double atom_CMAS[4];
  double atom_PIM[4];
  int all_type;
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
  else if(atomnum[Na] && atomnum[Cl] && atomnum[Si] && atomnum[O]){ //takizawa
    all_type = 4;
    atom_name = "NaCl-SiO2";
    atom_type[0] = "Na"; atom_mass[0] = 22.9897; atom_num[0] = atomnum[Na];
    atom_type[1] = "Cl"; atom_mass[1] = 35.446;  atom_num[1] = atomnum[Cl];
    atom_type[2] = "Si"; atom_mass[2] = 28.0855; atom_num[2] = atomnum[Si];
    atom_type[3] = "O";  atom_mass[3] = 15.9994; atom_num[3] = atomnum[O];

    atom_CIM[0] = 1.0; atom_CMAS[0] = atom_CIM[0]*0.4725; atom_PIM[0] = atom_CIM[0]*0.6;
    atom_CIM[1] =-1.0; atom_CMAS[1] = atom_CIM[1]*0.4725; atom_PIM[1] = atom_CIM[1]*0.6;
    atom_CIM[2] = 4.0; atom_CMAS[2] = atom_CIM[2]*0.4725; atom_PIM[2] = atom_CIM[2]*0.6;
    atom_CIM[3] =-2.0; atom_CMAS[3] = atom_CIM[3]*0.4725; atom_PIM[3] = atom_CIM[3]*0.6;
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
  else if(atomnum[Al] && atomnum[Si] && atomnum[O]){
    all_type = 3;
    atom_name = "Al2O3-SiO2";
    atom_type[0] = "Al"; atom_mass[0] = 26.9815; atom_num[0] = atomnum[Al];
    atom_type[1] = "Si"; atom_mass[1] = 28.0855; atom_num[1] = atomnum[Si];
    atom_type[2] = "O";  atom_mass[2] = 15.9994; atom_num[2] = atomnum[O];
  
    atom_CIM[0] = 3.0; atom_CMAS[0] = atom_CIM[0]*0.4725; atom_PIM[0] = atom_CIM[0]*0.6;
    atom_CIM[1] = 4.0; atom_CMAS[1] = atom_CIM[1]*0.4725; atom_PIM[1] = atom_CIM[1]*0.6;
    atom_CIM[2] =-2.0; atom_CMAS[2] = atom_CIM[2]*0.4725; atom_PIM[2] = atom_CIM[1]*0.6;
  }
  else if(atomnum[Al] && atomnum[O]){
    all_type = 2;
    atom_name = "Al2O3";
    atom_type[0] = "Al"; atom_mass[0] = 26.9815; atom_num[0] = atomnum[Al];
    atom_type[1] = "O";  atom_mass[1] = 15.9994; atom_num[1] = atomnum[O];

    atom_CIM[0] = 3.0; atom_CMAS[0] = atom_CIM[0]*0.4725; atom_PIM[0] = atom_CIM[0]*0.6;
    atom_CIM[1] =-2.0; atom_CMAS[1] = atom_CIM[1]*0.4725; atom_PIM[1] = atom_CIM[1]*0.6;
  }
  else if(atomnum[Si] && atomnum[O]){
    all_type = 2;
    atom_name = "SiO2";
    atom_type[0] = "Si"; atom_mass[0] = 28.0855; atom_num[0] = atomnum[Si];
    atom_type[1] = "O";  atom_mass[1] = 15.9994; atom_num[1] = atomnum[O];

    atom_CIM[0] = 4.0; atom_CMAS[0] = atom_CIM[0]*0.4725; atom_PIM[0] = atom_CIM[0]*0.6;
    atom_CIM[1] =-2.0; atom_CMAS[1] = atom_CIM[1]*0.4725; atom_PIM[1] = atom_CIM[0]*0.6;
  }
  else if(atomnum[Y] && atomnum[O]){
    all_type = 2;
    atom_name = "Y2O3";
    atom_type[0] = "Y";  atom_mass[0] = 88.9059; atom_num[0] = atomnum[Y];
    atom_type[1] = "O";  atom_mass[1] = 15.9994; atom_num[1] = atomnum[O];

    atom_CIM[0] = 3.0; atom_CMAS[0] = atom_CIM[0]*0.4725; atom_PIM[0] = atom_CIM[0]*0.6;
    atom_CIM[1] =-2.0; atom_CMAS[1] = atom_CIM[1]*0.4725; atom_PIM[1] = atom_CIM[1]*0.6;
  }
  else if(atomnum[Na] && atomnum[Cl]){
    all_type = 2;
    atom_name = "NaCl";
    atom_type[0] = "Na";  atom_mass[0] = 22.9897; atom_num[0] = atomnum[Y];
    atom_type[1] = "Cl";  atom_mass[1] = 35.446; atom_num[1] = atomnum[O];

    atom_CIM[0] = 1.0; atom_CMAS[0] = atom_CIM[0]*0.4725; atom_PIM[0] = atom_CIM[0]*0.6;
    atom_CIM[1] =-1.0; atom_CMAS[1] = atom_CIM[1]*0.4725; atom_PIM[1] = atom_CIM[1]*0.6;
  }
  else if(atomnum[Si]){
    all_type = 1;
    atom_name = "Si"; atom_mass[0] = 28.0855; atom_num[0] = atomnum[Si];
    atom_type[0] = atom_name;
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
  


  //make under part by takizawa

  //cutoff(-n)

  if(all_type+1 < atomA){
    cerr << " Err: Atom type " << atomA << " is unavailable!" << endl;
    return 1;
  }
  if(all_type+1 < atomB){
    cerr << " Err: Atom type " << atomB << " is unavailable!" << endl;
    return 1;
  }

  string center;
  if(all_type+1 == atomA){
   center = "all";
  }else{
   center = atom_type[atomA - 1];
  }

  string around;
  if(all_type+1 == atomB){
   around = "all";
  }else{
   around = atom_type[atomB - 1];
  }

  double cutDEF = cutoffC;

  cerr << "Center atom    : " << center << endl;
  cerr << "Around atom    : " << around << endl;
  cerr << "Cutoff Distance: " << cutDEF << " [Angstrom]" << endl;

  
  //cell error(cuboid)
  if(leng[0] / 2 < cutDEF || leng[1] / 2 < cutDEF || leng[2] / 2 < cutDEF ){
    cerr << " Err: A unit cell is very small!" << endl;
    return 1;
  }

  //max (cuboid)
  double max[3] = {0.0, 0.0, 0.0};
  double max2[3] = {0.0, 0.0, 0.0};
  double min[3] = {leng[0], leng[1], leng[2]};
  double min2[3] = {leng[0], leng[1], leng[2]};

  if(center == "all" || around == "all"){
    for(int i=0;i<all_atom;i++){
      if(max[0]<zero_atom[i].x * leng[0])  max[0] = zero_atom[i].x * leng[0];
      if(max[1]<zero_atom[i].y * leng[1])  max[1] = zero_atom[i].y * leng[1];
      if(max[2]<zero_atom[i].z * leng[2])  max[2] = zero_atom[i].z * leng[2];
      if(zero_atom[i].x * leng[0]<min[0])  min[0] = zero_atom[i].x * leng[0];
      if(zero_atom[i].y * leng[1]<min[1])  min[1] = zero_atom[i].y * leng[1];
      if(zero_atom[i].z * leng[2]<min[2])  min[2] = zero_atom[i].z * leng[2];
    }
  }else{
    for(int i=0;i<all_atom;i++){
      if(zero_atom[i].name == center){
        if(max[0]<zero_atom[i].x * leng[0])  max[0] = zero_atom[i].x * leng[0];
        if(max[1]<zero_atom[i].y * leng[1])  max[1] = zero_atom[i].y * leng[1];
        if(max[2]<zero_atom[i].z * leng[2])  max[2] = zero_atom[i].z * leng[2];
        if(zero_atom[i].x * leng[0]<min[0])  min[0] = zero_atom[i].x * leng[0];
        if(zero_atom[i].y * leng[1]<min[1])  min[1] = zero_atom[i].y * leng[1];
        if(zero_atom[i].z * leng[2]<min[2])  min[2] = zero_atom[i].z * leng[2];
      }
      if(zero_atom[i].name == around && center != around){
        if(max2[0]<zero_atom[i].x * leng[0]) max2[0] = zero_atom[i].x * leng[0];
        if(max2[1]<zero_atom[i].y * leng[1]) max2[1] = zero_atom[i].y * leng[1];
        if(max2[2]<zero_atom[i].z * leng[2]) max2[2] = zero_atom[i].z * leng[2];
        if(zero_atom[i].x * leng[0]<min2[0]) min2[0] = zero_atom[i].x * leng[0];
        if(zero_atom[i].y * leng[1]<min2[1]) min2[1] = zero_atom[i].y * leng[1];
        if(zero_atom[i].z * leng[2]<min2[2]) min2[2] = zero_atom[i].z * leng[2];
      }
    }
    for(int i=0;i<3;i++){
      if(max[i] < max2[i]) max[i] = max2[i];
      if(min[i] > min2[i]) min[i] = min2[i];
    }
  }

  double INTERFACE = 66.87;
  min[2] = INTERFACE - 2.5;
  max[2] = INTERFACE + 2.5;

  //boundary (cuboid)
  const string p  = "p";
  string bx = "p";
  string by = "p";
  string bz = "p";

  //binning
  double dr = 0.05; // input
  int ratio = 0; //0 ~ cutDEF/dr
  int ratiomax = ceil(cutDEF / dr);
  int coord[ratiomax];
  int Rate;
  for(int i=0;i<ratiomax;i++) coord[i] = 0;

  //coordination
  const int neighbormax = 20;
  int neighbor[neighbormax];
  for(int i=0;i<neighbormax;i++) neighbor[i] = 0;

  //distance
  double dx, dy, dz;
  double r, rx, ry, rz, rxy, ryz, rxz, rxyz;

  //record
  int centercount = 0;
  int aroundcount = 0; 

  for(int i=0;i<all_atom;i++){
    if(zero_atom[i].name != center && center != "all") continue;
    if((zero_atom[i].x * leng[0] - min[0] < cutDEF || max[0] - zero_atom[i].x * leng[0] < cutDEF) && bx != p) continue; //surface effect(cuboid)
    if((zero_atom[i].y * leng[1] - min[1] < cutDEF || max[1] - zero_atom[i].y * leng[1] < cutDEF) && by != p) continue;
    if((zero_atom[i].z * leng[2] - min[2] < cutDEF || max[2] - zero_atom[i].z * leng[2] < cutDEF) && bz != p) continue;
    if((zero_atom[i].x * leng[0] < min[0] || max[0] < zero_atom[i].x * leng[0]) && bx == p) continue;
    if((zero_atom[i].y * leng[1] < min[1] || max[1] < zero_atom[i].y * leng[1]) && by == p) continue;
    if((zero_atom[i].z * leng[2] < min[2] || max[2] < zero_atom[i].z * leng[2]) && bz == p) continue;
    centercount++;
    aroundcount = 0;
    for(int j=0;j<all_atom;j++){
      if(zero_atom[j].name != around && around != "all") continue;
      if(i == j) continue;
      dx = fabs(zero_atom[i].x - zero_atom[j].x) * leng[0];
      dy = fabs(zero_atom[i].y - zero_atom[j].y) * leng[1];
      dz = fabs(zero_atom[i].z - zero_atom[j].z) * leng[2]; 
      r  = sunsq(dx, dy, dz);
      if(r < cutDEF){
        Rate = floor(r/dr);
        ++coord[Rate];
        ++aroundcount;
        continue;
        cerr << Rate << " def" << endl;
      }
      if(bx == p){
        rx = sunsq(leng[0]-dx, dy, dz);
        if(rx < cutDEF){
          Rate = floor(rx/dr);
          ++coord[Rate];
          ++aroundcount;
          cerr << Rate << " x" << endl;
          continue;
	}
        if(by == p){
          rxy = sunsq(leng[0]-dx, leng[1]-dy, dz);
          if(rxy < cutDEF){
            Rate = floor(rxy/dr);
            ++coord[Rate];
            ++aroundcount;
            cerr << Rate << " xy" << endl;
            continue;
          }
          if(bz == p){
            rxyz = sunsq(leng[0]-dx, leng[1]-dy, leng[2]-dz);
            if(rxyz < cutDEF){
              Rate = floor(rxyz/dr);
              ++coord[Rate];
              ++aroundcount;
              cerr << Rate << " xyz" << endl;
              continue;
            }
          }
        }
        if(bz == p){
          rxz = sunsq(leng[0]-dx, dy, leng[2]-dz);
          if(rxz < cutDEF){
            Rate = floor(rxz/dr);
            ++coord[Rate];
            ++aroundcount;
            cerr << Rate << " xz" << endl;
            continue;
          }
        }
      }
      if(by == p){
        ry = sunsq(dx, leng[1]-dy, dz);
        if(ry < cutDEF){
          Rate = floor(ry/dr);
          ++coord[Rate];
          ++aroundcount;
          cerr << Rate << " y" << endl;
          continue;
        }
        if(bz == p){
          ryz = sunsq(dx, leng[1]-dy, leng[2]-dz);
          if(ryz < cutDEF){
            Rate = floor(ryz/dr);
            ++coord[Rate];
            ++aroundcount;
            cerr << Rate << " yz" << endl;
            continue;
          }
        }
      }
      if(bz == p){
        rz = sunsq(dx, dy, leng[2]-dz);
        if(rz < cutDEF){
          Rate = floor(rz/dr);
          ++coord[Rate];
          ++aroundcount;
          cerr << Rate << " z" << endl;
          continue;
        }
      }
    }
    ++neighbor[aroundcount];
  }

  
  //output
  if(centercount == 0){
    cerr << " Err: This structure DON'T have center atoms!" << endl;
    return 1;
  }
  if(opts.atomratio == "rdf"){

    double rho = all_atom / ((leng[0]) * (leng[1]) * (leng[2])); //density of atoms [units / A^3] !cuboid only!
    double gDEF = 4 * acos(-1.0) *  dr * rho;
    double g = 0.000;

    cout << "# rdf " << center << "-" << around << " cutoff=" << cutDEF << " bin=" << dr << endl;
    if(bx == "p") cout << "# x periodric" << endl;
    else          cout << "# x " << min[0] << " " << max[0] << endl;
    if(by == "p") cout << "# y periodric" << endl;
    else          cout << "# y " << min[1] << " " << max[1] << endl;
    if(bz == "p") cout << "# z periodric" << endl;
    else          cout << "# z " << min[2] << " " << max[2] << endl;

    if(coord[0] != 0){
        cerr << " Error!" << endl;
        return 1;
    }

    for(int i = 1; i < ratiomax; i++){
      g=((double)coord[i]/(double)centercount)/(gDEF * pow(i*dr, 2.0));
      cout << i << " " << setw(4) << i*dr << " " << setw(4) << (i+1)*dr << " " << setw(4) << coord[i] << " " << g <<endl;
    }

  }else if(opts.atomratio == "num"){
    cout << "# Coordination " << center << "-" << around << " cutoff=" << cutDEF << " bin=" << dr << endl;
    if(bx == "p") cout << "# x periodric" << endl;
    else          cout << "# x " << min[0] << " " << max[0] << endl;
    if(by == "p") cout << "# y periodric" << endl;
    else          cout << "# y " << min[1] << " " << max[1] << endl;
    if(bz == "p") cout << "# z periodric" << endl;
    else          cout << "# z " << min[2] << " " << max[2] << endl;

    for(int i = 0; i < neighbormax; i++){
      cout << i << " " << neighbor[i] <<endl;
    }
  }
  cerr << "centercount:" << centercount << endl;
  delete[] zero_atom;
}

