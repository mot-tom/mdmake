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

#define MAXELEMENT 22
typedef enum {Si, Ge, C, Sn, Mo, S, O, SX, Ga, As, In, Au, Ni, Ti, Al, Mg, Sr, Y, N, Zr, B, Be, Na, Cl} AtomIonType;

struct info        // information of atom
{
  string name;     // atom name
  double x, y, z;  // coordinate of atom
  int lay;         // layer number of atom
  string option;   // options of atom
};

extern int shape(string shape, double x, double y, double z);

int stack(int argc, char **argv, optpara opts, string atom_lmp){
  const int max_atom = 300000;         // maximum number of atom

  int nx = 1, ny = 1, nz = 1;

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

  double leng1[3] = {0.0, 0.0, 0.0};    // length of each side of crystal  
  double tilt1[3][3] = {{0.0, 0.0, 0.0}, 
                       {0.0, 0.0, 0.0},  // x
                       {0.0, 0.0, 0.0}}; // tilt length of each side of crystal

  double leng2[3] = {0.0, 0.0, 0.0};    // length of each side of crystal  
  double tilt2[3][3] = {{0.0, 0.0, 0.0}, 
                       {0.0, 0.0, 0.0},  // x
                       {0.0, 0.0, 0.0}}; // tilt length of each side of crystal

  int atomnum[MAXELEMENT];             // number of atom, respectively
  for(int i=0;i<MAXELEMENT;i++){
    atomnum[i] = 0;
  }

  ifstream file_fin1(opts.stack.c_str());
  ifstream file_fin2(opts.atomcomp.c_str());
  
  if(!file_fin1)
  {
    cerr << "Can't open " << opts.stack << "!" << endl;
    return 1;
  }
  if(!file_fin2)
  {
    cerr << "Can't open " << opts.atomcomp << "!" << endl;
    return 1;
  }

  string temp_line, file_type;
  int line;
  
  cerr << "Read \"" << opts.stack << "\"." << endl;
  line = 0;
  while(getline(file_fin1, temp_line)){
    if(temp_line.find("Gear(") != string::npos && line == 0){file_type = "mdli";}
    line++;
  }
  file_fin1.close();
  if(file_type == "mdli") cerr << "\"" << opts.stack << "\" is input file for Mdlabo." << endl;
  else{
    cerr << "\"" << opts.stack << "\" is wrong file type!" << endl;
    return 1;
  }

  cerr << "Read \"" << opts.atomcomp << "\"." << endl;
  line = 0; file_type = "";
  while(getline(file_fin2, temp_line)){
    if(temp_line.find("Gear(") != string::npos && line == 0){file_type = "mdli";}
    line++;
  }
  file_fin2.close();
  if(file_type == "mdli") cerr << "\"" << opts.atomcomp << "\" is input file for Mdlabo." << endl;
  else{
    cerr << "\"" << opts.atomcomp << "\" is wrong file type!" << endl;
    return 1;
  }

  int blank = 0, lc = 0;
  string title1, title2, first;
  info* zero_atom = new info[max_atom];
  int all_atom = 0;
  double ax, ay, az;

  file_fin1.open(opts.stack.c_str());
    while(getline(file_fin1, temp_line)){
      istringstream line_iss(temp_line.c_str());

      if(blank == 2)
      {
        line_iss >> ax >> ay >> az;

        if(ax != 0.0 && lc == 0)
          leng1[0] = ax;
        if(ay != 0.0 && lc == 0)
          tilt1[0][1] = ay;
        if(az != 0.0 && lc == 0)
          tilt1[0][2] = az;

        if(ax != 0.0 && lc == 1)
          tilt1[1][0] = ax;
        if(ay != 0.0 && lc == 1)
          leng1[1] = ay;
        if(az != 0.0 && lc == 1)
          tilt1[1][2] = az;

        if(ax != 0.0 && lc == 2)
          tilt1[2][0] = ax;
        if(ay != 0.0 && lc == 2)
          tilt1[2][1] = ay;
        if(az != 0.0 && lc == 2)      
          leng1[2] = az;

        if(lc < 3){
          cerr << "|ul" << lc << "0 ul" << lc << "1 ul" << lc << "2|";
          if(lc == 1) cerr << "=";
          else cerr << " ";
          if(lc == 0)
            cerr << showpoint << "|" << setw(6) << leng1[0]    << " " << tilt1[0][1] << " " << tilt1[0][2] << "|";
          else if(lc == 1)
            cerr << showpoint << "|" << setw(6) << tilt1[1][0] << " " << leng1[1]    << " " << tilt1[1][2] << "|";
          else if(lc == 2)
            cerr << showpoint << "|" << setw(6) << tilt1[2][0] << " " << tilt1[2][1] << " " << leng1[2] << "|";
          cerr << endl;
         }
        lc++;
      }

      if(blank == 3 && temp_line != "")
     {
        if(all_atom+1 > max_atom)
        {
          cout << "Dynamic memory space is not enough for atom!" << endl;
          return 1;
        }
      
        line_iss >> first;
      
        if(all_atom == 0)
        {
          if(first != "Si" && first != "SX" && first != "Ge" && first != "Sn" && first != "C"
          && first != "Mo" && first != "S"  && first != "S1" && first != "S2" && first != "O"
          && first != "Ga" && first != "As" && first != "In" && first != "Au" && first != "Ni"
          && first != "Ti" && first != "Al" && first != "Mg" && first != "Sr" && first != "Y"
          && first != "N"  && first != "Zr" && first != "B"  && first != "Be" && first != "Na" && first != "Cl")
          {
            cout << "Setting of atom name is wrong! (" << first << ")" << endl;
            return 1;
          }
        }
      
        if(first != "V")
        {
          zero_atom[all_atom].name = first;
          line_iss >> ax >> ay >> az;
          zero_atom[all_atom].x =ax*leng1[0], zero_atom[all_atom].y = ay*leng1[1], zero_atom[all_atom].z = az*leng1[2];
          getline(line_iss, zero_atom[all_atom].option);
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
          if(first == "Na") //takizawa
            atomnum[Na]++;
          if(first == "Cl") //takizawa
            atomnum[Cl]++;
        }
      }

      if(temp_line == "")
        blank++;

      if(blank == 1) title1 = temp_line;
    }
  file_fin1.close();

  lc = 0;
  blank = 0;
  file_fin2.open(opts.atomcomp.c_str());
    while(getline(file_fin2, temp_line)){
      istringstream line_iss(temp_line.c_str());

      if(blank == 2)
      {
        line_iss >> ax >> ay >> az;

        if(ax != 0.0 && lc == 0)
          leng2[0] = ax;
        if(ay != 0.0 && lc == 0)
          tilt2[0][1] = ay;
        if(az != 0.0 && lc == 0)
          tilt2[0][2] = az;

        if(ax != 0.0 && lc == 1)
          tilt2[1][0] = ax;
        if(ay != 0.0 && lc == 1)
          leng2[1] = ay;
        if(az != 0.0 && lc == 1)
          tilt2[1][2] = az;

        if(ax != 0.0 && lc == 2)
          tilt2[2][0] = ax;
        if(ay != 0.0 && lc == 2)
          tilt2[2][1] = ay;
        if(az != 0.0 && lc == 2)      
          leng2[2] = az;

        if(lc < 3){
          cerr << "|ul" << lc << "0 ul" << lc << "1 ul" << lc << "2|";
          if(lc == 1) cerr << "=";
          else cerr << " ";
          if(lc == 0)
            cerr << showpoint << "|" << setw(6) << leng2[0]    << " " << tilt2[0][1] << " " << tilt2[0][2] << "|";
          else if(lc == 1)
            cerr << showpoint << "|" << setw(6) << tilt2[1][0] << " " << leng2[1]    << " " << tilt2[1][2] << "|";
          else if(lc == 2)
            cerr << showpoint << "|" << setw(6) << tilt2[2][0] << " " << tilt2[2][1] << " " << leng2[2] << "|";
          cerr << endl;
         }
        lc++;
      }

      if(blank == 3 && temp_line != "")
     {
        if(all_atom+1 > max_atom)
        {
          cout << "Dynamic memory space is not enough for atom!" << endl;
          return 1;
        }
      
        line_iss >> first;
      
        if(all_atom == 0)
        {
          if(first != "Si" && first != "SX" && first != "Ge" && first != "Sn" && first != "C"
          && first != "Mo" && first != "S"  && first != "S1" && first != "S2" && first != "O"
          && first != "Ga" && first != "As" && first != "In" && first != "Au" && first != "Ni"
          && first != "Ti" && first != "Al" && first != "Mg" && first != "Sr" && first != "Y"
          && first != "N"  && first != "Zr" && first != "B"  && first != "Be" && first != "Na" && first != "Cl")
          {
            cout << "Setting of atom name is wrong! (" << first << ")" << endl;
            return 1;
          }
        }
      
        if(first != "V")
        {
          zero_atom[all_atom].name = first;
          line_iss >> ax >> ay >> az;
          zero_atom[all_atom].x =ax*leng2[0]+leng1[0], zero_atom[all_atom].y = ay*leng1[1], zero_atom[all_atom].z = az*leng1[2];
          getline(line_iss, zero_atom[all_atom].option);
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
          if(first == "Na") //takizawa
            atomnum[Na]++;
          if(first == "Cl") //takizawa
            atomnum[Cl]++;
        }
      }

      if(temp_line == "")
        blank++;

      if(blank == 1) title2 = temp_line;
    }
  file_fin2.close();

  leng1[0] *= nx;
  leng2[0] *= nx;
  leng1[1] *= ny;
  leng1[2] *= nz;
  tilt1[0][1] *= nx;
  tilt1[0][2] *= nx;
  tilt1[1][0] *= ny;
  tilt1[1][2] *= ny;
  tilt1[2][0] *= nz;
  tilt1[2][1] *= nz;

  string atom_name;
  string atom_type[4];
  double atom_mass[4];
  int atom_num[4];
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
  }
  else if(atomnum[In] && atomnum[As]){
    all_type = 2;
    atom_name = "InAs";
    atom_type[0] = "In"; atom_mass[0] = 114.818; atom_num[0] = atomnum[Ga];
    atom_type[1] = "As"; atom_mass[1] = 74.9216; atom_num[1] = atomnum[As];
  }
  else if(atomnum[Si] && atomnum[Ge] && atomnum[Sn]){
    all_type = 3;
    atom_name = "GeSiSn";
    atom_type[0] = "Ge"; atom_mass[0] = 72.63;   atom_num[0] = atomnum[Ge];
    atom_type[1] = "Si"; atom_mass[1] = 28.0855; atom_num[1] = atomnum[Si];
    atom_type[2] = "Sn"; atom_mass[2] = 118.710; atom_num[2] = atomnum[Sn];
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
  }
  else if(atomnum[Si] && atomnum[O]){
    all_type = 2;
    atom_name = "SiO2";
    atom_type[0] = "Si"; atom_mass[0] = 28.0855; atom_num[0] = atomnum[Si];
    atom_type[1] = "O";  atom_mass[1] = 15.9994; atom_num[1] = atomnum[O];
  }
  else if(atomnum[Ge] && atomnum[O]){
    all_type = 2;
    atom_name = "GeO2";
    atom_type[0] = "Ge"; atom_mass[0] = 72.63;   atom_num[0] = atomnum[Ge];
    atom_type[1] = "O";  atom_mass[1] = 15.9994; atom_num[1] = atomnum[O];
  }
  else if(atomnum[Sn] && atomnum[O]){
    all_type = 2;
    atom_name = "SnO2";
    atom_type[0] = "Sn"; atom_mass[0] = 118.710; atom_num[0] = atomnum[Sn];
    atom_type[1] = "O";  atom_mass[1] = 15.9994; atom_num[1] = atomnum[O];
  }
  else if(atomnum[Na] && atomnum[Cl]){
    all_type = 2;
    atom_name = "NaCl";
    atom_type[0] = "Na";  atom_mass[0] = 22.9897; atom_num[0] = atomnum[Na];
    atom_type[1] = "Cl";  atom_mass[1] = 35.446; atom_num[1] = atomnum[Cl];
  }
  else if(atomnum[Ti] && atomnum[O]){
    all_type = 2;
    atom_name = "TiO2";
    atom_type[0] = "Ti"; atom_mass[0] = 47.867;  atom_num[0] = atomnum[Ti];
    atom_type[1] = "O";  atom_mass[1] = 15.9994; atom_num[1] = atomnum[O];
  }
  else if(atomnum[Zr] && atomnum[O]){
    all_type = 2;
    atom_name = "ZrO2";
    atom_type[0] = "Zr"; atom_mass[0] = 91.224;  atom_num[0] = atomnum[Zr];
    atom_type[1] = "O";  atom_mass[1] = 15.9994; atom_num[1] = atomnum[O];
  }
  else if(atomnum[Mg] && atomnum[O]){
    all_type = 2;
    atom_name = "MgO";
    atom_type[0] = "Mg"; atom_mass[0] = 24.3050; atom_num[0] = atomnum[Mg];
    atom_type[1] = "O";  atom_mass[1] = 15.9994; atom_num[1] = atomnum[O];
  }
  else if(atomnum[Sr] && atomnum[O]){
    all_type = 2;
    atom_name = "SrO";
    atom_type[0] = "Sr"; atom_mass[0] = 87.62;   atom_num[0] = atomnum[Mg];
    atom_type[1] = "O";  atom_mass[1] = 15.9994; atom_num[1] = atomnum[O];
  }
  else if(atomnum[Be] && atomnum[O]){
    all_type = 2;
    atom_name = "BeO";
    atom_type[0] = "Be"; atom_mass[0] = 9.01218; atom_num[0] = atomnum[Be];
    atom_type[1] = "O";  atom_mass[1] = 15.9994; atom_num[1] = atomnum[O];
  }
  else if(atomnum[B] && atomnum[N]){
    all_type = 2;
    atom_name = "BN";
    atom_type[0] = "B";  atom_mass[0] = 10.811;  atom_num[0] = atomnum[B];
    atom_type[1] = "N";  atom_mass[1] = 14.0067; atom_num[1] = atomnum[N];
  }
  else if(atomnum[Al] && atomnum[N]){
    all_type = 2;
    atom_name = "AlN";
    atom_type[0] = "Al"; atom_mass[0] = 26.9815; atom_num[0] = atomnum[Al];
    atom_type[1] = "N";  atom_mass[1] = 14.0067; atom_num[1] = atomnum[N];
  }
  else if(atomnum[Al] && atomnum[O]){
    all_type = 2;
    atom_name = "Al2O3";
    atom_type[0] = "Al"; atom_mass[0] = 26.9815; atom_num[0] = atomnum[Al];
    atom_type[1] = "O";  atom_mass[1] = 15.9994; atom_num[1] = atomnum[O];
  }
  else if(atomnum[Y] && atomnum[O]){
    all_type = 2;
    atom_name = "Y2O3";
    atom_type[0] = "Y";  atom_mass[0] = 88.9059; atom_num[0] = atomnum[Y];
    atom_type[1] = "O";  atom_mass[1] = 15.9994; atom_num[1] = atomnum[O];
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

    double lengall0 = leng1[0]+leng2[0];
    cout << "Gear(5) RAND=10 " << opts.potential << " T=300.0 STEP=(0,1000100)" << endl;
    cout << "BLOCK FIXANGLE W=INF Q=INF LOGSTEP=10 PRINTVELOCITY dt=0.1fs scratch(100)" << endl << endl;
    cout << title1 << " + " << title2 << endl << endl;
    cout << lengall0 << " " << tilt1[0][1] << " " << tilt1[0][2] << opts.pbc_x << endl;
    cout << tilt1[1][0] << " " << leng1[1] << " " << tilt1[1][2] << opts.pbc_y << endl;
    cout << tilt1[2][0] << " " << tilt1[2][1] << " " << leng1[2] << opts.pbc_z << endl << endl;
    for(int j=0; j<nx; j++) for(int k=0; k<ny; k++) for(int l=0; l<nz; l++)
    for(int i=0;i<all_atom;i++){
      double ax = zero_atom[i].x/lengall0+(double)j/nx, ay = zero_atom[i].y/leng1[1]+(double)k/ny, az = zero_atom[i].z/leng1[2]+(double)l/nz;
      if(shape(opts.shape, ax, ay, az)){
        cout << zero_atom[i].name << " " << ax << " " << ay << " " << az;
        if(zero_atom[i].option != "") cout << zero_atom[i].option;
        cout << endl;
      }
    }

  delete[] zero_atom;
}

