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

#define MAXELEMENT 27
//takizawa
typedef enum {Si, Ge, C, Sn, Mo, S, O, SX, Ga, As, In, Au, Ni, Ti, Al, Mg, Sr, Y, N, Zr, B, Be, Ar, Sb, Te, Na, Cl} AtomIonType;

struct info        // information of atom
{
  string name;     // atom name
  double x, y, z;  // coordinate of atom
  int lay;         // layer number of atom
  string option;   // options of atom
};

extern int shape(string shape, double x, double y, double z);

int convert(int argc, char **argv, optpara opts, string atom_lmp){
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
    if(temp_line.find("#######################################################################") != string::npos && line == 0){file_type = "mdlo";}
    if(temp_line.find("ITEM: TIMESTEP") != string::npos && line == 0){file_type = "lmpo";}
    line++;
  }
  file_fin.close();

  if(file_type == "mdli") cerr << "\"" << opts.atomcomp << "\" is input file for Mdlabo." << endl;
  else if(file_type == "mdlo") cerr << "\"" << opts.atomcomp << "\" is output file of Mdlabo." << endl;
  else if(file_type == "lmpo"){ //takizawa
    cerr << "\"" << opts.atomcomp << "\" is output file of Lammps." << endl;
    if(atom_lmp != "Si" && atom_lmp != "Ge" && atom_lmp != "Sn" && atom_lmp != "C" && atom_lmp != "Ni" && atom_lmp != "Au"
    && atom_lmp != "SiGe" && atom_lmp != "GeSn"  && atom_lmp != "GeSiSn" && atom_lmp != "GeSbTe" && atom_lmp != "SiC"
    && atom_lmp != "SiO2" && atom_lmp != "GeO2" && atom_lmp != "SnO2" && atom_lmp != "SiO2Ar" && atom_lmp != "SiAr" 
    && atom_lmp != "Si3N4" && atom_lmp != "TiO2" && atom_lmp != "ZrO2" && atom_lmp != "MoS2"
    && atom_lmp != "TiN" && atom_lmp != "AlN" && atom_lmp != "BeO" && atom_lmp != "SrO" && atom_lmp != "MgO"
    && atom_lmp != "Y2O3" && atom_lmp != "Al2O3" && atom_lmp != "NaCl" &&  atom_lmp != "NaCl-SiO2"
    && atom_lmp != "Al2O3-SiO2" && atom_lmp != "TiO2-SiO2" && atom_lmp != "SrO-SiO2" && atom_lmp != "MgO-SiO2"){ //takizawa
      cerr << "\"" << opts.atomcomp << "\" is wrong atom type for lammps!" << endl;
      return 1;
    }
  }
  else{
    cerr << "\"" << opts.atomcomp << "\" is wrong file type!" << endl;
    return 1;
  }

  if(opts.atomratio == "mdl") cerr << "\"" << opts.atomcomp << "\" convert to mdl file." << endl;
  else if(opts.atomratio == "mdlZtoX") cerr << "\"" << opts.atomcomp << "\" convert to ZtoX rotated mdl file." << endl;
  else if(opts.atomratio == "mdlYtoX") cerr << "\"" << opts.atomcomp << "\" convert to YtoX rotated mdl file." << endl;
  else if(opts.atomratio == "mdlCUT") cerr << "\"" << opts.atomcomp << "\" convert to mdl file and delete atoms which are out of box." << endl;
  else if(opts.atomratio == "mdlFILM") cerr << "\"" << opts.atomcomp << "\" convert to mdl file resized 3 times to z direction." << endl;
  else if(opts.atomratio == "mdlWIRE") cerr << "\"" << opts.atomcomp << "\" convert to mdl file resized 3 times to all direction." << endl;
  else if(opts.atomratio == "mdlFIX") cerr << "\"" << opts.atomcomp << "\" convert to mdlFIX file." << endl;
  else if(opts.atomratio == "mdlXO2") cerr << "\"" << opts.atomcomp << "\" convert to mdlXO2 file." << endl;
  else if(opts.atomratio == "mdlCLEAN") cerr << "\"" << opts.atomcomp << "\" convert to mdl file and delete atoms away from film." << endl; //takizawa
  else if(opts.atomratio == "mdlZRES") cerr << "\"" << opts.atomcomp << "\" convert to resized mdl file." << endl; //takizawa
  else if(opts.atomratio == "mdlGROUP") cerr << "\"" << opts.atomcomp << "\" convert to group_mdl file." << endl; //takizawa
  else if(opts.atomratio == "mdlINTERFACE") cerr << "\"" << opts.atomcomp << "\" convert to mdl file and delete atoms near by interface." << endl; //takizawa
  else if(opts.atomratio == "lmp") cerr << "\"" << opts.atomcomp << "\" convert to lmp file." << endl;
  else if(opts.atomratio == "lmpCIM") cerr << "\"" << opts.atomcomp << "\" convert to lmp file with CIM potential." << endl;
  else if(opts.atomratio == "lmpCMAS") cerr << "\"" << opts.atomcomp << "\" convert to lmp file with CMAS potential." << endl;
  else if(opts.atomratio == "lmpAr") cerr << "\"" << opts.atomcomp << "\" convert to Ar ion added lmp file." << endl;
  else if(opts.atomratio == "lmpLAY") cerr << "\"" << opts.atomcomp << "\" convert to lmp file with oxidation layer property." << endl;
  else if(opts.atomratio == "xyz") cerr << "\"" << opts.atomcomp << "\" convert to xyz file." << endl;
  else if(opts.atomratio == "vasp") cerr << "\"" << opts.atomcomp << "\" convert to vasp file." << endl;
  else if(opts.atomratio == "alm") cerr << "\"" << opts.atomcomp << "\" convert to alm file." << endl;
  else{
    cerr << "Output type of \"" << opts.atomratio << "\" is wrong!" << endl;
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
          tilt[1][2] = zero_atom[0].z;

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
          && first != "N"  && first != "Zr" && first != "B"  && first != "Be" && first != "Ar"
	  && first != "Sb" && first != "Te" && first != "Na" && first != "Cl") //takizawa
          {
            cout << "Setting of atom name is wrong! (" << first << ")" << endl;
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
  }
  else if(file_type == "mdlo"){
    while(temp_line.find("****  Final Structure ****") == string::npos){
      getline(file_fin, temp_line);
    }
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
          tilt[1][2] = zero_atom[0].z;

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
          && first != "N"  && first != "Zr" && first != "B"  && first != "Be" && first != "Ar"
	  && first != "Sb" && first != "Te" && first != "Na" && first != "Cl") //takizawa
          {
            cout << "Setting of atom name is wrong! (" << first << ")" << endl;
            return 1;
          }
        }
      
        if(first != "V")
        {
          zero_atom[all_atom].name = first;
          line_iss >> zero_atom[all_atom].x >> zero_atom[all_atom].y >> zero_atom[all_atom].z;
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
  }
  else if(file_type == "lmpo"){
    title = "Imput file made from lammps result";
    getline(file_fin, temp_line);
    getline(file_fin, temp_line);
    getline(file_fin, temp_line);
    getline(file_fin, temp_line);
    if(1){
      istringstream line_iss(temp_line.c_str());
      line_iss >> all_atom;
    }

    if(all_atom+1 > max_atom)
    {
      cout << "Dynamic memory space is not enough for atom!" << endl;
      return 1;
    }

    getline(file_fin, temp_line);
    getline(file_fin, temp_line);
    if(1){
      istringstream line_iss(temp_line.c_str());
      line_iss >> zero_atom[0].x >> zero_atom[0].y >> zero_atom[0].z;
      leng[0] = zero_atom[0].y - zero_atom[0].x;
      tilt[0][1] = zero_atom[0].z * 0.5;
      tilt[1][0] = zero_atom[0].z * 0.5;
    }

    getline(file_fin, temp_line);
    if(1){
      istringstream line_iss(temp_line.c_str());
      line_iss >> zero_atom[0].x >> zero_atom[0].y >> zero_atom[0].z;
      leng[1] = zero_atom[0].y - zero_atom[0].x;
      tilt[0][2] = zero_atom[0].z * 0.5;
      tilt[2][0] = zero_atom[0].z * 0.5;
    }

    getline(file_fin, temp_line);
    if(1){
      istringstream line_iss(temp_line.c_str());
      line_iss >> zero_atom[0].x >> zero_atom[0].y >> zero_atom[0].z;
      leng[2] = zero_atom[0].y - zero_atom[0].x;
      tilt[1][2] = zero_atom[0].z * 0.5;
      tilt[2][1] = zero_atom[0].z * 0.5;
    }
    cerr << "|ul00 ul01 ul02| ";
    cerr << showpoint << "|" << setw(6) << leng[0]    << " " << tilt[0][1] << " " << tilt[0][2] << "|" << endl;
    cerr << "|ul10 ul11 ul12|=";
    cerr << showpoint << "|" << setw(6) << tilt[1][0] << " " << leng[1]    << " " << tilt[1][2] << "|" << endl;
    cerr << "|ul20 ul21 ul22| ";
    cerr << showpoint << "|" << setw(6) << tilt[2][0] << " " << tilt[2][1] << " " << leng[2] << "|" << endl;

    getline(file_fin, temp_line);
    for(int i=0; i<all_atom; i++){
      int label = 0, lay = 0, atom_type;
      getline(file_fin, temp_line);
      istringstream line_iss(temp_line.c_str());
      line_iss >> label >> atom_type;
      line_iss >> zero_atom[label-1].x >> zero_atom[label-1].y >> zero_atom[label-1].z >> lay;
      ostringstream os;

      if(atom_lmp == "MoS2") {
        if(atom_type == 1) {atomnum[S]++; zero_atom[label-1].name = "S";}
        if(atom_type == 2) {atomnum[S]++; zero_atom[label-1].name = "S";}
        if(atom_type == 3) {atomnum[Mo]++; zero_atom[label-1].name = "Mo";}
      }
      else if(atom_lmp == "GeSiSn") {
        if(atom_type == 1) {atomnum[Ge]++; zero_atom[label-1].name = "Ge";}
        if(atom_type == 2) {atomnum[Si]++; zero_atom[label-1].name = "Si";}
        if(atom_type == 3) {atomnum[Sn]++; zero_atom[label-1].name = "Sn";}
      }
      else if(atom_lmp == "GeSbTe") { //takizawa
        if(atom_type == 1) {atomnum[Ge]++; zero_atom[label-1].name = "Ge";}
        if(atom_type == 2) {atomnum[Sb]++; zero_atom[label-1].name = "Sb";}
        if(atom_type == 3) {atomnum[Te]++; zero_atom[label-1].name = "Te";}
      }
      else if(atom_lmp == "SiC") {
        if(atom_type == 1) {atomnum[Si]++; zero_atom[label-1].name = "Si";}
        if(atom_type == 2) {atomnum[C]++; zero_atom[label-1].name = "C";}
      }
      else if(atom_lmp == "SiGe") {
        if(atom_type == 1) {atomnum[Si]++; zero_atom[label-1].name = "Si";}
        if(atom_type == 2) {atomnum[Ge]++; zero_atom[label-1].name = "Ge";}
      }
      else if(atom_lmp == "GeSn") {
        if(atom_type == 1) {atomnum[Ge]++; zero_atom[label-1].name = "Ge";}
        if(atom_type == 2) {atomnum[Sn]++; zero_atom[label-1].name = "Sn";}
      }
      else if(atom_lmp == "Al2O3-SiO2") {
        if(atom_type == 1) {atomnum[Al]++; zero_atom[label-1].name = "Al";}
        if(atom_type == 2) {atomnum[Si]++; zero_atom[label-1].name = "Si";}
        if(atom_type == 3) {atomnum[O]++; zero_atom[label-1].name = "O";}
      }
      else if(atom_lmp == "Al2O3") {
        if(atom_type == 1) {atomnum[Al]++; zero_atom[label-1].name = "Al";}
        if(atom_type == 2) {atomnum[O]++; zero_atom[label-1].name = "O";}
      }
      else if(atom_lmp == "Y2O3") {
        if(atom_type == 1) {atomnum[Y]++; zero_atom[label-1].name = "Y";}
        if(atom_type == 2) {atomnum[O]++; zero_atom[label-1].name = "O";}
      }
      else if(atom_lmp == "MgO-SiO2") {
        if(atom_type == 1) {atomnum[Mg]++; zero_atom[label-1].name = "Mg";}
        if(atom_type == 2) {atomnum[Si]++; zero_atom[label-1].name = "Si";}
        if(atom_type == 3) {atomnum[O]++; zero_atom[label-1].name = "O";}
      }
      else if(atom_lmp == "MgO") {
        if(atom_type == 1) {atomnum[Mg]++; zero_atom[label-1].name = "Mg";}
        if(atom_type == 2) {atomnum[O]++; zero_atom[label-1].name = "O";}
      }
      else if(atom_lmp == "SrO-SiO2") {
        if(atom_type == 1) {atomnum[Sr]++; zero_atom[label-1].name = "Sr";}
        if(atom_type == 2) {atomnum[Si]++; zero_atom[label-1].name = "Si";}
        if(atom_type == 3) {atomnum[O]++; zero_atom[label-1].name = "O";}
      }
      else if(atom_lmp == "SrO") {
        if(atom_type == 1) {atomnum[Sr]++; zero_atom[label-1].name = "Sr";}
        if(atom_type == 2) {atomnum[O]++; zero_atom[label-1].name = "O";}
      }
      else if(atom_lmp == "BeO") {
        if(atom_type == 1) {atomnum[Be]++; zero_atom[label-1].name = "Be";}
        if(atom_type == 2) {atomnum[O]++; zero_atom[label-1].name = "O";}
      }
      else if(atom_lmp == "AlN") {
        if(atom_type == 1) {atomnum[Al]++; zero_atom[label-1].name = "Al";}
        if(atom_type == 2) {atomnum[N]++; zero_atom[label-1].name = "N";}
      }
      else if(atom_lmp == "TiN") {
        if(atom_type == 1) {atomnum[Ti]++; zero_atom[label-1].name = "Ti";}
        if(atom_type == 2) {atomnum[N]++; zero_atom[label-1].name = "N";}
      }
      else if(atom_lmp == "TiO2-SiO2") {
        if(atom_type == 1) {atomnum[Ti]++; zero_atom[label-1].name = "Ti";}
        if(atom_type == 2) {atomnum[Si]++; zero_atom[label-1].name = "Si";}
        if(atom_type == 3) {atomnum[O]++; zero_atom[label-1].name = "O";}
      }
      else if(atom_lmp == "TiO2") {
        if(atom_type == 1) {atomnum[Ti]++; zero_atom[label-1].name = "Ti";}
        if(atom_type == 2) {atomnum[O]++; zero_atom[label-1].name = "O";}
      }
      else if(atom_lmp == "ZrO2") {
        if(atom_type == 1) {atomnum[Zr]++; zero_atom[label-1].name = "Zr";}
        if(atom_type == 2) {atomnum[O]++; zero_atom[label-1].name = "O";}
      }
      else if(atom_lmp == "Si3N4") {
        if(atom_type == 1) {atomnum[Si]++; zero_atom[label-1].name = "Si";}
        if(atom_type == 2) {atomnum[N]++; zero_atom[label-1].name = "N";}
        if(lay == 255) {zero_atom[label-1].lay = lay; zero_atom[label-1].option = " F";}
        else if(lay != 0) {zero_atom[label-1].lay = lay; os << " G " << lay; zero_atom[label-1].option = os.str();}
      }
      else if(atom_lmp == "SiO2Ar") {
        if(atom_type == 1) {atomnum[Si]++; zero_atom[label-1].name = "Si";}
        if(atom_type == 2) {atomnum[O]++; zero_atom[label-1].name = "O";}
        if(atom_type == 3) {atomnum[Ar]++; zero_atom[label-1].name = "Ar";}
        if(lay == 255) {zero_atom[label-1].lay = lay; zero_atom[label-1].option = " F";}
        else if(lay != 0) {zero_atom[label-1].lay = lay; os << " G " << lay; zero_atom[label-1].option = os.str();}
      }
      else if(atom_lmp == "SiAr") {
        if(atom_type == 1) {atomnum[Si]++; zero_atom[label-1].name = "Si";}
        if(atom_type == 2) {atomnum[Ar]++; zero_atom[label-1].name = "Ar";}
        if(lay == 255) {zero_atom[label-1].lay = lay; zero_atom[label-1].option = " F";}
        else if(lay != 0) {zero_atom[label-1].lay = lay; os << " G " << lay; zero_atom[label-1].option = os.str();}
      }
      else if(atom_lmp == "SiO2") {
        if(atom_type == 1) {atomnum[Si]++; zero_atom[label-1].name = "Si";}
        if(atom_type == 2) {atomnum[O]++; zero_atom[label-1].name = "O";}
        if(lay == 255) {zero_atom[label-1].lay = lay; zero_atom[label-1].option = " F";}
        else if(lay != 0) {zero_atom[label-1].lay = lay; os << " G " << lay; zero_atom[label-1].option = os.str();}
      }
      else if(atom_lmp == "GeO2") {
        if(atom_type == 1) {atomnum[Ge]++; zero_atom[label-1].name = "Ge";}
        if(atom_type == 2) {atomnum[O]++; zero_atom[label-1].name = "O";}
        if(lay == 255) {zero_atom[label-1].lay = lay; zero_atom[label-1].option = " F";}
        else if(lay != 0) {zero_atom[label-1].lay = lay; os << " G " << lay; zero_atom[label-1].option = os.str();}
      }
      else if(atom_lmp == "SnO2") {
        if(atom_type == 1) {atomnum[Sn]++; zero_atom[label-1].name = "Sn";}
        if(atom_type == 2) {atomnum[O]++; zero_atom[label-1].name = "O";}
        if(lay == 255) {zero_atom[label-1].lay = lay; zero_atom[label-1].option = " F";}
        else if(lay != 0) {zero_atom[label-1].lay = lay; os << " G " << lay; zero_atom[label-1].option = os.str();}
      }
      else if(atom_lmp == "NaCl-SiO2") {
        if(atom_type == 1) {atomnum[Na]++; zero_atom[label-1].name = "Na";}
        if(atom_type == 2) {atomnum[Cl]++; zero_atom[label-1].name = "Cl";}
        if(atom_type == 3) {atomnum[Si]++; zero_atom[label-1].name = "Si";}
        if(atom_type == 4) {atomnum[O]++; zero_atom[label-1].name = "O";}
        if(lay == 255) {zero_atom[label-1].lay = lay; zero_atom[label-1].option = " F";}
        else if(lay != 0) {zero_atom[label-1].lay = lay; os << " G " << lay; zero_atom[label-1].option = os.str();}
      }
      else if(atom_lmp == "NaCl") {
        if(atom_type == 1) {atomnum[Na]++; zero_atom[label-1].name = "Na";}
        if(atom_type == 2) {atomnum[Cl]++; zero_atom[label-1].name = "Cl";}
        if(lay == 255) {zero_atom[label-1].lay = lay; zero_atom[label-1].option = " F";}
        else if(lay != 0) {zero_atom[label-1].lay = lay; os << " G " << lay; zero_atom[label-1].option = os.str();}
      }
      else if(atom_lmp == "Si") {
        if(atom_type == 1) {atomnum[Si]++; zero_atom[label-1].name = "Si";}
        if(lay == 255) {zero_atom[label-1].lay = lay; zero_atom[label-1].option = " F";}
        else if(lay != 0) {zero_atom[label-1].lay = lay; os << " G " << lay; zero_atom[label-1].option = os.str();}
      }
      else if(atom_lmp == "Ge") {
        if(atom_type == 1) {atomnum[Ge]++; zero_atom[label-1].name = "Ge";}
        if(lay == 255) {zero_atom[label-1].lay = lay; zero_atom[label-1].option = " F";}
        else if(lay != 0) {zero_atom[label-1].lay = lay; os << " G " << lay; zero_atom[label-1].option = os.str();}
      }
      else if(atom_lmp == "Sn") {
        if(atom_type == 1) {atomnum[Sn]++; zero_atom[label-1].name = "Sn";}
      }
      else if(atom_lmp == "C") {
        if(atom_type == 1) {atomnum[C]++; zero_atom[label-1].name = "C";}
      }
      else if(atom_lmp == "Au") {
        if(atom_type == 1) {atomnum[Au]++; zero_atom[label-1].name = "Au";}
      }
      else if(atom_lmp == "Ni") {
        if(atom_type == 1) {atomnum[Ni]++; zero_atom[label-1].name = "Ni";}
      }
    }
  }
  file_fin.close();
  cerr << "check1" << endl;

  if(opts.atomratio != "mdlGROUP"){ //takizawa
  	leng[0] *= nx;
  	leng[1] *= ny;
  	leng[2] *= nz;
  	tilt[0][1] *= nx;
  	tilt[0][2] *= nx;
  	tilt[1][0] *= ny;
  	tilt[1][2] *= ny;
  	tilt[2][0] *= nz;
  	tilt[2][1] *= nz;
  }

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
    atom_type[0] = "In"; atom_mass[0] = 114.818; atom_num[0] = atomnum[In];
    atom_type[1] = "As"; atom_mass[1] = 74.9216; atom_num[1] = atomnum[As];

    atom_CIM[0] = 3.0; atom_CMAS[0] = atom_CIM[0]*0.4725; atom_PIM[0] = atom_CIM[0]*0.6;
    atom_CIM[1] =-3.0; atom_CMAS[1] = atom_CIM[1]*0.4725; atom_PIM[1] = atom_CIM[0]*0.6;
  }
  else if(atomnum[Si] && atomnum[Ge] && atomnum[Sn]){ //takizawa
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
    atom_type[0] = "Ge"; atom_mass[0] = 72.63;  atom_num[0] = atomnum[Ge];
    atom_type[1] = "Sb"; atom_mass[1] = 121.76; atom_num[1] = atomnum[Sb];
    atom_type[2] = "Te"; atom_mass[2] = 127.60; atom_num[2] = atomnum[Te];
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
  else if(atomnum[Ti] && atomnum[Si] && atomnum[O]){ //takizawa
    all_type = 3;
    atom_name = "TiO2-SiO2";
    atom_type[0] = "Ti"; atom_mass[0] = 47.867;  atom_num[0] = atomnum[Ti];
    atom_type[1] = "Si"; atom_mass[1] = 28.0855; atom_num[1] = atomnum[Si];
    atom_type[2] = "O";  atom_mass[2] = 15.9994; atom_num[2] = atomnum[O];

    atom_CIM[0] = 4.0; atom_CMAS[0] = atom_CIM[0]*0.4725; atom_PIM[0] = atom_CIM[0]*0.6;
    atom_CIM[1] = 4.0; atom_CMAS[1] = atom_CIM[1]*0.4725; atom_PIM[1] = atom_CIM[1]*0.6;
    atom_CIM[2] =-2.0; atom_CMAS[2] = atom_CIM[2]*0.4725; atom_PIM[2] = atom_CIM[2]*0.6;
  }
  else if(atomnum[Mg] && atomnum[Si] && atomnum[O]){ //takizawa
    all_type = 3;
    atom_name = "MgO-SiO2";
    atom_type[0] = "Mg"; atom_mass[0] = 24.3050; atom_num[0] = atomnum[Mg];
    atom_type[1] = "Si"; atom_mass[1] = 28.0855; atom_num[1] = atomnum[Si];
    atom_type[2] = "O";  atom_mass[2] = 15.9994; atom_num[2] = atomnum[O]; 

    atom_CIM[0] = 2.0; atom_CMAS[0] = atom_CIM[0]*0.4725; atom_PIM[0] = atom_CIM[0]*0.6;
    atom_CIM[1] = 4.0; atom_CMAS[1] = atom_CIM[1]*0.4725; atom_PIM[1] = atom_CIM[1]*0.6;
    atom_CIM[2] =-2.0; atom_CMAS[2] = atom_CIM[2]*0.4725; atom_PIM[2] = atom_CIM[1]*0.6;
  }
  else if(atomnum[Sr] && atomnum[Si] && atomnum[O]){ //takizawa
    all_type = 3;
    atom_name = "SrO-SiO2";
    atom_type[0] = "Sr"; atom_mass[0] = 87.62;   atom_num[0] = atomnum[Sr];
    atom_type[1] = "Si"; atom_mass[1] = 28.0855; atom_num[1] = atomnum[Si];
    atom_type[2] = "O";  atom_mass[2] = 15.9994; atom_num[2] = atomnum[O];

    atom_CIM[0] = 2.0; atom_CMAS[0] = atom_CIM[0]*0.4725; atom_PIM[0] = atom_CIM[0]*0.6;
    atom_CIM[1] = 4.0; atom_CMAS[1] = atom_CIM[1]*0.4725; atom_PIM[1] = atom_CIM[1]*0.6;
    atom_CIM[2] =-2.0; atom_CMAS[2] = atom_CIM[2]*0.4725; atom_PIM[2] = atom_CIM[2]*0.6;
  }
  else if(atomnum[Si] && atomnum[SX] && atomnum[O]){ //takizawa
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
  else if(atomnum[Na] && atomnum[Cl]){
    all_type = 2;
    atom_name = "NaCl";
    atom_type[0] = "Na";  atom_mass[0] = 22.9897; atom_num[0] = atomnum[Na];
    atom_type[1] = "Cl";  atom_mass[1] = 35.446; atom_num[1] = atomnum[Cl];

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

  cerr << "crystal type " << atom_name << endl;

  string XO2_atom = "SX";
  double XO2_ratio = 0.97;
  double FIX_ratio = 0.00;
  double FILM_ratio = 3; //takizawa
  if(opts.atomratio == "mdl"){
    cout << "Gear(5) RAND=10 " << opts.potential << " T=300.0 STEP=(0,1000100)" << endl;
    cout << "BLOCK FIXANGLE W=INF Q=INF LOGSTEP=10 PRINTVELOCITY dt=0.1fs scratch(100)" << endl << endl;
    cout << title << endl << endl;
    cout << leng[0] << " " << tilt[0][1] << " " << tilt[0][2] << opts.pbc_x << endl;
    cout << tilt[1][0] << " " << leng[1] << " " << tilt[1][2] << opts.pbc_y << endl;
    cout << tilt[2][0] << " " << tilt[2][1] << " " << leng[2] << opts.pbc_z << endl << endl;
    for(int j=0; j<nx; j++) for(int k=0; k<ny; k++) for(int l=0; l<nz; l++)
    for(int i=0;i<all_atom;i++){
      double ax = (zero_atom[i].x+j)/nx, ay = (zero_atom[i].y+k)/ny, az = (zero_atom[i].z+l)/nz;
      if(shape(opts.shape, ax, ay, az)){
        cout << zero_atom[i].name << " " << ax << " " << ay << " " << az;
        if(zero_atom[i].option != "") cout << zero_atom[i].option;
        cout << endl;
      }
    }
  }
  else if(opts.atomratio == "mdlYtoX"){
    cout << "Gear(5) RAND=10 " << opts.potential << " T=300.0 STEP=(0,1000100)" << endl;
    cout << "BLOCK FIXANGLE W=INF Q=INF LOGSTEP=10 PRINTVELOCITY dt=0.1fs scratch(100)" << endl << endl;
    cout << title << endl << endl;
    cout << leng[1] << " " << tilt[1][0] << " " << tilt[1][2] << opts.pbc_x << endl;
    cout << tilt[0][1] << " " << leng[0] << " " << tilt[0][2] << opts.pbc_y << endl;
    cout << tilt[2][1] << " " << tilt[2][0] << " " << leng[2] << opts.pbc_z << endl << endl;
    for(int j=0; j<nx; j++) for(int k=0; k<ny; k++) for(int l=0; l<nz; l++)
    for(int i=0;i<all_atom;i++){
      double ax = (zero_atom[i].y+j)/nx, ay = (zero_atom[i].x+k)/ny, az = (zero_atom[i].z+l)/nz;
      if(shape(opts.shape, ax, ay, az)){
        cout << zero_atom[i].name << " " << ax << " " << ay << " " << az;
        if(zero_atom[i].option != "") cout << zero_atom[i].option;
        cout << endl;
      }
    }
  }
  else if(opts.atomratio == "mdlZtoX"){
    cout << "Gear(5) RAND=10 " << opts.potential << " T=300.0 STEP=(0,1000100)" << endl;
    cout << "BLOCK FIXANGLE W=INF Q=INF LOGSTEP=10 PRINTVELOCITY dt=0.1fs scratch(100)" << endl << endl;
    cout << title << endl << endl;
    cout << leng[2] << " " << tilt[2][1] << " " << tilt[2][0] << opts.pbc_x << endl;
    cout << tilt[1][2] << " " << leng[1] << " " << tilt[1][0] << opts.pbc_y << endl;
    cout << tilt[0][2] << " " << tilt[0][1] << " " << leng[0] << opts.pbc_z << endl << endl;
    for(int j=0; j<nx; j++) for(int k=0; k<ny; k++) for(int l=0; l<nz; l++)
    for(int i=0;i<all_atom;i++){
      double ax = (zero_atom[i].z+j)/nx, ay = (zero_atom[i].y+k)/ny, az = (zero_atom[i].x+l)/nz;
      if(shape(opts.shape, ax, ay, az)){
        cout << zero_atom[i].name << " " << ax << " " << ay << " " << az;
        if(zero_atom[i].option != "") cout << zero_atom[i].option;
        cout << endl;
      }
    }
  }
  else if(opts.atomratio == "mdlCUT"){
    cout << "OPT(P) RAND=10 " << opts.potential << " T=300.0 STEP=(0,20000)" << endl;
    cout << "BLOCK W=INF FIXANGLE Q=1.0" << endl << endl;
    cout << title << endl << endl;
    cout << leng[0] << " " << tilt[0][1] << " " << tilt[0][2] << opts.pbc_x << endl;
    cout << tilt[1][0] << " " << leng[1] << " " << tilt[1][2] << opts.pbc_y << endl;
    cout << tilt[2][0] << " " << tilt[2][1] << " " << leng[2] << opts.pbc_z << endl << endl;
    for(int j=0; j<nx; j++) for(int k=0; k<ny; k++) for(int l=0; l<nz; l++)
    for(int i=0;i<all_atom;i++){
      double ax = (zero_atom[i].x+j)/nx, ay = (zero_atom[i].y+k)/ny, az = (zero_atom[i].z+l)/nz;
      if(shape(opts.shape, ax, ay, az) && az >= 0 && az <= 1){
        cout << zero_atom[i].name << " " << ax << " " << ay << " " << az << endl;
        if(zero_atom[i].option != "") cout << zero_atom[i].option;
        cout << endl;
      }
    }
  }
  if(opts.atomratio == "mdlFILM"){
    cout << "Gear(5) RAND=10 " << opts.potential << " T=300.0 STEP=(0,1000100)" << endl;
    cout << "BLOCK FIXANGLE W=INF Q=INF LOGSTEP=10 PRINTVELOCITY dt=0.1fs scratch(100)" << endl << endl;
    cout << title << endl << endl;
    cout << leng[0] << " " << tilt[0][1] << " " << tilt[0][2]*3 << opts.pbc_x << endl;
    cout << tilt[1][0] << " " << leng[1] << " " << tilt[1][2]*3 << opts.pbc_y << endl;
    cout << tilt[2][0]*3 << " " << tilt[2][1]*3 << " " << leng[2]*3 << opts.pbc_z << endl << endl;
    for(int j=0; j<nx; j++) for(int k=0; k<ny; k++) for(int l=0; l<nz; l++)
    for(int i=0;i<all_atom;i++){
      double ax = (zero_atom[i].x+j)/nx, ay = (zero_atom[i].y+k)/ny, az = (zero_atom[i].z+l)/(nz*3)+1.0/3.0;
      if(shape(opts.shape, ax, ay, az)){
        cout << zero_atom[i].name << " " << ax << " " << ay << " " << az;
        if(zero_atom[i].option != "") cout << zero_atom[i].option;
        cout << endl;
      }
    }
  }
  if(opts.atomratio == "mdlWIRE"){
    cout << "Gear(5) RAND=10 " << opts.potential << " T=300.0 STEP=(0,1000100)" << endl;
    cout << "BLOCK FIXANGLE W=INF Q=INF LOGSTEP=10 PRINTVELOCITY dt=0.1fs scratch(100)" << endl << endl;
    cout << title << endl << endl;
    cout << leng[0] << " " << tilt[0][1]*3 << " " << tilt[0][2]*3 << opts.pbc_x << endl;
    cout << tilt[1][0]*3 << " " << leng[1]*3 << " " << tilt[1][2]*3 << opts.pbc_y << endl;
    cout << tilt[2][0]*3 << " " << tilt[2][1]*3 << " " << leng[2]*3 << opts.pbc_z << endl << endl;
    for(int j=0; j<nx; j++) for(int k=0; k<ny; k++) for(int l=0; l<nz; l++)
    for(int i=0;i<all_atom;i++){
      double ax = (zero_atom[i].x+j)/nx, ay = (zero_atom[i].y+k)/(ny*3)+1.0/3.0, az = (zero_atom[i].z+l)/(nz*3)+1.0/3.0;
      if(shape(opts.shape, ax, ay, az)){
        cout << zero_atom[i].name << " " << ax << " " << ay << " " << az;
        if(zero_atom[i].option != "") cout << zero_atom[i].option;
        cout << endl;
      }
    }
  }
  else if(opts.atomratio == "mdlFIX"){
    cout << "OPT(P) RAND=10 " << opts.potential << " T=300.0 STEP=(0,20000)" << endl;
    cout << "BLOCK W=INF FIXANGLE Q=1.0" << endl << endl;
    cout << title << endl << endl;
    cout << leng[0] << " " << tilt[0][1] << " " << tilt[0][2] << opts.pbc_x << endl;
    cout << tilt[1][0] << " " << leng[1] << " " << tilt[1][2] << opts.pbc_y << endl;
    cout << tilt[2][0] << " " << tilt[2][1] << " " << leng[2] << opts.pbc_z << endl << endl;
    for(int j=0; j<nx; j++) for(int k=0; k<ny; k++) for(int l=0; l<nz; l++)
    for(int i=0;i<all_atom;i++){
      double ax = (zero_atom[i].x+j)/nx, ay = (zero_atom[i].y+k)/ny, az = (zero_atom[i].z+l)/nz;
      if(shape(opts.shape, ax, ay, az)){
        cout << zero_atom[i].name << " ";
        cout << ax << " " << ay << " " << az;
        if(az<FIX_ratio) cout << " F";
        if(zero_atom[i].option != "") cout << zero_atom[i].option;
        cout << endl;
      }
    }
  }
  else if(opts.atomratio == "mdlXO2"){
    cout << "OPT(P) RAND=10 " << opts.potential << " T=300.0 STEP=(0,20000)" << endl;
    cout << "BLOCK W=INF FIXANGLE Q=1.0" << endl << endl;
    cout << title << endl << endl;
    cout << leng[0] << " " << tilt[0][1] << " " << tilt[0][2] << opts.pbc_x << endl;
    cout << tilt[1][0] << " " << leng[1] << " " << tilt[1][2] << opts.pbc_y << endl;
    cout << tilt[2][0] << " " << tilt[2][1] << " " << leng[2] << opts.pbc_z << endl << endl;
    for(int j=0; j<nx; j++) for(int k=0; k<ny; k++) for(int l=0; l<nz; l++)
    for(int i=0;i<all_atom;i++){
      double ax = (zero_atom[i].x+j)/nx, ay = (zero_atom[i].y+k)/ny, az = (zero_atom[i].z+l)/nz;
      if(shape(opts.shape, ax, ay, az)){
        if(az>XO2_ratio && zero_atom[i].name == "Si") cout << XO2_atom << " ";
        else cout << zero_atom[i].name << " ";
        cout << ax << " " << ay << " " << az;
        if(zero_atom[i].option != "") cout << zero_atom[i].option;
        cout << endl;
      }
    }
  }

  else if(opts.atomratio == "mdlCLEAN"){ //takizawa
    int ok;
    int group[all_atom];
    int group_now    = 1;
    int group_max    = 0;
    int cut          = pow(4.0, 2);
    for(int i=0; i<all_atom; i++){
      group[i] = 0;
     }
    for(int i=0; i<all_atom; i++){
     for(int j=i+1; j<all_atom; j++){
       ok = 0;
       double dx = fabs(zero_atom[i].x - zero_atom[j].x) * leng[0];
       double dy = fabs(zero_atom[i].y - zero_atom[j].y) * leng[1];
       double dz = fabs(zero_atom[i].z - zero_atom[j].z) * leng[2];

       double r   = pow(dx, 2.0) + pow(dy, 2.0) + pow(dz, 2.0);
       double rx  = pow(leng[0] - dx, 2.0) + pow(dy, 2.0) + pow(dz, 2.0);
       double ry  = pow(dx, 2.0) + pow(leng[1] - dy, 2.0) + pow(dz, 2.0);
       double rxy  = pow(leng[0] - dx, 2.0) + pow(leng[1] - dy, 2.0) + pow(dz, 2.0);

       if(dx<=leng[0]/2.0 && dy<=leng[1]/2.0 && r<cut && r!=0) ok = 1;
       if(dx>leng[0]/2.0 && dy<=leng[1]/2.0 && rx<cut) ok = 1;
       if(dx<=leng[0]/2.0 && dy>leng[1]/2.0 && ry<cut) ok = 1;
       if(dx>leng[0]/2.0 && dy>leng[1]/2.0 && rxy<cut) ok = 1;

       if(ok == 1){
           if(group[i] == 0 && group[j] == 0){
               group[i] = group_now;
               group[j] = group_now;
               ++group_now;
           }else if(group[i] == 0 && group[j] != 0){
               group[i] = group[j];
           }else if(group[i] != 0 && group[j] == 0){
               group[j] = group[i];
           }else if(group[i] != 0 && group[j] != 0 && group[i] < group[j]){
               for(int k=0; k<all_atom; k++){
                    if(group[k] == group[i]) group[k] = group[j];
                   }
           }else if(group[i] != 0 && group[j] != 0 && group[i] > group[j]){
               for(int k=0; k<all_atom; k++){
                    if(group[k] == group[j]) group[k] = group[i];
                   }
              }
         }
       }
     }
    int group_number[group_now];
    for(int i=0; i<group_now; i++){
      group_number[i] = 0;
     }
    int group_number_max = 0;
    for(int i=0; i<all_atom; i++){
      if(group[i] != 0){
         group_number[group[i]] = group_number[group[i]] + 1;
        }
     }
    for(int i=1; i<group_now; i++){
      if(group_number_max < group_number[i]){
         group_number_max = group_number[i];
         group_max        = i;
        }
     }
    cerr << "max group have " << group_number_max << endl;

    cout << "Gear(5) RAND=10 " << opts.potential << " T=300.0 STEP=(0,1000100)" << endl;
    cout << "BLOCK FIXANGLE W=INF Q=INF LOGSTEP=10 PRINTVELOCITY dt=0.1fs scratch(100)" << endl << endl;
    cout << title << endl << endl;
    cout << leng[0] << " " << tilt[0][1] << " " << tilt[0][2] << opts.pbc_x << endl;
    cout << tilt[1][0] << " " << leng[1] << " " << tilt[1][2] << opts.pbc_y << endl;
    cout << tilt[2][0] << " " << tilt[2][1] << " " << leng[2] << opts.pbc_z << endl << endl;
    for(int j=0; j<nx; j++) for(int k=0; k<ny; k++) for(int l=0; l<nz; l++)
    for(int i=0;i<all_atom;i++){
      double ax = (zero_atom[i].x+j)/nx, ay = (zero_atom[i].y+k)/ny, az = (zero_atom[i].z+l)/nz;
      if(shape(opts.shape, ax, ay, az) && group[i] == group_max){
        cout << zero_atom[i].name << " " << ax << " " << ay << " " << az;
        if(zero_atom[i].option != "") cout << zero_atom[i].option;
        cout << endl;
      }
    }
  }
  else if(opts.atomratio == "mdlZRES"){ //takizawa
    double max = 0;
    double min = 1;
    for(int i=0;i<all_atom;i++){
       if(max < zero_atom[i].z) max = zero_atom[i].z;
       if(min > zero_atom[i].z) min = zero_atom[i].z;
     }
    cout << "Gear(5) RAND=10 " << opts.potential << " T=300.0 STEP=(0,1000100)" << endl;
    cout << "BLOCK FIXANGLE W=INF Q=INF LOGSTEP=10 PRINTVELOCITY dt=0.1fs scratch(100)" << endl << endl;
    cout << title << endl << endl;
    cout << leng[0] << " " << tilt[0][1] << " " << tilt[0][2] << opts.pbc_x << endl;
    cout << tilt[1][0] << " " << leng[1] << " " << tilt[1][2] << opts.pbc_y << endl;
    cout << tilt[2][0] << " " << tilt[2][1] << " " << leng[2] * (max - min) * FILM_ratio << opts.pbc_z << endl << endl;
    for(int j=0; j<nx; j++) for(int k=0; k<ny; k++) for(int l=0; l<nz; l++)
    for(int i=0;i<all_atom;i++){
      double ax = (zero_atom[i].x+j)/nx, ay = (zero_atom[i].y+k)/ny; 
      double az = ((zero_atom[i].z+l)/nz - min) / ((max - min) * FILM_ratio) + (FILM_ratio - 1.0 ) / (FILM_ratio * 2.0);
      if(shape(opts.shape, ax, ay, az)){
        cout << zero_atom[i].name << " " << ax << " " << ay << " " << az;
        if(zero_atom[i].option != "") cout << zero_atom[i].option;
        cout << endl;
      }
    }
  }

  else if(opts.atomratio == "mdlGROUP"){ //takizawa
    int ok;
    double dx, dy, dz, r, rx, ry, rz, rxy, ryz, rxz, rxyz;
    int group[all_atom];
    int group_now    = 1;
    double cutG         = pow(3.0, 2); //2.35
    for(int i=0; i<all_atom; i++){
      group[i] = 0;
    }
    for(int i=0; i<all_atom; i++){
     if(zero_atom[i].name != "Ge") continue;
     for(int j=i+1; j<all_atom; j++){
       if(zero_atom[j].name != "Ge") continue;
       ok = 0;
       dx = fabs(zero_atom[i].x - zero_atom[j].x) * leng[0];
       dy = fabs(zero_atom[i].y - zero_atom[j].y) * leng[1];
       dz = fabs(zero_atom[i].z - zero_atom[j].z) * leng[2];

       r    = pow(dx, 2.0) + pow(dy, 2.0) + pow(dz, 2.0);
       rx   = pow(leng[0] - dx, 2.0) + pow(dy, 2.0) + pow(dz, 2.0);
       ry   = pow(dx, 2.0) + pow(leng[1] - dy, 2.0) + pow(dz, 2.0);
       rz   = pow(dx, 2.0) + pow(dy, 2.0) + pow(leng[2] - dz, 2.0);
       rxy  = pow(leng[0] - dx, 2.0) + pow(leng[1] - dy, 2.0) + pow(dz, 2.0);
       ryz  = pow(dx, 2.0) + pow(leng[1] - dy, 2.0) + pow(leng[2] - dz, 2.0);
       rxz  = pow(leng[0] - dx, 2.0) + pow(dy, 2.0) + pow(leng[2] - dz, 2.0);
       rxyz = pow(leng[0] - dx, 2.0) + pow(leng[1] - dy, 2.0) + pow(leng[2] - dz, 2.0);

       if(r < cutG) ok = 1;
       else if(rz < cutG && nz==1) ok = 1;
       else if(ry < cutG && ny==1) ok = 1;
       else if(ryz < cutG && ny==1 && nz==1) ok = 1;
       else if(rx < cutG && nx==1) ok = 1;
       else if(rxz  <cutG && nx==1 && nz==1) ok = 1;
       else if(rxy  <cutG && nx==1 && ny==1) ok = 1;
       else if(rxyz <cutG && nx==1 && ny==1 && nz==1) ok = 1; 

/*
dx <= leng[0]/2.0 && dy <= leng[1]/2.0 && dz <= leng[2]/2.0 && 
dx <= leng[0]/2.0 && dy <= leng[1]/2.0 && dz >  leng[2]/2.0 && 
dx <= leng[0]/2.0 && dy >  leng[1]/2.0 && dz <= leng[2]/2.0 && 
dx <= leng[0]/2.0 && dy >  leng[1]/2.0 && dz >  leng[2]/2.0 && 
dx >  leng[0]/2.0 && dy <= leng[1]/2.0 && dz <= leng[2]/2.0 && 
dx >  leng[0]/2.0 && dy <= leng[1]/2.0 && dz >  leng[2]/2.0 && 
dx >  leng[0]/2.0 && dy >  leng[1]/2.0 && dz <= leng[2]/2.0 && 
dx >  leng[0]/2.0 && dy >  leng[1]/2.0 && dz >  leng[2]/2.0 &&
*/

       if(ok == 1){
           if(group[i] == 0 && group[j] == 0){
               group[i] = group_now;
               group[j] = group_now;
               ++group_now;
           }else if(group[i] == 0 && group[j] != 0){
               group[i] = group[j];
           }else if(group[i] != 0 && group[j] == 0){
               group[j] = group[i];
           }else if(group[i] != 0 && group[j] != 0 && group[i] < group[j]){
               for(int k=0; k<all_atom; k++){
                    if(group[k] == group[i]) group[k] = group[j];
                   }
           }else if(group[i] != 0 && group[j] != 0 && group[i] > group[j]){
               for(int k=0; k<all_atom; k++){
                    if(group[k] == group[j]) group[k] = group[i];
               }
           }
       }
     }
    }
    int renumber[all_atom];
    int numbercount[all_atom];
    int number = 1;
    for(int i=0; i<=all_atom; i++){
     renumber[i] = 0;
    }
    for(int i=0; i<=all_atom; i++){
     numbercount[i] = 0;
    }
    cout << "Gear(5) RAND=10 " << opts.potential << " T=300.0 STEP=(0,1000100)" << endl;
    cout << "BLOCK FIXANGLE W=INF Q=INF LOGSTEP=10 PRINTVELOCITY dt=0.1fs scratch(100)" << endl << endl;
    cout << title << endl << endl;
    cout << leng[0] << " " << tilt[0][1] << " " << tilt[0][2] << opts.pbc_x << endl;
    cout << tilt[1][0] << " " << leng[1] << " " << tilt[1][2] << opts.pbc_y << endl;
    cout << tilt[2][0] << " " << tilt[2][1] << " " << leng[2] << opts.pbc_z << endl << endl;
    for(int i=0;i<all_atom;i++){
      double ax = zero_atom[i].x, ay = zero_atom[i].y, az = zero_atom[i].z;
      if(shape(opts.shape, ax, ay, az)){
        cout << zero_atom[i].name << " " << ax << " " << ay << " " << az;
      	if(zero_atom[i].name == "Ge"){
	  if(group[i]==0){
	    cout << " G " << number;
            ++numbercount[number];
	    ++number;
	  }
	  else{
	    if(renumber[group[i]] != 0){
	      cout << " G " << renumber[group[i]];
              ++numbercount[renumber[group[i]]];
	    }
	    else{
	      cout << " G " << number;
	      renumber[group[i]] = number;
              ++numbercount[number];
	      ++number;
	    }
	  }
      	}
      	else{
         	cout << " F";
      	}
        cout << endl;
      }
    }
    stringstream ready;
    ready << opts.atomcomp << ".cdist"; 
    string name = ready.str();
    ofstream dist(name.c_str());
    for(int i=1; i<number; i++){
      dist << i << " " << numbercount[i] << endl;
    }
  }

  else if(opts.atomratio == "mdlINTERFACE"){ //takizawa
    double aisle_xy = 0.0; 
    double aisle_z = 10.0;   
    double aisle_ratio = 6; 

    double Amax=0.0;
    double Bmin=0.7;
    for(int i=0;i<all_atom;i++){
	if(Amax<zero_atom[i].z && zero_atom[i].z<0.7 && zero_atom[i].name=="Si") Amax=zero_atom[i].z;
	if(Bmin>zero_atom[i].z && zero_atom[i].z>0.3 && zero_atom[i].name!="Si" && zero_atom[i].name!="O") Bmin=zero_atom[i].z;
    }
    double INTERFACE=leng[2]*(Amax+Bmin)/2;

    int OK[all_atom];
    int count = 0;
    int Si_num = 0, O_numA = 0, K_num = 0, O_numB = 0;
    for(int i=0;i<all_atom;i++) {
	OK[i] = 1;
	count++;
	if(zero_atom[i].z*leng[2] > INTERFACE-aisle_z/2.0 && zero_atom[i].z*leng[2] < INTERFACE+aisle_z/2.0 && 
	   (zero_atom[i].x*leng[0] > aisle_xy || zero_atom[i].y*leng[1] > aisle_xy) && count%10 >= aisle_ratio) {
		OK[i] = 0;
		if(zero_atom[i].name == "Si") Si_num++;
		if(zero_atom[i].name == "O" && zero_atom[i].z*leng[2] <= INTERFACE) O_numA++;
		if(zero_atom[i].name != "Si" && zero_atom[i].name != "O") K_num++;
		if(zero_atom[i].name == "O" && zero_atom[i].z*leng[2] >= INTERFACE) O_numB++;
	}
    }

    int charge = Si_num*(int)atom_CIM[1] + (O_numA+O_numB)*(int)atom_CIM[2] + K_num*(int)atom_CIM[0];
    cerr << "Delete Area: xy= leng[xy]-" << aisle_xy << " z= "<< aisle_z << endl;
    cerr << "Delete Charge: " << charge << endl;
    cerr << "Metal:" << K_num << " Si:" << Si_num << " O:" << O_numA+O_numB << endl;
    int Si_delete = 0, O_deleteA = 0, K_delete = 0, O_deleteB = 0;
    if(charge != 0){
	int CIM_K = (int)atom_CIM[0], CIM_Si = (int)atom_CIM[1], CIM_O = -1 * (int)atom_CIM[2];
	while((O_numA+O_deleteA)%2 != 0) O_deleteA++;
	if     (Si_num*2 > (O_numA+O_deleteA))O_deleteA = Si_num*2 - O_numA;
	else if(Si_num*2 < (O_numA+O_deleteA))Si_delete =(O_numA+O_deleteA)/2-Si_num;		

	while((CIM_K*(K_num + K_delete)) % CIM_O != 0) K_delete++; 
	while((CIM_O*(O_numB + O_deleteB)) % CIM_K != 0) O_deleteB++;
	if     ((K_num + K_delete)*CIM_K > (O_numB + O_deleteB)*CIM_O) O_deleteB = CIM_K*(K_num + K_delete)/CIM_O - O_numB;
	else if((K_num + K_delete)*CIM_K < (O_numB + O_deleteB)*CIM_O)  K_delete = CIM_O*(O_numB + O_deleteB)/CIM_K - K_num;

	if((Si_num+Si_delete)*(int)atom_CIM[1] + (O_numA+O_numB+O_deleteA+O_deleteB)*(int)atom_CIM[2] + (K_num+K_delete)*(int)atom_CIM[0] != 0){
		cerr << "Program Error!: charge neutrality[delete]" << endl;
		return 1;
	}
	cerr << "Additional Delete Charge:" << endl;
	cerr << "Metal:" << K_delete << " Si:" << Si_delete << " O:" << O_deleteA+O_deleteB << endl;

	double aisle_d = 0.001;
	double aisle_dd = 0.001;
	while(Si_delete > 0 || O_deleteA > 0 || K_delete > 0 || O_deleteB > 0){
		for(int i=0;i<all_atom;i++) {
			if(zero_atom[i].x*leng[0] > aisle_xy && zero_atom[i].y*leng[1] > aisle_xy && OK[i] == 1){
				if(zero_atom[i].z*leng[2] > INTERFACE - aisle_z/2.0 - aisle_d && zero_atom[i].z*leng[2] < INTERFACE - aisle_z/2.0){
					if(zero_atom[i].name == "Si" && Si_delete > 0) OK[i] = 0, Si_delete--, Si_num++;
					if(zero_atom[i].name == "O" && O_deleteA > 0) OK[i] = 0, O_deleteA--, O_numA++;
				}
				if(zero_atom[i].z*leng[2] < INTERFACE + aisle_z/2.0 + aisle_d && zero_atom[i].z*leng[2] > INTERFACE + aisle_z/2.0){
					if(zero_atom[i].name != "Si" && zero_atom[i].name != "O" && K_delete > 0) OK[i] = 0, K_delete--, K_num++;
					if(zero_atom[i].name == "O" && O_deleteB > 0) OK[i] = 0, O_deleteB--, O_numB++;
				}		
			
			}
		}
		aisle_d = aisle_d + aisle_dd;
		if(aisle_d > aisle_dd * 100000){
			cerr << "Program Error!: Infinity Roop" << endl;
			return 1;
		}
	}
	cerr << "d: " << aisle_d << endl;
    }

    charge = Si_num*(int)atom_CIM[1] + (O_numA+O_numB)*(int)atom_CIM[2] + K_num*(int)atom_CIM[0];
    if(charge != 0){
	cerr << "Program Error!: charge neutrality[delete2]" << endl;
	return 1;
    }

    int Si_all = 0, O_all = 0, K_all = 0;
    for(int i=0;i<all_atom;i++){
	if(zero_atom[i].name == "Si" && OK[i] == 1) Si_all++;
	if(zero_atom[i].name == "O"  && OK[i] == 1)  O_all++;
	if(zero_atom[i].name != "Si" && zero_atom[i].name != "O" && OK[i] == 1) K_all++;	
    }
    if(Si_all*(int)atom_CIM[1] + O_all*(int)atom_CIM[2] + K_all*(int)atom_CIM[0] != 0){
	cerr << "Program Error!: charge neutrality[all]" << endl;
	return 1;
    }

    cout << "Gear(5) RAND=10 " << opts.potential << " T=300.0 STEP=(0,1000100) INTERFACE=" << INTERFACE << endl;
    cout << "BLOCK FIXANGLE W=INF Q=INF LOGSTEP=10 PRINTVELOCITY dt=0.1fs scratch(100)" << endl << endl;
    cout << title << endl << endl;
    cout << leng[0] << " " << tilt[0][1] << " " << tilt[0][2] << opts.pbc_x << endl;
    cout << tilt[1][0] << " " << leng[1] << " " << tilt[1][2] << opts.pbc_y << endl;
    cout << tilt[2][0] << " " << tilt[2][1] << " " << leng[2] << opts.pbc_z << endl << endl;
    for(int j=0; j<nx; j++) for(int k=0; k<ny; k++) for(int l=0; l<nz; l++)
    for(int i=0;i<all_atom;i++){
      double ax = (zero_atom[i].x+j)/nx, ay = (zero_atom[i].y+k)/ny, az = (zero_atom[i].z+l)/nz;
      if(shape(opts.shape, ax, ay, az) && OK[i]==1){
        cout << zero_atom[i].name << " " << ax << " " << ay << " " << az;
        if(zero_atom[i].option != "") cout << zero_atom[i].option;
        cout << endl;
      }
    }
  }

  else if(opts.atomratio == "lmp"){
    cout << title << endl << endl;
    cout << " " << all_atom*nx*ny*nz << " atoms" << endl;
    cout << " " << all_type << " atom types" << endl;
    cout << " 0.00 " << leng[0] << " xlo xhi" << endl;
    cout << " 0.00 " << leng[1] << " ylo yhi" << endl;
    cout << " 0.00 " << leng[2] << " zlo zhi" << endl;
    if(tilt[0][1] || tilt[0][2] || tilt[1][0] || tilt[1][2] || tilt[2][0] || tilt[2][1])
      cout << " " << tilt[0][1] + tilt[1][0] << " " << tilt[0][2] + tilt[2][0] << " " << tilt[1][2] + tilt[2][1] << " xy xz yz" << endl;
    cout << endl << " Masses" << endl<< endl;
    for(int i=0;i<all_type;i++) cout << i+1 << " " << atom_mass[i] << endl;
    cout << endl << " Atoms" << endl << endl;
    for(int j=0; j<nx; j++) for(int k=0; k<ny; k++) for(int l=0; l<nz; l++)
    for(int i=0;i<all_atom;i++){


      double ax = (zero_atom[i].x+j)/nx, ay = (zero_atom[i].y+k)/ny, az = (zero_atom[i].z+l)/nz;
      if(shape(opts.shape, ax, ay, az)){
        cout << i+1+(j*ny*nz*all_atom)+(k*nz*all_atom)+(l*all_atom) << " ";
        for(int j=0;j<all_type;j++) if(zero_atom[i].name == atom_type[j]) cout << j+1 << " ";
        cout << ax*leng[0] << " " << ay*leng[1] << " " << az*leng[2] << endl;
      }
    }
  }
  else if(opts.atomratio == "lmpCIM"){
    cout << title << endl << endl;
    cout << " " << all_atom*nx*ny*nz << " atoms" << endl;
    cout << " " << all_type << " atom types" << endl;
    cout << " 0.00 " << leng[0] << " xlo xhi" << endl;
    cout << " 0.00 " << leng[1] << " ylo yhi" << endl;
    cout << " 0.00 " << leng[2] << " zlo zhi" << endl;
    if(tilt[0][1] || tilt[0][2] || tilt[1][0] || tilt[1][2] || tilt[2][0] || tilt[2][1])
      cout << " " << tilt[0][1] + tilt[1][0] << " " << tilt[0][2] + tilt[2][0] << " " << tilt[1][2] + tilt[2][1] << " xy xz yz" << endl;
    cout << endl << " Masses" << endl<< endl;
    for(int i=0;i<all_type;i++) cout << i+1 << " " << atom_mass[i] << endl;
    cout << endl << " Atoms" << endl << endl;
    for(int j=0; j<nx; j++) for(int k=0; k<ny; k++) for(int l=0; l<nz; l++)
    for(int i=0;i<all_atom;i++){
      double ax = (zero_atom[i].x+j)/nx, ay = (zero_atom[i].y+k)/ny, az = (zero_atom[i].z+l)/nz;
      if(shape(opts.shape, ax, ay, az)){
        cout << i+1+(j*ny*nz*all_atom)+(k*nz*all_atom)+(l*all_atom) << " ";
        for(int j=0;j<all_type;j++) if(zero_atom[i].name == atom_type[j]) cout << j+1 << " " << atom_CIM[j] << " ";
        cout << ax*leng[0] << " " << ay*leng[1] << " " << az*leng[2] << endl;
      }
    }
  }
  else if(opts.atomratio == "lmpCMAS"){
    cout << title << endl << endl;
    cout << " " << all_atom*nx*ny*nz << " atoms" << endl;
    cout << " " << all_type << " atom types" << endl;
    cout << " 0.00 " << leng[0] << " xlo xhi" << endl;
    cout << " 0.00 " << leng[1] << " ylo yhi" << endl;
    cout << " 0.00 " << leng[2] << " zlo zhi" << endl;
    if(tilt[0][1] || tilt[0][2] || tilt[1][0] || tilt[1][2] || tilt[2][0] || tilt[2][1])
      cout << " " << tilt[0][1] + tilt[1][0] << " " << tilt[0][2] + tilt[2][0] << " " << tilt[1][2] + tilt[2][1] << " xy xz yz" << endl;
    cout << endl << " Masses" << endl<< endl;
    for(int i=0;i<all_type;i++) cout << i+1 << " " << atom_mass[i] << endl;
    cout << endl << " Atoms" << endl << endl;
    for(int j=0; j<nx; j++) for(int k=0; k<ny; k++) for(int l=0; l<nz; l++)
    for(int i=0;i<all_atom;i++){
      double ax = (zero_atom[i].x+j)/nx, ay = (zero_atom[i].y+k)/ny, az = (zero_atom[i].z+l)/nz;
      if(shape(opts.shape, ax, ay, az)){
        cout << i+1+(j*ny*nz*all_atom)+(k*nz*all_atom)+(l*all_atom) << " ";
        for(int j=0;j<all_type;j++) if(zero_atom[i].name == atom_type[j]) cout << j+1 << " " << atom_CMAS[j] << " ";
        cout << ax*leng[0] << " " << ay*leng[1] << " " << az*leng[2] << endl;
      }
    }
  }
  else if(opts.atomratio == "lmpAr"){
    cout << title << endl << endl;
    cout << " " << all_atom*nx*ny*nz << " atoms" << endl;
    cout << " " << all_type+1 << " atom types" << endl;
    cout << " 0.00 " << leng[0] << " xlo xhi" << endl;
    cout << " 0.00 " << leng[1] << " ylo yhi" << endl;
    cout << " 0.00 " << leng[2] << " zlo zhi" << endl;
    if(tilt[0][1] || tilt[0][2] || tilt[1][0] || tilt[1][2] || tilt[2][0] || tilt[2][1])
      cout << " " << tilt[0][1] + tilt[1][0] << " " << tilt[0][2] + tilt[2][0] << " " << tilt[1][2] + tilt[2][1] << " xy xz yz" << endl;
    cout << endl << " Masses" << endl<< endl;
    for(int i=0;i<all_type;i++) cout << i+1 << " " << atom_mass[i] << endl;
    cout << all_type+1 << " 39.948" << endl;
    cout << endl << " Atoms" << endl << endl;
    for(int j=0; j<nx; j++) for(int k=0; k<ny; k++) for(int l=0; l<nz; l++)
    for(int i=0;i<all_atom;i++){
      double ax = (zero_atom[i].x+j)/nx, ay = (zero_atom[i].y+k)/ny, az = (zero_atom[i].z+l)/nz;
      if(shape(opts.shape, ax, ay, az)){
        cout << i+1+(j*ny*nz*all_atom)+(k*nz*all_atom)+(l*all_atom) << " ";
        for(int j=0;j<all_type;j++) if(zero_atom[i].name == atom_type[j]) cout << j+1 << " ";
        cout << ax*leng[0] << " " << ay*leng[1] << " " << az*leng[2] << endl;
      }
    }
  }
  else if(opts.atomratio == "lmpLAY"){
    cout << title << endl << endl;
    cout << " " << all_atom*nx*ny*nz << " atoms" << endl;
    cout << " " << all_type << " atom types" << endl;
    cout << " 0.00 " << leng[0] << " xlo xhi" << endl;
    cout << " 0.00 " << leng[1] << " ylo yhi" << endl;
    cout << " 0.00 " << leng[2] << " zlo zhi" << endl;
    if(tilt[0][1] || tilt[0][2] || tilt[1][0] || tilt[1][2] || tilt[2][0] || tilt[2][1])
      cout << " " << tilt[0][1] + tilt[1][0] << " " << tilt[0][2] + tilt[2][0] << " " << tilt[1][2] + tilt[2][1] << " xy xz yz" << endl;
    cout << endl << " Masses" << endl<< endl;
    for(int i=0;i<all_type;i++) cout << i+1 << " " << atom_mass[i] << endl;
    cout << endl << " Atoms" << endl << endl;
    for(int j=0; j<nx; j++) for(int k=0; k<ny; k++) for(int l=0; l<nz; l++)
    for(int i=0;i<all_atom;i++){
      double ax = (zero_atom[i].x+j)/nx, ay = (zero_atom[i].y+k)/ny, az = (zero_atom[i].z+l)/nz;
      if(shape(opts.shape, ax, ay, az)){
        cout << i+1+(j*ny*nz*all_atom)+(k*nz*all_atom)+(l*all_atom) << " " << zero_atom[i].lay << " ";
        for(int j=0;j<all_type;j++) if(zero_atom[i].name == atom_type[j]) cout << j+1 << " ";
        cout << ax*leng[0] << " " << ay*leng[1] << " " << az*leng[2] << endl;
      }
    }
  }
  else if(opts.atomratio == "xyz"){
    cout << all_atom*nx*ny*nz << endl;
    cout << title << endl;
    for(int j=0; j<nx; j++) for(int k=0; k<ny; k++) for(int l=0; l<nz; l++)
    for(int i=0;i<all_atom;i++){
      double ax = (zero_atom[i].x+j)/nx, ay = (zero_atom[i].y+k)/ny, az = (zero_atom[i].z+l)/nz;
      if(shape(opts.shape, ax, ay, az)){
        cout << zero_atom[i].name << " " << ax*leng[0] << " " << ay*leng[1] << " " << az*leng[2] << endl;
      }
    }
  }
  else if(opts.atomratio == "vasp"){
    cout << title << endl << "1.0" << endl;
    cout << fixed << setprecision(10) << setw(20) << leng[0] << " " << setw(20) << tilt[0][1] << " " << setw(20) << tilt[0][2] << endl;
    cout << fixed << setprecision(10) << setw(20) << tilt[1][0] << " " << setw(20) << leng[1] << " " << setw(20) << tilt[1][2] << endl;
    cout << fixed << setprecision(10) << setw(20) << tilt[2][0] << " " << setw(20) << tilt[2][1] << " " << setw(20) << leng[2] << endl;
    for(int m=0; m<all_type; m++) cout << setw(5) << atom_type[m];
    cout << endl;
    for(int m=0; m<all_type; m++) cout << setprecision(10) << setw(5) << atom_num[m]*nx*ny*nz;
    cout << endl << "Direct" << endl;
    for(int m=0; m<all_type; m++)
    for(int j=0; j<nx; j++) for(int k=0; k<ny; k++) for(int l=0; l<nz; l++)
    for(int i=0;i<all_atom;i++){
      double ax = (zero_atom[i].x+j)/nx, ay = (zero_atom[i].y+k)/ny, az = (zero_atom[i].z+l)/nz;
      if(shape(opts.shape, ax, ay, az) && zero_atom[i].name == atom_type[m]){
        cout << setprecision(10) << setw(20) << ax << " " << setw(20) << ay << " " << setw(20) << az << endl;
      }
    }
  }
  else if(opts.atomratio == "alm"){
    double bohr = 0.5291772106712;
    double leng_max = leng[0];
    if(leng_max < leng[1]) leng_max = leng[1];
    if(leng_max < leng[2]) leng_max = leng[2];
/*
    cout << "&general" << endl;
    cout << "  PREFIX = " << atom_name << endl;
    cout << "  MODE = suggest" << endl;
    cout << "  NAT = " << all_atom << "; NKD = " << all_type << endl;
    cout << "  KD =";
    for(int j=0;j<all_type;j++) cout << " " << atom_type[j];
    cout << endl;
    cout << "  PERIODIC =";
    if(opts.pbc_x == " N") cout << " 0"; else cout << " 1";
    if(opts.pbc_y == " N") cout << " 0"; else cout << " 1";


    if(opts.pbc_z == " N") cout << " 0"; else cout << " 1";
    cout << endl;
    cout << "/" << endl << endl;
    cout << "&interaction" << endl;
    cout << "  NORDER = 2" << endl;
    cout << "/" << endl << endl;
*/
    cout << "&cell" << endl;
    cout << "  " << leng_max / bohr << endl;
    cout << "  " << leng[0]/leng_max << " " << tilt[0][1]/leng_max << " " << tilt[0][2]/leng_max << endl;
    cout << "  " << tilt[1][0]/leng_max << " " << leng[1]/leng_max << " " << tilt[1][2]/leng_max << endl;
    cout << "  " << tilt[2][0]/leng_max << " " << tilt[2][1]/leng_max << " " << leng[2]/leng_max << endl;
    cout << "/" << endl << endl;
    cout << "&position" << endl;
    for(int j=0; j<nx; j++) for(int k=0; k<ny; k++) for(int l=0; l<nz; l++)
    for(int i=0;i<all_atom;i++){
      double ax = (zero_atom[i].x+j)/nx, ay = (zero_atom[i].y+k)/ny, az = (zero_atom[i].z+l)/nz;
      if(shape(opts.shape, ax, ay, az)){
        cout << "  ";
        for(int j=0;j<all_type;j++) if(zero_atom[i].name == atom_type[j]) cout << j+1;
        cout << " " << ax << " " << ay << " " << az << endl;
      }
    }
  }

  delete[] zero_atom;
}

