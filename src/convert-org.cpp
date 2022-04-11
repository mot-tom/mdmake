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

#define MAXATOMTYPE 32
#define MAXBONDTYPE 32
#define MAXANGLTYPE 32
#define MAXTORSTYPE 32
#define MAXIMPRTYPE 32

struct info        // information of atom
{
  string name;// atom name
  int atom_type_id;
  double x, y, z;  // coordinate of atom
  int nei_atom[4];
  double nei_bond[4];
  int nei_num;
  double max_bond;
};

struct bond
{
  string name;
  int bond_type_id;
  int atom[2];
};

struct angle
{
  string name;
  int angle_type_id;
  int atom[3];
};

struct torsion
{
  string name;
  int torsion_type_id;
  int atom[4];
};

struct improper
{
  string name;
  int improper_type_id;
  int atom[4];
};

int convert_org(int argc, char **argv, optpara opts){
  const int max_atom = 100000;         // maximum number of atom
  const int max_bond = 50000;         // maximum number of bond
  const int max_angle = 50000;         // maximum number of angle
  const int max_torsion = 50000;         // maximum number of torsion
  const int max_improper = 50000;         // maximum number of angle

  double leng_h = 0.0;    // length of each side of crystal
  double leng_l = 0.0;    // length of each side of crystal 

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
    if(temp_line.find("%chk=") != string::npos && line == 0){file_type = "com";}
    line++;
  }
  file_fin.close();

  if(file_type == "com") cerr << "\"" << opts.atomcomp << "\" is input file for gaussian." << endl;
  else{
    cerr << "\"" << opts.atomcomp << "\" is wrong file type!" << endl;
    return 1;
  }

  if(opts.atomratio == "lmp") cerr << "\"" << opts.atomcomp << "\" convert to lmp file." << endl;
  else{
    cerr << "Output type of \"" << opts.atomratio << "\" is wrong!" << endl;
    return 1;
  }


  file_fin.open(opts.atomcomp.c_str());
  int blank = 0;
  string title, first;
  info* zero_atom = new info[max_atom];
  int all_atom = 0;
  if(file_type == "com"){
    while(getline(file_fin, temp_line)){
      istringstream line_iss(temp_line.c_str());

      if(blank == 1 && temp_line != "")
      {
        title = temp_line;
      }

      if(blank == 2 && temp_line != "")
      {
        line_iss >> first;
        if(first == "0") continue;
        if(first == "Si" || first == "C" || first == "O" || first == "H")
        {
          zero_atom[all_atom].name = first;
          line_iss >> zero_atom[all_atom].x >> zero_atom[all_atom].y >> zero_atom[all_atom].z;
          if(zero_atom[all_atom].x < leng_l) leng_l = zero_atom[all_atom].x;
          if(zero_atom[all_atom].x > leng_h) leng_h = zero_atom[all_atom].x;
          if(zero_atom[all_atom].y < leng_l) leng_l = zero_atom[all_atom].y;
          if(zero_atom[all_atom].y > leng_h) leng_h = zero_atom[all_atom].y;
          if(zero_atom[all_atom].z < leng_l) leng_l = zero_atom[all_atom].z;
          if(zero_atom[all_atom].z > leng_h) leng_h = zero_atom[all_atom].z;
          all_atom++;
        }
      }

      if(blank == 3 && temp_line != "")
      {
        int i;
        line_iss >> i;
        i--;
        while(!line_iss.eof())
        {
          line_iss >> zero_atom[i].nei_atom[zero_atom[i].nei_num] >> zero_atom[i].nei_bond[zero_atom[i].nei_num];
          int j = zero_atom[i].nei_atom[zero_atom[i].nei_num]-1;
          double k = zero_atom[i].nei_bond[zero_atom[i].nei_num];
          zero_atom[j].nei_atom[zero_atom[j].nei_num] = i+1;
          zero_atom[j].nei_bond[zero_atom[j].nei_num] = k;
          if(zero_atom[i].max_bond < k) zero_atom[i].max_bond = k;
          if(zero_atom[j].max_bond < k) zero_atom[j].max_bond = k;
          zero_atom[i].nei_num++;
          zero_atom[j].nei_num++;
        }
      }

      if(temp_line == "")
        blank++;
    }
  }
  file_fin.close();

  for(int i=0;i<all_atom;i++)
  {
    if(zero_atom[i].name == "C")
    {
      if(zero_atom[i].max_bond == 2.0) zero_atom[i].name = "CM";
      else if(zero_atom[i].max_bond == 1.5)  zero_atom[i].name = "CA";
      else zero_atom[i].name = "CT";
    }
    else if (zero_atom[i].name == "H")
    {
      int j = zero_atom[i].nei_atom[0]-1;
      if(zero_atom[j].name == "O" || zero_atom[j].name == "OH") zero_atom[i].name = "HO";
      else if(zero_atom[j].name == "CA") zero_atom[i].name = "HA";
      else if(zero_atom[j].name == "C" && zero_atom[j].max_bond == 1.5) zero_atom[i].name = "HA";
      else zero_atom[i].name = "HC";
    }
    else if(zero_atom[i].name == "O")
    {
      int j = zero_atom[i].nei_atom[0]-1, k = zero_atom[i].nei_atom[1]-1;
      if(zero_atom[j].max_bond == 2.0) zero_atom[i].name = "O2";
      else if(zero_atom[j].name == "H" || zero_atom[k].name == "H"
           || zero_atom[j].name == "HO" || zero_atom[k].name == "HO") zero_atom[i].name = "OH";
      else zero_atom[i].name = "OS";
    }
  }

  string atom_type[MAXATOMTYPE];
  int all_atom_type = 0;
  bond zero_bond[max_bond];
  int all_bond = 0;
  angle zero_angle[max_angle];
  int all_angle = 0;
  torsion zero_torsion[max_torsion];
  int all_torsion = 0;
  improper zero_improper[max_improper];
  int all_improper = 0;
  for(int i=0;i<all_atom;i++)
  {
    int find = 0;
    for(int j=0;j<all_atom_type;j++)
      if(atom_type[j] == zero_atom[i].name)
      {
        zero_atom[i].atom_type_id = all_atom_type;
        find++;
      }
    if(find == 0)
    {
      atom_type[all_atom_type] = zero_atom[i].name;
      all_atom_type++;
      zero_atom[i].atom_type_id = all_atom_type;
    }

    for(int j=0;j<zero_atom[i].nei_num;j++)
    {
      int k = zero_atom[i].nei_atom[j];
      if(i+1<k)
      {
        stringstream bond_name;
        bond_name << zero_atom[i].name << "-" << zero_atom[k-1].name;
        zero_bond[all_bond].name = bond_name.str();
        zero_bond[all_bond].atom[0] = i+1;
        zero_bond[all_bond].atom[1] = k;
        all_bond++;
      }
      for(int l=0;l<zero_atom[k-1].nei_num;l++)
      {
        int m = zero_atom[k-1].nei_atom[l];
        if(i+1<k && i+1<m)
        {
          stringstream angle_name;
          angle_name << zero_atom[i].name << "-" << zero_atom[k-1].name << "-" << zero_atom[m-1].name;
          zero_angle[all_angle].name = angle_name.str();
          zero_angle[all_angle].atom[0] = i+1;
          zero_angle[all_angle].atom[1] = k;
          zero_angle[all_angle].atom[2] = m;
          all_angle++;

          for(int n=0;n<zero_atom[i].nei_num;n++)
          {
            int o = zero_atom[i].nei_atom[n];
            if(o<m && o != k)
            {
              stringstream torsion_name;
              torsion_name << zero_atom[o-1].name << "-" << zero_atom[i].name << "-" << zero_atom[k-1].name << "-" << zero_atom[m-1].name;
              zero_torsion[all_torsion].name = torsion_name.str();
              zero_torsion[all_torsion].atom[0] = o;
              zero_torsion[all_torsion].atom[1] = i+1;
              zero_torsion[all_torsion].atom[2] = k;
              zero_torsion[all_torsion].atom[3] = m;
              all_torsion++;
            }
          }

          for(int n=0;n<zero_atom[k-1].nei_num;n++)
          {
            int o = zero_atom[k-1].nei_atom[n];
            if(m<o)
            {
              stringstream improper_name, improper_name2, improper_name3;
              improper_name << zero_atom[i].name << "-" << zero_atom[k-1].name << "-" << zero_atom[m-1].name << "-" << zero_atom[o-1].name;
              zero_improper[all_improper].name = improper_name.str();
              zero_improper[all_improper].atom[0] = i+1;
              zero_improper[all_improper].atom[1] = k;
              zero_improper[all_improper].atom[2] = m;
              zero_improper[all_improper].atom[3] = o;
              all_improper++;
              improper_name2 << zero_atom[m-1].name << "-" << zero_atom[k-1].name << "-" << zero_atom[i].name << "-" << zero_atom[o-1].name;
              zero_improper[all_improper].name = improper_name2.str();
              zero_improper[all_improper].atom[0] = m;
              zero_improper[all_improper].atom[1] = k;
              zero_improper[all_improper].atom[2] = i+1;
              zero_improper[all_improper].atom[3] = o;
              all_improper++;
              improper_name3 << zero_atom[i].name << "-" << zero_atom[k-1].name << "-" << zero_atom[o-1].name << "-" << zero_atom[m-1].name;
              zero_improper[all_improper].name = improper_name3.str();
              zero_improper[all_improper].atom[0] = i+1;
              zero_improper[all_improper].atom[1] = k;
              zero_improper[all_improper].atom[2] = o;
              zero_improper[all_improper].atom[3] = m;
              all_improper++;
            }
          }
        }
      }
      for(int l=0;l<zero_atom[i].nei_num;l++)
      {
        int m = zero_atom[i].nei_atom[l];
        if(i+1<k && k<m)
        {
          stringstream angle_name;
          angle_name << zero_atom[k-1].name << "-" << zero_atom[i].name << "-" << zero_atom[m-1].name;
          zero_angle[all_angle].name = angle_name.str();
          zero_angle[all_angle].atom[0] = k;
          zero_angle[all_angle].atom[1] = i+1;
          zero_angle[all_angle].atom[2] = m;
          all_angle++;
        }
        for(int n=0;n<zero_atom[k-1].nei_num;n++)
        {
          int o = zero_atom[k-1].nei_atom[n];
          if(i+1<k && o<m && o != i+1 && m != k)
          {
            stringstream torsion_name;
            torsion_name << zero_atom[o-1].name << "-" << zero_atom[k-1].name << "-" << zero_atom[i].name << "-" << zero_atom[m-1].name;
            zero_torsion[all_torsion].name = torsion_name.str();
            zero_torsion[all_torsion].atom[0] = o;
            zero_torsion[all_torsion].atom[1] = k;
            zero_torsion[all_torsion].atom[2] = i+1;
            zero_torsion[all_torsion].atom[3] = m;
            all_torsion++;
          }
        }

        for(int n=0;n<zero_atom[i].nei_num;n++)
        {
          int o = zero_atom[i].nei_atom[n];
          if(i+1<k && k<m && m<o)
          {
            stringstream improper_name, improper_name2, improper_name3;
            improper_name << zero_atom[k-1].name << "-" << zero_atom[i].name << "-" << zero_atom[m-1].name << "-" << zero_atom[o-1].name;
            zero_improper[all_improper].name = improper_name.str();
            zero_improper[all_improper].atom[0] = k;
            zero_improper[all_improper].atom[1] = i+1;
            zero_improper[all_improper].atom[2] = m;
            zero_improper[all_improper].atom[3] = o;
            all_improper++;
            improper_name2 << zero_atom[m-1].name << "-" << zero_atom[i].name << "-" << zero_atom[k-1].name << "-" << zero_atom[o-1].name;
            zero_improper[all_improper].name = improper_name2.str();
            zero_improper[all_improper].atom[0] = m;
            zero_improper[all_improper].atom[1] = i+1;
            zero_improper[all_improper].atom[2] = k;
            zero_improper[all_improper].atom[3] = o;
            all_improper++;
            improper_name3 << zero_atom[k-1].name << "-" << zero_atom[i].name << "-" << zero_atom[o-1].name << "-" << zero_atom[m-1].name;
            zero_improper[all_improper].name = improper_name3.str();
            zero_improper[all_improper].atom[0] = k;
            zero_improper[all_improper].atom[1] = i+1;
            zero_improper[all_improper].atom[2] = o;
            zero_improper[all_improper].atom[3] = m;
            all_improper++;
          }
        }

      }
    }
  }

  string bond_type[MAXBONDTYPE];
  int all_bond_type = 0;
  for(int i=0;i<all_bond;i++)
  {
    int find = 0;
    stringstream bond_name_inv;
    bond_name_inv << zero_atom[zero_bond[i].atom[1]-1].name << "-" << zero_atom[zero_bond[i].atom[0]-1].name;
    for(int j=0;j<all_bond_type;j++)
      if(bond_type[j] == zero_bond[i].name  || bond_type[j] == bond_name_inv.str())
      {
        zero_bond[i].bond_type_id = j+1;
        find++;
      }
    if(find == 0)
    {
      bond_type[all_bond_type] = zero_bond[i].name;
      all_bond_type++;
      zero_bond[i].bond_type_id = all_bond_type;
    }
  }

  string angle_type[MAXANGLTYPE];
  int all_angle_type = 0;
  for(int i=0;i<all_angle;i++)
  {
    int find = 0;
    stringstream angle_name_inv;
    angle_name_inv << zero_atom[zero_angle[i].atom[2]-1].name << "-" << zero_atom[zero_angle[i].atom[1]-1].name << "-"
                   << zero_atom[zero_angle[i].atom[0]-1].name;
    for(int j=0;j<all_angle_type;j++)
      if(angle_type[j] == zero_angle[i].name  || angle_type[j] == angle_name_inv.str())
      {
        zero_angle[i].angle_type_id = j+1;
        find++;
      }
    if(find == 0)
    {
      angle_type[all_angle_type] = zero_angle[i].name;
      all_angle_type++;
      zero_angle[i].angle_type_id = all_angle_type;
    }
  }

/*
for(int i=0;i<all_torsion;i++) cerr << i+1 << " " << zero_torsion[i].name //<< " " << zero_torsion[i].torsion_type_id 
                                    << " " << zero_torsion[i].atom[0] << " " << zero_torsion[i].atom[1]
                                    << " " << zero_torsion[i].atom[2] << " " << zero_torsion[i].atom[3] << endl;
cerr << endl;
for(int i=0;i<all_improper;i++) cerr << i+1 << " " << zero_improper[i].name //<< " " << zero_improper[i].improper_type_id 
                                     << " " << zero_improper[i].atom[0] << " " << zero_improper[i].atom[1]
                                     << " " << zero_improper[i].atom[2] << " " << zero_improper[i].atom[3] << endl;
cerr << endl;
*/

  for(int i=0;i<all_torsion;i++)
  {
    stringstream torsion_name;
    if(zero_atom[zero_torsion[i].atom[0]-1].name == "OS" && zero_atom[zero_torsion[i].atom[3]-1].name == "OS")
      torsion_name << zero_torsion[i].name;
    else if(zero_atom[zero_torsion[i].atom[0]-1].name == "OS" && zero_atom[zero_torsion[i].atom[3]-1].name == "OH")
      torsion_name << zero_torsion[i].name;
    else if(zero_atom[zero_torsion[i].atom[0]-1].name == "OH" && zero_atom[zero_torsion[i].atom[3]-1].name == "OS")
      torsion_name << zero_torsion[i].name;
    else if(zero_atom[zero_torsion[i].atom[0]-1].name == "OH" && zero_atom[zero_torsion[i].atom[3]-1].name == "OH")
      torsion_name << zero_torsion[i].name;

    else if(zero_atom[zero_torsion[i].atom[1]-1].name == "CT" && zero_atom[zero_torsion[i].atom[2]-1].name == "CT")
      torsion_name << "XX-CT-CT-XX";
    else if(zero_atom[zero_torsion[i].atom[1]-1].name == "CA" && zero_atom[zero_torsion[i].atom[2]-1].name == "CA")
      torsion_name << "XX-CA-CA-XX";
    else if(zero_atom[zero_torsion[i].atom[1]-1].name == "CM" && zero_atom[zero_torsion[i].atom[2]-1].name == "CM")
      torsion_name << "XX-CM-CM-XX";

    else if(zero_atom[zero_torsion[i].atom[1]-1].name == "CT" && zero_atom[zero_torsion[i].atom[2]-1].name == "CA")
      torsion_name << "XX-CT-CA-XX";
    else if(zero_atom[zero_torsion[i].atom[1]-1].name == "CA" && zero_atom[zero_torsion[i].atom[2]-1].name == "CT")
      torsion_name << "XX-CA-CT-XX";
    else if(zero_atom[zero_torsion[i].atom[1]-1].name == "CT" && zero_atom[zero_torsion[i].atom[2]-1].name == "CM")
      torsion_name << "XX-CT-CM-XX";
    else if(zero_atom[zero_torsion[i].atom[1]-1].name == "CM" && zero_atom[zero_torsion[i].atom[2]-1].name == "CT")
      torsion_name << "XX-CM-CT-XX";
    else if(zero_atom[zero_torsion[i].atom[1]-1].name == "CA" && zero_atom[zero_torsion[i].atom[2]-1].name == "CM")
      torsion_name << "XX-CA-CM-XX";
    else if(zero_atom[zero_torsion[i].atom[1]-1].name == "CM" && zero_atom[zero_torsion[i].atom[2]-1].name == "CA")
      torsion_name << "XX-CM-CA-XX";

    else torsion_name << zero_torsion[i].name;
    zero_torsion[i].name = torsion_name.str();
  }

  for(int i=0;i<all_improper;i++)
  {
    stringstream improper_name;
    if(zero_atom[zero_improper[i].atom[2]-1].name == "CA" && zero_atom[zero_improper[i].atom[3]-1].name == "HA")
      improper_name << "XX-XX-CA-HA";
    else if(zero_atom[zero_improper[i].atom[0]-1].name == "HA" && zero_atom[zero_improper[i].atom[1]-1].name == "CA")
      improper_name << "HA-CA-XX-XX";
    else if(zero_atom[zero_improper[i].atom[2]-1].name == "CM" && zero_atom[zero_improper[i].atom[3]-1].name == "HA")
      improper_name << "XX-XX-CM-HA";
    else if(zero_atom[zero_improper[i].atom[0]-1].name == "HA" && zero_atom[zero_improper[i].atom[1]-1].name == "CM")
      improper_name << "HA-CM-XX-XX";

    else improper_name << "NONE";
    zero_improper[i].name = improper_name.str();
  }

  string torsion_type[MAXTORSTYPE];
  int all_torsion_type = 0;
  for(int i=0;i<all_torsion;i++)
  {
    int find = 0;
    stringstream torsion_name_inv;
    torsion_name_inv << zero_atom[zero_torsion[i].atom[3]-1].name << "-" << zero_atom[zero_torsion[i].atom[2]-1].name << "-"
                     << zero_atom[zero_torsion[i].atom[1]-1].name << "-" << zero_atom[zero_torsion[i].atom[0]-1].name;
    for(int j=0;j<all_torsion_type;j++)
      if(torsion_type[j] == zero_torsion[i].name  || torsion_type[j] == torsion_name_inv.str())
      {
        zero_torsion[i].torsion_type_id = j+1;
        find++;
      }
    if(find == 0)
    {
      torsion_type[all_torsion_type] = zero_torsion[i].name;
      all_torsion_type++;
      zero_torsion[i].torsion_type_id = all_torsion_type;
    }
  }

  string improper_type[MAXTORSTYPE];
  int all_improper_type = 0;
  for(int i=0;i<all_improper;i++)
  {
    int find = 0;
    stringstream improper_name_inv;
    improper_name_inv << zero_atom[zero_improper[i].atom[3]-1].name << "-" << zero_atom[zero_improper[i].atom[2]-1].name << "-"
                      << zero_atom[zero_improper[i].atom[1]-1].name << "-" << zero_atom[zero_improper[i].atom[0]-1].name;
    for(int j=0;j<all_improper_type;j++)
      if(improper_type[j] == zero_improper[i].name  || improper_type[j] == improper_name_inv.str())
      {
        zero_improper[i].improper_type_id = j+1;
        find++;
      }
    if(find == 0)
    {
      improper_type[all_improper_type] = zero_improper[i].name;
      all_improper_type++;
      zero_improper[i].improper_type_id = all_improper_type;
    }
  }

  if(opts.atomratio == "lmp"){
    cout << title << endl << endl;
    cout << " " << all_atom << " atoms" << endl;
    cout << " " << all_bond << " bonds" << endl;
    cout << " " << all_angle << " angles" << endl;
    cout << " " << all_torsion+all_improper << " dihedrals" << endl;
    cout << " 0 impropers" << endl; 
    cout << " " << all_atom_type << " atom types" << endl;
    cout << " " << all_bond_type << " bond types" << endl;
    cout << " " << all_angle_type << " angle types" << endl;
    cout << " " << all_torsion_type+all_improper_type << " dihedral types" << endl;
    cout << " 0 improper types" << endl; 
    cout << " " << leng_l*3 << " " << leng_h*3 << " xlo xhi" << endl;
    cout << " " << leng_l*3 << " " << leng_h*3 << " ylo yhi" << endl;
    cout << " " << leng_l*3 << " " << leng_h*3 << " zlo zhi" << endl;
    cout << endl << " Masses" << endl<< endl;
    for(int i=0;i<all_atom_type;i++) cout << i+1 << " " << "mass" << " # " << atom_type[i] << endl;
    cout << endl << " PairIJ Coeffs # bond_style lennard/mdf, pair_coeff A B r_m r_c" << endl<< endl;
    for(int i=0;i<all_atom_type;i++)
      for(int j=i;j<all_atom_type;j++) 
        cout << i+1 << " " << j+1 << " " << "A" << " " << "B" << " 8.0 9.0 # " << atom_type[i] << "-" << atom_type[j] << endl;
    cout << endl << " Bond Coeffs # bond_style harmonic, bond_coeff K r" << endl<< endl;
    for(int i=0;i<all_bond_type;i++) cout << i+1 << " " << "K" << " " << "r" << " # " << bond_type[i] << endl;
    cout << endl << " Angle Coeffs # angle_style harmonic, angle_coeff K th" << endl<< endl;
    for(int i=0;i<all_angle_type;i++) cout << i+1 << " " << "K" << " " << "th" << " # " << angle_type[i] << endl;
    cout << endl << " Dihedral Coeffs # dihedral_style charmm, dihedral_coeff (V/2)*paths n gamma 0.0" << endl<< endl;
    for(int i=0;i<all_torsion_type;i++) cout << i+1
                                             << " " << "(V/2)*paths" << " " << "n" << " " << "gamma" << " 0.0 # " << torsion_type[i] << " torsion" << endl;
    for(int i=0;i<all_improper_type;i++) cout << i+1+all_torsion_type
                                              << " " << "(V/2)*paths" << " " << "n" << " " << "gamma" << " 0.0 # " << improper_type[i] << " improper" << endl;
    cout << endl << " Atoms # atom_style full, atm-id mol-id atm-type q x y z" << endl << endl;
    for(int i=0;i<all_atom;i++)
      cout << i+1 << " " << "255" << " " << zero_atom[i].atom_type_id << " " << "0.0" << " "
                         << zero_atom[i].x << " " << zero_atom[i].y << " " << zero_atom[i].z << " # " << zero_atom[i].name << endl;
    cout << endl << " Bonds" << endl << endl;
    for(int i=0;i<all_bond;i++)
      cout << i+1 << " " << zero_bond[i].bond_type_id << " "
                         << zero_bond[i].atom[0] << " " << zero_bond[i].atom[1] << " # " << zero_bond[i].name << endl;
    cout << endl << " Angles" << endl << endl;
    for(int i=0;i<all_angle;i++)
      cout << i+1 << " " << zero_angle[i].angle_type_id << " "
                         << zero_angle[i].atom[0] << " " << zero_angle[i].atom[1] << " " << zero_angle[i].atom[2] << " # " << zero_angle[i].name << endl;
    cout << endl << " Dihedrals" << endl << endl;
    for(int i=0;i<all_torsion;i++)
      cout << i+1 << " " << zero_torsion[i].torsion_type_id << " "
                         << zero_torsion[i].atom[0] << " " << zero_torsion[i].atom[1] << " "
                         << zero_torsion[i].atom[2] << " " << zero_torsion[i].atom[3] << " # " << zero_torsion[i].name << " torsion" << endl;
    for(int i=0;i<all_improper;i++)
      cout << all_torsion+i+1 << " " << zero_improper[i].improper_type_id+all_torsion_type << " "
                         << zero_improper[i].atom[0] << " " << zero_improper[i].atom[1] << " "
                         << zero_improper[i].atom[2] << " " << zero_improper[i].atom[3] << " # " << zero_improper[i].name << " improper" << endl;
  }

}

