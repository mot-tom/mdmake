#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <math.h>
#include "opt.h"
using namespace std;

// Outline of this program 
/*
  This program is for oxidize "NW" whose z axis is long axis.
  You have to make "layer(number)" in the same directory.
*/

double round(double round, double digit);   // function to round off the number

struct info_oxi        // information of atom
{
  string name;     // atom name
  string state;    // state of (O) atom; "good":can add; "bad":can't add
  string status;   // status of (O) atoms; under("yet") / on ("new") oxide layer currently
  double x, y, z;  // coordinate of atom
  int layer;       // layer of atom
  int check;       // whether to add (O) atoms in the current loop; "1":add; "0":don't add
};

int oxi_layer(int argc, char **argv, optpara opts){
        int nx = 1, ny = 1, nz = 1;
        if(opts.atomcomp == "zl-SiO2" || opts.atomcomp == "zl-GeO2") {nx = 4; ny = 4; nz = 2;}
        else {nx = 4; ny = 4; nz = 4;}

        if(opts.cellnum.length() > 0) {
          istringstream cell_iss(opts.cellnum.c_str());
          cell_iss >> nx >> ny >> nz;
          if(nx == 0 || ny == 0 || nz == 0) {
            cerr << "Cell number is wrong" << endl;
            return 1;
           }}

  /******************************/
  
  const int oxi_times = 10;       // number of times of adding O atoms  
  const int oxi_check = 5;        // check for oxidation finished
  int max_oxi_lay = 8;      // maximum layers to oxidize
  const int max_add = 10000;     // number of allocating dynamic memory
  
  string input = "stable.txt";    // name of the optimized input file
  string output = "add_O.txt";    // name of the output file
  string log = "log.txt";         // log of oxidation process
  
  /******************************/

  double lattice;       // lattice constant of Si / Ge [A]
  double dx, dy, dz;    // squared distance of each direction between atom i and atom j [A^2]  
  double cut_Si_Ge;     // squared cutoff distance of Si-Si / Ge-Ge [A^2]
  double cut_O;         // squared cutoff distance of O-O [A^2]
  double r, rx, ry, rz; // squared distance of two atoms [A^2] (rz is considered boundary condition)
  double rxy;
  
  int oxi_lay;          // oxide layer currently
  int times;            // time of oxidation process currently
  int dir_num=2;
  
  string atom_name;     // atom name (Si / Ge)
  string layer_dir;     // layer direction
  string md_st[3];      // line stored structure in the MDlabo
  string empty = "";    // initialize stringstream

  int debug = 0;

  ostringstream os;

  if(opts.atomcomp == "zl-SiO2")
  {
    atom_name = "Si";
    layer_dir = "zl";
    max_oxi_lay = nz-1;
    os << "mdmake -c zl-Si -n \"" << nx << " " << ny << " " << nz << "\" >str.txt 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);

    if(debug) os << "mdmake --convert mdlFILM -c str.txt >" << input;
    else os << "mdmake --convert mdlFILM -c str.txt >" << input << " 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);      
  }
  else if(opts.atomcomp == "zl1-SiO2")
  {
    atom_name = "Si";
    layer_dir = "zl";
    max_oxi_lay = nz-1;
    os << "mdmake -c zl1-Si -n \"" << nx << " " << ny << " " << nz << "\" >str.txt 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);

    if(debug) os << "mdmake --convert mdlFILM -c str.txt >" << input;
    else os << "mdmake --convert mdlFILM -c str.txt >" << input << " 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);      
  }
  else if(opts.atomcomp == "zl-3T-SiO2")
  {
    atom_name = "Si";
    layer_dir = "zl";
    max_oxi_lay = nz-1;
    os << "mdmake -c zl-3T-Si -n \"" << nx << " " << ny << " " << nz << "\" >str.txt 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);

    if(debug) os << "mdmake --convert mdlFILM -c str.txt >" << input;
    else os << "mdmake --convert mdlFILM -c str.txt >" << input << " 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);      
  }
  else if(opts.atomcomp == "zl1-3T-SiO2")
  {
    atom_name = "Si";
    layer_dir = "zl";
    max_oxi_lay = nz-1;
    os << "mdmake -c zl1-3T-Si -n \"" << nx << " " << ny << " " << nz << "\" >str.txt 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);

    if(debug) os << "mdmake --convert mdlFILM -c str.txt >" << input;
    else os << "mdmake --convert mdlFILM -c str.txt >" << input << " 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);      
  }
  else if(opts.atomcomp == "zl-3H-SiO2")
  {
    atom_name = "Si";
    layer_dir = "zl";
    max_oxi_lay = nz-1;
    os << "mdmake -c zl-3H-Si -n \"" << nx << " " << ny << " " << nz << "\" >str.txt 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);

    if(debug) os << "mdmake --convert mdlFILM -c str.txt >" << input;
    else os << "mdmake --convert mdlFILM -c str.txt >" << input << " 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);      
  }
  else if(opts.atomcomp == "zl1-3H-SiO2")
  {
    atom_name = "Si";
    layer_dir = "zl";
    max_oxi_lay = nz-1;
    os << "mdmake -c zl1-3H-Si -n \"" << nx << " " << ny << " " << nz << "\" >str.txt 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);

    if(debug) os << "mdmake --convert mdlFILM -c str.txt >" << input;
    else os << "mdmake --convert mdlFILM -c str.txt >" << input << " 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);      
  }
  else if(opts.atomcomp == "zl-3HO-SiO2")
  {
    atom_name = "Si";
    layer_dir = "zl";
    max_oxi_lay = nz-1;
    os << "mdmake -c zl-3HO-Si -n \"" << nx << " " << ny << " " << nz << "\" >str.txt 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);

    if(debug) os << "mdmake --convert mdlFILM -c str.txt >" << input;
    else os << "mdmake --convert mdlFILM -c str.txt >" << input << " 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);      
  }
  else if(opts.atomcomp == "zl1-3HO-SiO2")
  {
    atom_name = "Si";
    layer_dir = "zl";
    max_oxi_lay = nz-1;
    os << "mdmake -c zl1-3HO-Si -n \"" << nx << " " << ny << " " << nz << "\" >str.txt 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);

    if(debug) os << "mdmake --convert mdlFILM -c str.txt >" << input;
    else os << "mdmake --convert mdlFILM -c str.txt >" << input << " 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);      
  }
  else if(opts.atomcomp == "zl-3TO-SiO2")
  {
    atom_name = "Si";
    layer_dir = "zl";
    max_oxi_lay = nz-1;
    os << "mdmake -c zl-3TO-Si -n \"" << nx << " " << ny << " " << nz << "\" >str.txt 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);

    if(debug) os << "mdmake --convert mdlFILM -c str.txt >" << input;
    else os << "mdmake --convert mdlFILM -c str.txt >" << input << " 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);      
  }
  else if(opts.atomcomp == "zl1-3TO-SiO2")
  {
    atom_name = "Si";
    layer_dir = "zl";
    max_oxi_lay = nz-1;
    os << "mdmake -c zl1-3TO-Si -n \"" << nx << " " << ny << " " << nz << "\" >str.txt 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);

    if(debug) os << "mdmake --convert mdlFILM -c str.txt >" << input;
    else os << "mdmake --convert mdlFILM -c str.txt >" << input << " 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);      
  }
  else if(opts.atomcomp == "zl-GeO2")
  {
    atom_name = "Ge";
    layer_dir = "zl";
    max_oxi_lay = nz-1;
    os << "mdmake -c zl-Ge -l -n \"" << nx << " " << ny << " " << nz << "\" >str.txt 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);

    os << "mdmake --convert mdlFILM -c str.txt -l >" << input << " 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);      
  }
  else
  {
    os << "mdmake --convert mdl -c " << opts.atomcomp << " -n \"" << nx << " " << ny << " " << nz << "\" >" << input << " 2>/dev/null";
    cout << os.str() << endl;
    system(os.str().c_str());
    os.str(empty);
  }

      if(atom_name!="Si" && atom_name!="Ge" && atom_name!="SiO" && atom_name!="GeO")
      {
        cout << "Setting of atom name is wrong!" << endl;
        return 1;
      }

      // determin lattice constant and cut off distance
      if(atom_name=="Si" && atom_name!="SiO")
      {
        lattice = 5.43;
        cut_Si_Ge = pow(2.8, 2.0);
        cut_O = pow(0.8, 2.0);
      }
      if(atom_name=="Ge" && atom_name!="GeO")
      {
        lattice = 5.65;
        cut_Si_Ge = pow(3.0, 2.0);
        cut_O = pow(1.0, 2.0);
      }

  if(atom_name == "Si")
  {
    if(debug) os << "mdmake --convert lmpLAY -c " << input << " >stable.lmp";
    else os << "mdmake --convert lmpLAY -c " << input << " >stable.lmp 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);

    os << "in.stable.Si";
    ofstream fout_in_stable(os.str().c_str());
    os.str(empty);

        fout_in_stable << "package      omp " << opts.omp << endl;
        fout_in_stable << "units        metal" << endl;
        fout_in_stable << "boundary     p p p" << endl;
        fout_in_stable << "atom_style   molecular" << endl;
        fout_in_stable << "atom_modify  sort 10000 1.0" << endl;
        fout_in_stable << "read_data    stable.lmp" << endl;
        fout_in_stable << "pair_style   sw/wfnho" << endl;
        fout_in_stable << "pair_coeff   * * SiO.sw_wfnho Si" << endl;
        fout_in_stable << "neighbor     4.0 bin" << endl;
        fout_in_stable << "neigh_modify every 1 delay 0 check yes" << endl;
        fout_in_stable << "group        fix_lay molecule 255" << endl;
        fout_in_stable << "fix          1 all nve" << endl;
        fout_in_stable << "fix          2 all box/relax aniso 0.0 fixedpoint 0.0 0.0 0.0" << endl;
        fout_in_stable << "fix          3 fix_lay setforce 0 0 0 " << endl;
        fout_in_stable << "min_style    cg" << endl;
        fout_in_stable << "minimize     1e-25 1e-25 50000 100000" << endl;
        fout_in_stable << "dump         1 all custom 1 stable.final id type xs ys zs mol" << endl;
        fout_in_stable << "dump_modify  1 sort id" << endl;
        fout_in_stable << "run          0" << endl;
        fout_in_stable << "undump       1" << endl << endl;
        fout_in_stable.close();

    os << "mpirun -np " << opts.omp << " lmp_custom <in.stable.Si >/dev/null 2>&1";
    system(os.str().c_str());
    os.str(empty);

    os << "mdmake --convert mdl --lammps Si -c stable.final >" << input << " 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);
    atom_name = "SiO";
  }
  else if(atom_name == "Ge")
  {
    os << "cat " << input << " >str.txt";
    system(os.str().c_str());
    os.str(empty);

    os << "head.stable";
    ofstream fout_in_stable(os.str().c_str());
    os.str(empty);

    fout_in_stable << "OPT(P) RAND=10 ESWGE T=300.0 STEP=(0,20000)" << endl;
    fout_in_stable << "BLOCK W=INF FIXANGLE Q=1.0" << endl;
    fout_in_stable.close();

    os << "tail -n +3 str.txt | cat head.stable - >" << input;
    system(os.str().c_str());
    os.str(empty);

    os << "mdlabo < " << input << "> stable.final";
    system(os.str().c_str());
    os.str(empty);

    os << "mdmake --convert mdl -c stable.final -l >" << input << " 2>/dev/null";
    system(os.str().c_str());
    os.str(empty);

    atom_name = "GeO";
  }
  if(atom_name == "SiO")
  {
    os << "in.stable.SiO";
    ofstream fout_in_stable(os.str().c_str());
    os.str(empty);

        fout_in_stable << "package      omp " << opts.omp << endl;
        fout_in_stable << "units        metal" << endl;
        fout_in_stable << "boundary     p p p" << endl;
        fout_in_stable << "atom_style   molecular" << endl;
        fout_in_stable << "atom_modify  sort 10000 1.0" << endl;
        fout_in_stable << "read_data    stable.lmp" << endl;
        fout_in_stable << "pair_style   sw/wfnho" << endl;
        fout_in_stable << "pair_coeff   * * SiO.sw_wfnho Si O" << endl;
        fout_in_stable << "neighbor     4.0 bin" << endl;
        fout_in_stable << "neigh_modify every 1 delay 0 check yes" << endl;
        fout_in_stable << "group        fix_lay molecule 255" << endl;
        fout_in_stable << "fix          1 all nve" << endl;
        fout_in_stable << "fix          2 all box/relax aniso 0.0 fixedpoint 0.0 0.0 0.0" << endl;
        fout_in_stable << "fix          3 fix_lay setforce 0 0 0" << endl;
        fout_in_stable << "min_style    cg" << endl;
        fout_in_stable << "minimize     1e-25 1e-25 50000 100000" << endl;
        fout_in_stable << "dump         1 all custom 1 stable.final id type xs ys zs mol" << endl;
        fout_in_stable << "dump_modify  1 sort id" << endl;
        fout_in_stable << "run          0" << endl;
        fout_in_stable << "undump       1" << endl << endl;
        fout_in_stable.close();

    os << "in.md.SiO";
    ofstream fout_in_md(os.str().c_str());
    os.str(empty);

        fout_in_md << "package      omp " << opts.omp << endl;
        fout_in_md << "units        metal" << endl;
        fout_in_md << "boundary     p p p" << endl;
        fout_in_md << "atom_style   molecular" << endl;
        fout_in_md << "atom_modify  sort 10000 1.0" << endl;
        fout_in_md << "read_data    md.lmp" << endl;
        fout_in_md << "pair_style   sw/wfnho" << endl;
        fout_in_md << "pair_coeff   * * SiO.sw_wfnho Si O" << endl;
        fout_in_md << "neighbor     4.0 bin" << endl;
        fout_in_md << "neigh_modify every 1 delay 0 check yes" << endl;
        fout_in_md << "group        fix_lay molecule 255" << endl;
        fout_in_md << "group        mov_lay molecule != 255" << endl;
        fout_in_md << "fix          1 all npt temp 1000 1000 0.05 x 0.0 0.0 0.5 y 0.0 0.0 0.5 z 0.0 0.0 0.5" << endl;
        fout_in_md << "fix          2 fix_lay setforce 0 0 0" << endl;
        fout_in_md << "velocity     mov_lay create 1000 458127641 mom yes rot yes dist gaussian" << endl;
        fout_in_md << "timestep     0.0001" << endl;
        fout_in_md << "dump         1 all xyz 100 oxi.xyz" << endl; //takizawa
        fout_in_md << "dump_modify  1 sort id element Si O append no" << endl; //takizawa
        fout_in_md << "run          20000" << endl;
        fout_in_md << "undump       1" << endl; //takizawa
        fout_in_md << "dump         2 all custom 1 md.final id type xs ys zs mol" << endl;
        fout_in_md << "dump_modify  2 sort id" << endl;
        fout_in_md << "run          0" << endl;
        fout_in_md << "undump       2" << endl << endl;
        fout_in_md.close();
  }
  else if(atom_name == "GeO")
  {
    os << "head.md";
    ofstream fout_in_md(os.str().c_str());
    os.str(empty);

    fout_in_md << "Gear(5) PRECISE(12) RAND=10 ESW T=1000.0 STEP=(0,20000)" << endl;
    fout_in_md << "BLOCK W=INF FIXANGLE Q=1.0 dt=0.1fs" << endl;
    fout_in_md.close();
  }
  
  system("mkdir layer2");
  os << "cp -P " << input << " layer2/stable0_" << oxi_times << "-" << oxi_check << ".txt";
  system(os.str().c_str());
  os.str(empty);

  cout << "Oxidation Start!" << endl;

  for(oxi_lay=2; oxi_lay<=max_oxi_lay; oxi_lay++)
  {   
    ofstream fout_log(log.c_str());
    if(!fout_log)
    {
      cout << "Can't open " << log << "!" << endl;
      return 1;
    }

    if(oxi_lay > dir_num)
    {
      os << "mkdir layer" << oxi_lay;
      system(os.str().c_str());
      os.str(empty);
    }

    for(int times=1; times<=oxi_times+oxi_check; times++)
    {
      string temp;                // line of string in the file
      string first;               // first section of each line

      int count = 0;              // number of atoms read from the file currently
      int blank = 0;              // number of blanks in the file currently

      int Si_count = 0;           // number of Si atoms
      int Ge_count = 0;           // number of Ge atoms
      int Si_Ge_count = 0;        // number of Si+Ge atoms
      int O_count = 0;            // number of O atoms
      int all_count = 0;          // number of all atoms
      int add_count = 0;          // number of atoms needed for oxidation of the layer
      int yet_count = 0;          // number of "yet" state atoms
      int new_count = 0;          // number of "new" state atoms

      int add_num = 0;            // number of adding atoms in current oxidation process
      int add_yet_num = 0;        // number of adding "yet" state atoms in current oxidation process   
      int add_new_num = 0;        // number of adding "new" state atoms in current oxidation process

      double lx, ly, lz;          // length of each side of box in the MDlabo
      double ratio;               // ratio of adding O atoms

      // allocate dynamic memory
      info_oxi* add_atom = new info_oxi[max_add];   // pointer for added O atom
 
      fout_log << "****************************************" << endl;
      fout_log << "Layer:" << oxi_lay << "  ";    
      cout << "Layer:" << oxi_lay << "  ";
      if(times<=oxi_times)
      {
        fout_log << "Times:" << times << "/" << oxi_times << endl << endl;
        cout << "Times:" << times << "/" << oxi_times << "  " << flush;
      }
      else
      {
        fout_log << "Check:" << times-oxi_times << "/" << oxi_check << endl << endl;
        cout << "Check:" << times-oxi_times << "/" << oxi_check << "  " << flush;
      }

      // read atom name / condition of MDlabo / number of each atom
      ifstream fin(input.c_str());
      if(!fin)
      {
        cout << "Can't open " << input << '!' << endl;
        return 1;
      }
      
      while(getline(fin, temp))
      {
        istringstream is(temp);
        info_oxi atom;                // information of atom temporarily

        if(blank==2 && temp!="")
        {
          is >> atom.x >> atom.y >> atom.z;

          if(atom.x!=0.0)
          {
            lx = atom.x;
            md_st[0] = temp;
          }
          if(atom.y!=0.0)
          {
            ly = atom.y;
            md_st[1] = temp;
          }
          if(atom.z!=0.0)
          {      
            lz = atom.z;
            md_st[2] = temp;
          }
        }
        
        if(blank==3 && temp!="")
        {
          is >> first;
                   
          if(first=="Si")
            Si_count++;
          if(first=="Ge")
            Ge_count++;          
          if(first.find("O")!=string::npos)
            O_count++;
        }
        
        if(temp=="")
          blank++;
      }
      fin.close();

      Si_Ge_count = Si_count + Ge_count;
      all_count = Si_Ge_count + O_count;

      fout_log << "all:" << all_count << "  Si:" << Si_count << "  Ge:" << Ge_count << "  " << "O:" << O_count << endl;
      fout_log << "x:" << round(lx/10.0, 2.0) << "[nm]" << "  " << "y:" << round(ly/10.0, 2.0) << "[nm]"<< "  "
               << "z:" << round(lz/10.0, 2.0) << "[nm]" << endl << endl;

      info_oxi atom[all_count];            // information of all atoms

      ofstream fout(output.c_str());
      if(!fout)
      {
        cout << "Can't open " << output << '!' << endl;
        return 1;
      }

      if(atom_name=="Si" || atom_name=="SiO")
        fout << "OPT(P) RAND=10 ESW T=300.0 STEP=(0,20000)" << endl;
      if(atom_name=="Ge" || atom_name=="GeO")
        fout << "OPT(P) RAND=10 ESWGE T=300.0 STEP=(0,20000)" << endl;
      fout << "BLOCK W=INF FIXANGLE Q=1.0" << endl << endl;
      fout << "oxidation" << endl << endl;

      fout << md_st[0] << endl;
      fout << md_st[1] << endl;
      fout << md_st[2] << endl << endl;   

      // read coordinate of each atoms
      blank = 0;      
      ifstream fin2(input.c_str());
      if(!fin2)
      {
        cout << "Can't open " << input << '!' << endl;
        return 1;
      }

      while(getline(fin2, temp))
      {
        istringstream is(temp);
        if(blank==3 && temp!="")
        {
          fout << temp << endl;
          is >> first;
          if(first=="Si" || first=="Ge")
          { 
            string ch[6];    // "F":atom fixed / "M":atom marked / "L":atom label / "G":atom layer / number:value of "L" or "G"
            atom[count].name = first;
            is >> atom[count].x >> atom[count].y >> atom[count].z
               >> ch[0] >> ch[1] >> ch[2] >> ch[3] >> ch[4] >> ch[5];
            for(int i=0; i<5; i++)
            {
              if(ch[i].find("G")!=string::npos)
              {
                istringstream is_lay(ch[i+1].c_str());
                is_lay >> atom[count].layer;
              }
              else if(ch[i].find("F")!=string::npos)
              {
                atom[count].layer = 255;
              }
            }
            count++;
          }
          if(first.find("O")!=string::npos)
          { 
            atom[count].name = first;
            is >> atom[count].x >> atom[count].y >> atom[count].z;
            count++;
          }
        }
        if(temp=="")
          blank++;
      }
      fin2.close();

      // add O atoms to midpoint of Si-Si
      for(int i=0; i<Si_Ge_count-1; i++)
      {
        for(int j=i+1; j<Si_Ge_count; j++)
        {
          int scene = 0;     // value for deside "status"; 1:"yet"; 2:"new"
          if(atom[i].layer<oxi_lay && atom[j].layer<oxi_lay && abs(atom[i].layer-atom[j].layer)<=1)
            scene = 1;
          if((atom[i].layer==oxi_lay && atom[j].layer==oxi_lay-1) || (atom[i].layer==oxi_lay && atom[j].layer==oxi_lay) || 
             (atom[i].layer==oxi_lay-1 && atom[j].layer==oxi_lay))
            scene = 2;
          
          if(scene==1 || scene==2)
          {	    
            if(max_add==add_count)
              cout << "Dynamic memory space for O atom is not enough!" << endl;
            
            dx = fabs(atom[i].x - atom[j].x) * lx;
            dy = fabs(atom[i].y - atom[j].y) * ly;
            dz = fabs(atom[i].z - atom[j].z) * lz;
            
            r   = pow(dx, 2.0) + pow(dy, 2.0) + pow(dz, 2.0);
            rx  = pow(lx-dx, 2.0) + pow(dy, 2.0) + pow(dz, 2.0);
            ry  = pow(dx, 2.0) + pow(ly-dy, 2.0) + pow(dz, 2.0);
            rz  = pow(dx, 2.0) + pow(dy, 2.0) + pow(lz-dz, 2.0);
            rxy  = pow(lx-dx, 2.0) + pow(ly-dy, 2.0) + pow(dz, 2.0);
            
            add_atom[add_count].name = "O";
            add_atom[add_count].check = 0;
            
            if(scene==1)
              add_atom[add_count].status = "yet";
            if(scene==2)
              add_atom[add_count].status = "new";

            if(dx<=lx/2.0 && dy<=ly/2.0 && r<cut_Si_Ge && r!=0 && layer_dir == "zl")
            {
              add_atom[add_count].x = (atom[i].x + atom[j].x) / 2.0;
              add_atom[add_count].y = (atom[i].y + atom[j].y) / 2.0;
              add_atom[add_count].z = (atom[i].z + atom[j].z) / 2.0;
              add_atom[add_count].state = "good";
              add_count++;
            }
            if(dx>lx/2.0 && dy<=ly/2.0 && rx<cut_Si_Ge && layer_dir == "zl")
            {
              if(atom[i].x + atom[j].x > 1.0)
                add_atom[add_count].x = (atom[i].x + atom[j].x - 1.0) / 2.0;
              else
                add_atom[add_count].x = (atom[i].x + atom[j].x + 1.0) / 2.0;

              add_atom[add_count].z = (atom[i].z + atom[j].z) / 2.0;
              add_atom[add_count].y = (atom[i].y + atom[j].y) / 2.0;
              add_atom[add_count].state = "good";
              add_count++;
            }
            if(dx<=lx/2.0 && dy>ly/2.0 && ry<cut_Si_Ge && layer_dir == "zl")
            {
              if(atom[i].y + atom[j].y > 1.0)
                add_atom[add_count].y = (atom[i].y + atom[j].y - 1.0) / 2.0;
              else
                add_atom[add_count].y = (atom[i].y + atom[j].y + 1.0) / 2.0;
              
              add_atom[add_count].x = (atom[i].x + atom[j].x) / 2.0;
              add_atom[add_count].z = (atom[i].z + atom[j].z) / 2.0;
              add_atom[add_count].state = "good";
              add_count++;
            }
            if(dx>lx/2.0 && dy>ly/2.0 && rxy<cut_Si_Ge && layer_dir == "zl")
            {
              if(atom[i].x + atom[j].x > 1.0)
                add_atom[add_count].x = (atom[i].x + atom[j].x - 1.0) / 2.0;
              else
                add_atom[add_count].x = (atom[i].x + atom[j].x + 1.0) / 2.0;

              if(atom[i].y + atom[j].y > 1.0)
                add_atom[add_count].y = (atom[i].y + atom[j].y - 1.0) / 2.0;
              else
                add_atom[add_count].y = (atom[i].y + atom[j].y + 1.0) / 2.0;

              add_atom[add_count].z = (atom[i].z + atom[j].z) / 2.0;
              add_atom[add_count].state = "good";
              add_count++;
            }
          }
        }
      }

      // whether to check O atoms can or can't add 
      for(int i=0; i<add_count; i++)
      {
        for(int j=0; j<all_count; j++)
        {  
          dx = fabs(add_atom[i].x - atom[j].x) * lx;
          dy = fabs(add_atom[i].y - atom[j].y) * ly;
          dz = fabs(add_atom[i].z - atom[j].z) * lz;

          r   = pow(dx, 2.0) + pow(dy, 2.0) + pow(dz, 2.0);
          rx  = pow(lx-dx, 2.0) + pow(dy, 2.0) + pow(dz, 2.0);
          ry  = pow(dx, 2.0) + pow(ly-dy, 2.0) + pow(dz, 2.0);
          rz  = pow(dx, 2.0) + pow(dy, 2.0) + pow(lz-dz, 2.0);
          rxy  = pow(lx-dx, 2.0) + pow(ly-dy, 2.0) + pow(dz, 2.0);

          if(dx<=lx/2.0 && dy<=ly/2.0 && r<cut_O && r!=0 && layer_dir == "zl")
            add_atom[i].state = "bad";
          if(dx>lx/2.0 && dy<=ly/2.0 && rx<cut_O && layer_dir == "zl")
            add_atom[i].state = "bad";
          if(dx<=lx/2.0 && dy>ly/2.0 && ry<cut_O && layer_dir == "zl")
            add_atom[i].state = "bad";
          if(dx>lx/2.0 && dy>ly/2.0 && rxy<cut_O && layer_dir == "zl")
            add_atom[i].state = "bad";
        }
        for(int j=0; j<all_count; j++)
        {  
          dx = fabs(add_atom[i].x - atom[j].x) * lx;
          dy = fabs(add_atom[i].y - atom[j].y) * ly;
          dz = fabs(add_atom[i].z - atom[j].z) * lz;

          r   = pow(dx, 2.0) + pow(dy, 2.0) + pow(dz, 2.0);
          rx  = pow(lx-dx, 2.0) + pow(dy, 2.0) + pow(dz, 2.0);
          ry  = pow(dx, 2.0) + pow(ly-dy, 2.0) + pow(dz, 2.0);
          rz  = pow(dx, 2.0) + pow(dy, 2.0) + pow(lz-dz, 2.0);
          rxy  = pow(lx-dx, 2.0) + pow(ly-dy, 2.0) + pow(dz, 2.0);

          if(dx<=lx/2.0 && dy<=ly/2.0 && r<cut_O && r!=0 && layer_dir == "zl")
            add_atom[i].state = "bad";
          if(dx>lx/2.0 && dy<=ly/2.0 && rx<cut_O && layer_dir == "zl")
            add_atom[i].state = "bad";
          if(dx<=lx/2.0 && dy>ly/2.0 && ry<cut_O && layer_dir == "zl")
            add_atom[i].state = "bad";
          if(dx>lx/2.0 && dy>ly/2.0 && rxy<cut_O && layer_dir == "zl")
            add_atom[i].state = "bad";
        }
        for(int j=0; j<all_count; j++)
        {  
          dx = fabs(add_atom[i].x - atom[j].x) * lx;
          dy = fabs(add_atom[i].y - atom[j].y) * ly;
          dz = fabs(add_atom[i].z - atom[j].z) * lz;

          r   = pow(dx, 2.0) + pow(dy, 2.0) + pow(dz, 2.0);
          rx  = pow(lx-dx, 2.0) + pow(dy, 2.0) + pow(dz, 2.0);
          ry  = pow(dx, 2.0) + pow(ly-dy, 2.0) + pow(dz, 2.0);
          rz  = pow(dx, 2.0) + pow(dy, 2.0) + pow(lz-dz, 2.0);
          rxy  = pow(lx-dx, 2.0) + pow(ly-dy, 2.0) + pow(dz, 2.0);

          if(dx<=lx/2.0 && dy<=ly/2.0 && r<cut_O && r!=0 && layer_dir == "zl")
            add_atom[i].state = "bad";
          if(dx>lx/2.0 && dy<=ly/2.0 && rx<cut_O && layer_dir == "zl")
            add_atom[i].state = "bad";
          if(dx<=lx/2.0 && dy>ly/2.0 && ry<cut_O && layer_dir == "zl")
            add_atom[i].state = "bad";
          if(dx>lx/2.0 && dy>ly/2.0 && rxy<cut_O && layer_dir == "zl")
            add_atom[i].state = "bad";
        }        for(int j=0; j<all_count; j++)
        {  
          dx = fabs(add_atom[i].x - atom[j].x) * lx;
          dy = fabs(add_atom[i].y - atom[j].y) * ly;
          dz = fabs(add_atom[i].z - atom[j].z) * lz;

          r   = pow(dx, 2.0) + pow(dy, 2.0) + pow(dz, 2.0);
          rx  = pow(lx-dx, 2.0) + pow(dy, 2.0) + pow(dz, 2.0);
          ry  = pow(dx, 2.0) + pow(ly-dy, 2.0) + pow(dz, 2.0);
          rz  = pow(dx, 2.0) + pow(dy, 2.0) + pow(lz-dz, 2.0);
          rxy  = pow(lx-dx, 2.0) + pow(ly-dy, 2.0) + pow(dz, 2.0);

          if(dx<=lx/2.0 && dy<=ly/2.0 && r<cut_O && r!=0 && layer_dir == "zl")
            add_atom[i].state = "bad";
          if(dx>lx/2.0 && dy<=ly/2.0 && rx<cut_O && layer_dir == "zl")
            add_atom[i].state = "bad";
          if(dx<=lx/2.0 && dy>ly/2.0 && ry<cut_O && layer_dir == "zl")
            add_atom[i].state = "bad";
          if(dx>lx/2.0 && dy>ly/2.0 && rxy<cut_O && layer_dir == "zl")
            add_atom[i].state = "bad";
        }
        for(int j=i+1; j<add_count; j++)
        {  
          dx = fabs(add_atom[i].x - add_atom[j].x) * lx;
          dy = fabs(add_atom[i].y - add_atom[j].y) * ly;
          dz = fabs(add_atom[i].z - add_atom[j].z) * lz;

          r   = pow(dx, 2.0) + pow(dy, 2.0) + pow(dz, 2.0);
          rx  = pow(lx-dx, 2.0) + pow(dy, 2.0) + pow(dz, 2.0);
          ry  = pow(dx, 2.0) + pow(ly-dy, 2.0) + pow(dz, 2.0);
          rz  = pow(dx, 2.0) + pow(dy, 2.0) + pow(lz-dz, 2.0);
          rxy  = pow(lx-dx, 2.0) + pow(ly-dy, 2.0) + pow(dz, 2.0);

          if(dx<=lx/2.0 && dy<=ly/2.0 && r<cut_O && r!=0 && layer_dir == "zl")
            add_atom[i].state = "bad";
          if(dx>lx/2.0 && dy<=ly/2.0 && rx<cut_O && layer_dir == "zl")
            add_atom[i].state = "bad";
          if(dx<=lx/2.0 && dy>ly/2.0 && ry<cut_O && layer_dir == "zl")
            add_atom[i].state = "bad";
          if(dx>lx/2.0 && dy>ly/2.0 && rxy<cut_O && layer_dir == "zl")
            add_atom[i].state = "bad";
        }
      }
      for(int i=0; i<add_count; i++)
      {
        if(add_atom[i].status=="yet" && add_atom[i].state=="good")
          yet_count++;
        if(add_atom[i].status=="new" && add_atom[i].state=="good")
          new_count++;
      }

      // determin number of O atoms to add in the current oxidation process
      if(times<oxi_times)
      {
        ratio = (int)((yet_count + new_count) / (oxi_times - times +1));
        if((yet_count+new_count)%(oxi_times-times+1)!=0)
          ratio++;
      }

      // choose adding O atoms according to percent of adding O atoms
      srand((unsigned int)time(NULL));
      if(times<oxi_times)
      {
        int i;
        while(add_num<ratio)
        {
          i = rand()%add_count;
          if(add_atom[i].state=="good" && add_atom[i].check==0)
          {
            add_atom[i].check = 1;
            add_num++;
            if(add_num==ratio) break;
          }
        }
      }

      // output information of all adding O atoms
      for(int i=0; i<add_count; i++)
      {        
        if(times<oxi_times && add_atom[i].state=="good" && add_atom[i].check==1)
        {
          fout << add_atom[i].name << " " << add_atom[i].x << " " << add_atom[i].y << " " << add_atom[i].z << endl;
          if(add_atom[i].status=="yet")
            add_yet_num++;
          if(add_atom[i].status=="new")
            add_new_num++;
        }
        if(times>=oxi_times && add_atom[i].state=="good")
        {
          fout << add_atom[i].name << " " << add_atom[i].x << " " << add_atom[i].y << " " << add_atom[i].z << endl;
          if(add_atom[i].status=="yet")
            add_yet_num++;
          if(add_atom[i].status=="new")
            add_new_num++;
        }
      }
      fout.close();

      if(atom_name=="Si")
      {
        atom_name="SiO";
      }
      else if(atom_name=="Ge")
      {
        atom_name="GeO";
      }



      fout_log << "~ Oxidation Start! ~" << endl << endl;
      fout_log << "add yet atom:" << add_yet_num << "/" << yet_count << endl;
      fout_log << "add new atom:" << add_new_num << "/" << new_count << endl;
      fout_log << "add all atom:" << add_yet_num + add_new_num << "/" << yet_count + new_count << endl;
      fout_log << "all:" << all_count + add_yet_num + add_new_num << "  "
               << "Si/Ge:" << Si_Ge_count << "  "
               << "O:" << O_count + add_yet_num + add_new_num<< endl << endl;
      
      fout_log << "Optimization \"" << output << "\" Start!  ";

      if(atom_name=="SiO")
      {
        os << "mdmake --convert lmpLAY -c " << output << " >stable.lmp 2>/dev/null";
        system(os.str().c_str());
        os.str(empty);
        os << "mpirun -np " << opts.omp << " lmp_custom <in.stable.SiO >/dev/null 2>&1";
        system(os.str().c_str());
        os.str(empty);
        os << "mdmake --convert mdlCLEAN --lammps SiO2 -c stable.final >opt.txt 2>/dev/null";
        system(os.str().c_str());
        os.str(empty);
        cout << "X" << flush;
        fout_log << "Finish!" << endl;
        fout_log << "MD Calculate \"opt.txt\".  Start!  ";
        os << "mdmake --convert lmpLAY -c opt.txt >md.lmp 2>/dev/null";
        system(os.str().c_str());
        os.str(empty);
        os << "mpirun -np " << opts.omp << " lmp_custom <in.md.SiO >/dev/null 2>&1";
        system(os.str().c_str());
        os.str(empty);
        os << "mdmake --convert mdlCLEAN --lammps SiO2 -c md.final >md.txt 2>/dev/null";
        system(os.str().c_str());
        os.str(empty);
        cout << "X" << flush;
        fout_log << "Finish!" << endl;
        fout_log << "Optimize \"md.txt\".       Start!  ";
        os << "mdmake --convert lmpLAY -c md.txt >stable.lmp 2>/dev/null";
        system(os.str().c_str());
        os.str(empty);
        os << "mpirun -np " << opts.omp << " lmp_custom <in.stable.SiO >/dev/null 2>&1";
        system(os.str().c_str());
        os.str(empty);
        os << "mdmake --convert mdlCLEAN --lammps SiO2 -c stable.final >" << input << " 2>/dev/null";
        system(os.str().c_str());
        os.str(empty);
        fout_log << "Finish!" << endl;
        cout << "X  Finish!" << endl;
        fout_log << "\"stable.txt\" is made!" << endl << endl;
        fout_log << "~ Oxidation Finish! ~" << endl;
        fout_log << "****************************************" << endl << endl;
      }
      else if(atom_name=="GeO")
      {
        os << "tail -n +3 " << output << " | cat head.stable - >stable.mdl";
        system(os.str().c_str());
        os.str(empty);
        os << "mdlabo < stable.mdl > stable.final";
        system(os.str().c_str());
        os.str(empty);
        os << "mdmake --convert mdlCLEAN -c stable.final -l >opt.txt 2>/dev/null";
        system(os.str().c_str());
        os.str(empty);
        cout << "X" << flush;
        fout_log << "Finish!" << endl;
        fout_log << "MD Calculate \"opt.txt\".  Start!  ";
        os << "tail -n +3 opt.txt | cat head.md - >md.mdl";
        system(os.str().c_str());
        os.str(empty);
        os << "mdlabo < md.mdl > md.final";
        system(os.str().c_str());
        os.str(empty);
        os << "mdmake --convert mdlCLEAN -c md.final -l >md.txt 2>/dev/null";
        system(os.str().c_str());
        os.str(empty);
        cout << "X" << flush;
        fout_log << "Finish!" << endl;
        fout_log << "Optimize \"md.txt\".       Start!  ";
        os << "tail -n +3 md.txt | cat head.stable - >stable.mdl";
        system(os.str().c_str());
        os.str(empty);
        os << "mdlabo < stable.mdl > stable.final";
        system(os.str().c_str());
        os.str(empty);
        os << "mdmake --convert mdlCLEAN -c stable.final -l >" << input << " 2>/dev/null";
        system(os.str().c_str());
        os.str(empty);
        fout_log << "Finish!" << endl;
        cout << "X  Finish!" << endl;
        fout_log << "\"stable.txt\" is made!" << endl << endl;
        fout_log << "~ Oxidation Finish! ~" << endl;
        fout_log << "****************************************" << endl << endl;
      }

      // copy file of each process
      os << "cp ./add_O.txt ./layer" << oxi_lay << "/add_O_" << times << "_" << oxi_times << "-" << oxi_check << ".txt";
      system(os.str().c_str());
      os.str(empty);
      os << "cp ./opt.txt ./layer" << oxi_lay << "/opt" << times << "_" << oxi_times << "-" << oxi_check << ".txt";
      system(os.str().c_str());
      os.str(empty);
      os << "cp ./md.txt ./layer" << oxi_lay << "/md" << times << "_" << oxi_times << "-" << oxi_check << ".txt";
      system(os.str().c_str());
      os.str(empty);
      os << "cp ./stable.txt ./layer" << oxi_lay << "/stable" << times << "_" << oxi_times << "-" << oxi_check << ".txt";
      system(os.str().c_str());
      os.str(empty);
      os << "cp ./oxi.xyz ./layer" << oxi_lay << "/oxi" << times << "_" << oxi_times << "-" << oxi_check << ".xyz";
      system(os.str().c_str());
      os.str(empty);
      os << "rm ./oxi.xyz";
      system(os.str().c_str());
      os.str(empty);

      // free dynamic memory
      delete[] add_atom;
    }
    fout_log.close();                  // close the file
    os << "cp ./log.txt ./layer" << oxi_lay << "/log.txt";
    system(os.str().c_str());
    os.str(empty);
  }

  return 0;
}

// round off the number to "digit" places of decimals
double round(double round, double digit)
{
  return (double)((int)(round * pow(10.0, digit) + 0.5)) / pow(10.0, digit);
}

