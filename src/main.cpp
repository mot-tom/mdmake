#include <unistd.h>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <string>
#include "opt.h"
using namespace std;

static option options[] =
  {
    {"help", no_argument, NULL, 'h'},
    {"composition", required_argument, NULL, 'c'},
    {"compratio", required_argument, NULL, 'r'},
    {"seed", required_argument, NULL, 'z'},
    {"stressed", required_argument, NULL, 's'},
    {"strained", required_argument, NULL, 'e'},
    {"potential", required_argument, NULL, 'p'},
    {"cellnum", required_argument, NULL, 'n'},
    {"layer", optional_argument, NULL, 'l'},
    {"wire", optional_argument, NULL, 'w'},
    {"freepbc", no_argument, NULL, 'f'},
    {"threed", required_argument, NULL, 't'},
    {"lammps", required_argument, NULL, 'v'},
    {"convert", required_argument, NULL, 'C'},
    {"convert-org", required_argument, NULL, 'O'},
    {"oxi", no_argument, NULL, 'X'},
    {"omp", required_argument, NULL, 'M'},
    {"stack", required_argument, NULL, 'S'},
    {"swap", required_argument, NULL, 'W'}, //takizawa
    {"cluster", required_argument, NULL, 'L'}, //takizawa
    {"rdf", required_argument, NULL, 'R'}, //takizawa
    {0, 0, 0, 0}
  };

int confinementdia(int iargc, char* argv[], optpara opts);
int tetragonaldia(int iargc, char* argv[], optpara opts);
int orthodia(int iargc, char* argv[], optpara opts);
int orthodia2(int iargc, char* argv[], optpara opts);
int quartz(int argc, char **argv, optpara opts);
int tridymite(int argc, char **argv, optpara opts);
int rutile(int argc, char **argv, optpara opts);
int rocksalt(int argc, char **argv, optpara opts);
int corundum(int argc, char **argv, optpara opts);
int wurtzite(int argc, char **argv, optpara opts);
int dichalcogenide(int argc, char **argv, optpara opts);
int zincblende(int iargc, char* argv[], optpara opts);
int diamond(int iargc, char* argv[], optpara opts);
int fcc(int iargc, char* argv[], optpara opts);

int oxi_amo(int iargc, char* argv[], optpara opts);
int oxi_layer(int iargc, char* argv[], optpara opts);
int diamond_layer(int iargc, char* argv[], optpara opts);
int diamond_oxide(int iargc, char* argv[], optpara opts);
int diamond_comp(int iargc, char* argv[], optpara opts);
int tetra_layer(int iargc, char* argv[], optpara opts);
int tetra_comp(int iargc, char* argv[], optpara opts);
int wurtz_layer(int iargc, char* argv[], optpara opts);
int ortho_layer(int iargc, char* argv[], optpara opts);
int ortho2_layer(int iargc, char* argv[], optpara opts);
int super_lattice(int iargc, char* argv[], optpara opts);

int diamond_lmp(int iargc, char* argv[], optpara opts, string lmpname);

int stack(int iargc, char* argv[], optpara opts, string lmpname);
int convert(int iargc, char* argv[], optpara opts, string lmpname);
int convert_org(int iargc, char* argv[], optpara opts);
int rdf(int iargc, char* argv[], optpara opts); //takizawa

int swap(int iargc, char* argv[], optpara opts); //takizawa
int cluster(int iargc, char* argv[], optpara opts); //takizawa


extern int shape(string shape, double x, double y, double z);
double sqrwav(double w, double d);
double triwav(double w, double d);

int main(int argc, char **argv)
{
       cerr << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
       cerr << "                Crystal maker for mdlabo version 2.00" << endl << endl; //takizawa
       cerr << "                Motohiro Tomita & Junya Takizawa 2022" << endl;
       cerr << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

  string dir;
  int opt, index;
  int disp = 0;
  int lmp = 0;
  int conv = 0;
  int conv_o = 0;
  int oxi = 0;
  int stk = 0;
  int wap = 0; //takizawa
  int cls = 0; //takizawa
  int rdf_i = 0; //takizawa
  string lmpname = "";
  optpara opts;
  opts.atomcomp = "Si";
  opts.atomratio = "";
  opts.seed = 10;
  opts.stress = "";
  opts.strain = "";
  opts.potential = "ESWGE";
  opts.cellnum = "";
  opts.pbc_x = "";
  opts.pbc_y = "";
  opts.pbc_z = "";
  opts.shape = "cuboid";
  opts.omp = 1;
  opts.stack = "";
  opts.method = ""; //takizawa
  opts.periodic = ""; //takizawa

  

  while((opt = getopt_long(argc, argv, "hc:r:z:s:e:p:n:l::w::ft:v:C:O:XM:S:", options, &index)) != -1){
    switch(opt){
      case 'h':
        cerr << "Help of options." << endl;
        cerr << "        -c: Input atom composition." << endl;
        cerr << "        -r: Input atom ratio." << endl;
        cerr << "    --seed: Input seed of rand." << endl;
        cerr << "        -s: Induce stress. [GPa]" << endl;
        cerr << "        -e: Induce strain. [%]" << endl;
        cerr << "        -p: Input md potential." << endl;
        cerr << "        -n: Input cell number. ex. \"4 4 4\"" << endl;
        cerr << "        -l: PBC layered." << endl;
        cerr << "        -w: PBC wired." << endl;
        cerr << "        -f: PBC free." << endl;
        cerr << "        -t: Input 3D shape." << endl;
        cerr << "  --lammps: Output for lammps." << endl;
        cerr << " --convert: Convert file type with -c option." << endl;
        cerr << "             ex. --convert lmp -c <filename>" << endl;
        cerr << "   --stack: Stacking 2 structure files to x direction." << endl;
        cerr << "             ex. --stack <filename 1> -c <filename 2>" << endl;
        cerr << "     --oxi: Oxidation program." << endl;
        cerr << "     --omp: Number of OpenMP for optimizing MD calculation." << endl;
        cerr << "    --swap: Swap coordinations of two atoms." << endl; //takizawa
        cerr << " --cluster: make atomic clusters by compounds." << endl; //takizawa
        cerr << endl << "Atom compositions for option -c." << endl;
        cerr << "           Diamond: Si, Ge, Sn, C, SiGe, SiC, GeSiSn" << endl;
        cerr << "                    3H-Si, 3H-SiGe-alloy, 3H-GeSiSn-alloy" << endl;
        cerr << "                    3T-Si, 3T-SiGe-alloy, 3T-GeSiSn-alloy" << endl;
        cerr << "                    3HO-Si, 3HO-SiGe-alloy, 3HO-GeSiSn-alloy" << endl;
        cerr << "                    3TO-Si, 3TO-SiGe-alloy, 3TO-GeSiSn-alloy" << endl;
        cerr << "                    Si-in-Ge, Ge-in-Si" << endl;
        cerr << "        Zincblende: 3C-SiC, 3C-SiGe, 3H-SiC, GaAs, InAs, 3C-BN" << endl;
        cerr << "          Wurtzite: 2H-SiC, 4H-SiC, 6H-SiC, 8H-SiC, 10H-SiC," << endl;
        cerr << "                    12H-SiC, 18H-SiC, 15R-SiC, 2H-Si, 4H-Si" << endl;
        cerr << "                    2H-AlN, 2H-BN, 2H-BeO" << endl;
        cerr << "     Diamond oxide: bc-SiO2, bc2-SiO2, bc3-SiO2, ac-SiO2, " << endl;
        cerr << "                    bc-GeO2, bc2-GeO2, bc3-GeO2, ac-GeO2" << endl;
        cerr << "      Quartz oxide: bq-SiO2, aq-SiO2, bq-GeO2, aq-GeO2" << endl;
        cerr << "   Tridymite oxide: bt-SiO2, at-SiO2, bt-GeO2, at-GeO2" << endl;
        cerr << "      Rutile oxide: st-SiO2, st-GeO2, ru-SnO2, ru-TiO2" << endl;
        cerr << "  Diamond compound: g-Si3N4, C-ZrO2, C-SnO2, C-SiO2, Y2O3" << endl;
        cerr << "                    an-TiO2, T-ZrO2" << endl;
        cerr << " Rocksalt compound: NaCl, MgO, SrO, TiN, GeSbTe" << endl; //takizawa
        cerr << "         Corundume: a-Al2O3" << endl;
        cerr << "     Super lattice: sl1-SiGe (stack every 1 layer)" << endl;
        cerr << "                    sl2-SiGe (stack every 2 layers)" << endl;
        cerr << "                    slr-SiGe (stack layer randomly)" << endl;
        cerr << "         Amorphous: am-Si, am-Ge" << endl;
        cerr << "                    am-aq-SiO2, am-bq-SiO2, am-C-SiO2" << endl;
        cerr << "                    am-at-SiO2, am-bt-SiO2, am-st-SiO2, " << endl;
        cerr << "                    am-ac-SiO2, am-bc-SiO2, am-bc-GeO2" << endl;
        cerr << "                    am-ru-TiO2, am-a-Al2O3" << endl;
        cerr << "                    am-MgO, am-SrO" << endl;
        cerr << " Oxidation program: zl-SiO2, zl-3H-SiO2, zl-3T-SiO2" << endl;
        cerr << "                    zl-3HO-SiO2, zl-GeO2 (both surfaces)" << endl;
        cerr << "                    zl1-SiO2, zl1-3T-SiO2, zl1-3H-SiO2" << endl;
        cerr << "                    zl1-3HO-SiO2 (only front side)" << endl;
        cerr << "    Dichalcogenide: MoS2" << endl;
        cerr << "             Metal: b-Sn, Ni, Au" << endl;
        cerr << endl << "3D shapes for option -t." << endl;
        cerr << "         sphere: Sphere shape." << endl;
        cerr << "       cylinder: Cylinder to x direction." << endl;
        cerr << "     cylinder_y: Cylinder to y direction." << endl;
        cerr << "     cylinder_z: Cylinder to z direction." << endl;
        cerr << " wavedcylinder2: Waved cylinder of period 2 to x direction." << endl;
        cerr << " wavedcylinder4: Waved cylinder of period 4 to x direction." << endl;
        cerr << " sawedcylinder2: Sawed cylinder of period 2 to x direction." << endl;
        cerr << " sawedcylinder4: Sawed cylinder of period 4 to x direction." << endl;
        cerr << "    wavedprism2: Waved prism of period 2 to x direction." << endl;
        cerr << "    wavedprism4: Waved prism of period 4 to x direction." << endl;
        cerr << "    sawedprism2: Sawed prism of period 2 to x direction." << endl;
        cerr << "    sawedprism4: Sawed prism of period 4 to x direction." << endl;
        cerr << endl << "Convert output file types for option --convert." << endl;
        cerr << "         mdl: Structure file for Mdlabo." << endl;
        cerr << "     mdlZtoX: Structure file for Mdlabo with exchanging Z and X axis." << endl;
        cerr << "     mdlYtoX: Structure file for Mdlabo with exchanging Y and X axis." << endl;
        cerr << "      mdlCUT: Structure file for Mdlabo" << endl;
        cerr << "              with deleting atoms which are out of simulation box." << endl;
        cerr << "     mdlFILM: Structure file for Mdlabo with converting thin film." << endl;
        cerr << "     mdlWIRE: Structure file for Mdlabo with converting thin wire." << endl;
        cerr << "    mdlCLEAN: Structure file for Mdlabo with deleting atoms away from film." << endl;
        cerr << "     mdlZRES: Structure file for Mdlabo with resizing." << endl; //takizawa
        cerr << "    mdlGROUP: Structure file for Mdlabo with grouping." << endl; //takizawa
        cerr << "mdlINTERFACE: Structure file for Mdlabo with deleting atoms near by interface." << endl; //takizawa
        cerr << "         lmp: Structure file for LAMMPS." << endl;
        cerr << "       lmpAr: Structure file for LAMMPS with adding Ar atom type." << endl;
        cerr << "      lmpLAY: Structure file for LAMMPS with oxidation layer property." << endl;
        cerr << "         xyz: Structure file for VMD." << endl;
        cerr << "         alm: Structure file for Alamode." << endl;
        cerr << "        vasp: Structure file for VASP." << endl;
        cerr << endl << "Convert input file types for option --convert and --stack." << endl;
        cerr << endl << "Calculate RDF:--rdf rdf   Calculate Coordination number:--rdf num." << endl; //takizawa
        cerr << " Mdlabo input file, Mdlabo output file, LAMMPS output file" << endl;
        cerr << endl;
        return 0;
      case 'c':
        opts.atomcomp = optarg;
        break;
      case 'S':
        opts.stack = optarg;
        stk = 1;
        break;
      case 'W': //takizawa
        opts.method = optarg;
        wap = 1;
        break;
      case 'L': //takizawa
        opts.periodic = optarg;
        cls = 1;
        break;
      case 'R': //takizawa
        opts.atomratio = optarg;
        rdf_i = 1;
        break;
      case 'r':
        opts.atomratio = optarg;
        break;
      case 'z':
        opts.seed = atoi(optarg);
        break;
      case 's':
        opts.stress = optarg;
        break;
      case 'e':
        opts.strain = optarg;
        break;
      case 'p':
        opts.potential = optarg;
        break;
      case 'n':
        opts.cellnum = optarg;
        break;
      case 'l':
        if(optarg) dir = optarg;
        if(dir == "x" || dir == "X")
          opts.pbc_x = " N";
        else if(dir == "y" || dir == "Y")
          opts.pbc_y = " N";
        else
          opts.pbc_z = " N";
        break;
      case 'w':
        if(optarg) dir = optarg;
        if(dir == "y" || dir == "Y")
          { opts.pbc_x = " N"; opts.pbc_z = " N"; }
        else if(dir == "z" || dir == "Z")
          { opts.pbc_x = " N"; opts.pbc_y = " N"; }
        else
          { opts.pbc_y = " N"; opts.pbc_z = " N"; }
        break;
      case 'f':
        opts.pbc_x = " N";
        opts.pbc_y = " N";
        opts.pbc_z = " N";
        break;
      case 't':
        opts.shape = optarg;
        break;
      case 'v':
        lmpname = optarg;
        lmp = 1;
        break;
      case 'C':
        opts.atomratio = optarg;
        conv = 1;
        break;
      case 'O':
        opts.atomratio = optarg;
        conv_o = 1;
        break;
      case 'X':
        oxi = 1;
        break;
      case 'M':
        opts.omp = atoi(optarg);
        break;
      default:
        cerr << "Wrong option" << endl;
        return 1;
        break;
     }
  }

  if(argv[optind]) {
    cerr << "Wrong option" << endl;
    return 1;
  }
  if(opts.stress.length() > 0 && opts.strain.length() > 0) {
    cerr << "Wrong option" << endl;
    return 1;
  }

  if(conv) convert(argc, argv, opts, lmpname);
  else if(stk) stack(argc, argv, opts, lmpname);
  else if(conv_o) convert_org(argc, argv, opts);
  else if(oxi) oxi_layer(argc, argv, opts);
  else if(wap) swap(argc, argv, opts); //takizawa
  else if(cls) cluster(argc, argv, opts); //takizawa
  else if(rdf_i) rdf(argc, argv, opts); //takizawa
  else if((opts.atomcomp == "Si" || opts.atomcomp == "Ge" || opts.atomcomp == "C" || opts.atomcomp == "Sn"
  || opts.atomcomp == "SiGe" || opts.atomcomp == "SiC" || opts.atomcomp == "GeSiSn")
  && lmp)
    diamond_lmp(argc, argv, opts, lmpname);
  else if(opts.atomcomp == "Si" || opts.atomcomp == "Ge" || opts.atomcomp == "C" || opts.atomcomp == "Sn"
  || opts.atomcomp == "SiGe" || opts.atomcomp == "SiC" || opts.atomcomp == "GeSiSn")
    diamond(argc, argv, opts);
  else if(opts.atomcomp == "slr-SiGe" || opts.atomcomp == "sl1-SiGe" || opts.atomcomp == "sl2-SiGe")
    super_lattice(argc, argv, opts);
  else if(opts.atomcomp == "zl-3H-Si" || opts.atomcomp == "zl1-3H-Si")
    wurtz_layer(argc, argv, opts);
  else if(opts.atomcomp == "zl-3T-Si" || opts.atomcomp == "zl1-3T-Si")
    tetra_layer(argc, argv, opts);
  else if(opts.atomcomp == "zl-3HO-Si" || opts.atomcomp == "zl1-3HO-Si")
    ortho_layer(argc, argv, opts);
  else if(opts.atomcomp == "zl-3TO-Si" || opts.atomcomp == "zl1-3TO-Si")
    ortho2_layer(argc, argv, opts);
  else if(opts.atomcomp == "zl-Si" || opts.atomcomp == "zl1-Si" || opts.atomcomp == "zl-Ge")
    diamond_layer(argc, argv, opts);
  else if(opts.atomcomp == "zl-SiO2" || opts.atomcomp == "zl1-SiO2" || opts.atomcomp == "zl-GeO2"
  || opts.atomcomp == "zl-3H-SiO2" || opts.atomcomp == "zl1-3H-SiO2"
  || opts.atomcomp == "zl-3T-SiO2" || opts.atomcomp == "zl1-3T-SiO2"
  || opts.atomcomp == "zl-3HO-SiO2" || opts.atomcomp == "zl1-3HO-SiO2"
  || opts.atomcomp == "zl-3TO-SiO2" || opts.atomcomp == "zl1-3TO-SiO2")
    oxi_layer(argc, argv, opts);
  else if(opts.atomcomp == "am-Si" || opts.atomcomp == "am-Ge"
  || opts.atomcomp == "am-aq-SiO2" || opts.atomcomp == "am-bq-SiO2"
  || opts.atomcomp == "am-at-SiO2" || opts.atomcomp == "am-bt-SiO2" || opts.atomcomp == "am-st-SiO2"
  || opts.atomcomp == "am-ac-SiO2" || opts.atomcomp == "am-bc-SiO2" || opts.atomcomp == "am-bc-GeO2"
  || opts.atomcomp == "am-ru-TiO2" || opts.atomcomp == "am-a-Al2O3"
  || opts.atomcomp == "am-MgO" || opts.atomcomp == "am-SrO")
    oxi_amo(argc, argv, opts);
  else if(opts.atomcomp == "b-Sn" || opts.atomcomp == "T-ZrO2" || opts.atomcomp == "an-TiO2")
    tetra_comp(argc, argv, opts);
  else if(opts.atomcomp == "g-Si3N4" || opts.atomcomp == "Y2O3" || opts.atomcomp == "C-ZrO2"
  || opts.atomcomp == "C-SnO2")
    diamond_comp(argc, argv, opts);
  else if(opts.atomcomp == "bc-SiO2" || opts.atomcomp == "bc2-SiO2" || opts.atomcomp == "bc3-SiO2" || opts.atomcomp == "ac-SiO2"
  || opts.atomcomp == "bc-GeO2" || opts.atomcomp == "bc2-GeO2"|| opts.atomcomp == "bc3-GeO2" || opts.atomcomp == "ac-GeO2")
    diamond_oxide(argc, argv, opts);
  else if(opts.atomcomp == "bq-SiO2" || opts.atomcomp == "aq-SiO2"
  || opts.atomcomp == "bq-GeO2" || opts.atomcomp == "aq-GeO2")
    quartz(argc, argv, opts);
  else if(opts.atomcomp == "bt-SiO2" || opts.atomcomp == "at-SiO2"
  || opts.atomcomp == "bt-GeO2" || opts.atomcomp == "at-GeO2")
    tridymite(argc, argv, opts);
  else if(opts.atomcomp == "st-SiO2" || opts.atomcomp == "st-GeO2"
  || opts.atomcomp == "ru-SnO2" || opts.atomcomp == "ru-TiO2")
    rutile(argc, argv, opts);
  else if(opts.atomcomp == "NaCl" || opts.atomcomp == "MgO" || opts.atomcomp == "SrO" || opts.atomcomp == "TiN" || opts.atomcomp == "GeSbTe") //takizawa
    rocksalt(argc, argv, opts);
  else if(opts.atomcomp == "a-Al2O3")
    corundum(argc, argv, opts);
  else if(opts.atomcomp == "3C-SiC" || opts.atomcomp == "3C-SiCX"
       || opts.atomcomp == "3C-SiGe" || opts.atomcomp == "GaAs" || opts.atomcomp == "InAs" || opts.atomcomp == "3C-BN")
    zincblende(argc, argv, opts);
  else if(opts.atomcomp == "2H-SiC" || opts.atomcomp == "3H-SiC" || opts.atomcomp == "4H-SiC" || opts.atomcomp == "6H-SiC"
       || opts.atomcomp == "8H-SiC" || opts.atomcomp == "10H-SiC"|| opts.atomcomp == "12H-SiC"|| opts.atomcomp == "18H-SiC"
       || opts.atomcomp == "15R-SiC"|| opts.atomcomp == "2H-SiCX"
       || opts.atomcomp == "2H-Si"  || opts.atomcomp == "3H-Si"  || opts.atomcomp == "4H-Si"
       || opts.atomcomp == "3H-SiGe-alloy" || opts.atomcomp == "3H-GeSiSn-alloy"
       || opts.atomcomp == "2H-AlN" || opts.atomcomp == "2H-BN"  || opts.atomcomp == "2H-BeO")
    wurtzite(argc, argv, opts);
  else if(opts.atomcomp == "3T-Si"
       || opts.atomcomp == "3T-SiGe-alloy" || opts.atomcomp == "3T-GeSiSn-alloy")
    tetragonaldia(argc, argv, opts);
  else if(opts.atomcomp == "3HO-Si"
       || opts.atomcomp == "3HO-SiGe-alloy" || opts.atomcomp == "3HO-GeSiSn-alloy")
    orthodia(argc, argv, opts);
  else if(opts.atomcomp == "3TO-Si"
       || opts.atomcomp == "3TO-SiGe-alloy" || opts.atomcomp == "3TO-GeSiSn-alloy")
    orthodia2(argc, argv, opts);
  else if(opts.atomcomp == "Si-in-Ge" || opts.atomcomp == "Ge-in-Si")
    confinementdia(argc, argv, opts);
  else if(opts.atomcomp == "MoS2" || opts.atomcomp == "MoSX" || opts.atomcomp == "MoXS2")
    dichalcogenide(argc, argv, opts);
  else if(opts.atomcomp == "Ni" || opts.atomcomp == "Au")
    fcc(argc, argv, opts);
  else {
    cerr << "Wrong atom composition" << endl;
    return 1;
  }
  return 0;
}

extern int shape(string shape, double x, double y, double z){
  if(shape == "sphere" && ((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) + (z-0.5)*(z-0.5)) > 0.25) return 0;
  if(shape == "cylinder" && ((y-0.5)*(y-0.5) + (z-0.5)*(z-0.5)) > 0.25) return 0;
  if(shape == "cylinder_y" && ((x-0.5)*(x-0.5) + (z-0.5)*(z-0.5)) > 0.25) return 0;
  if(shape == "cylinder_z" && ((y-0.5)*(y-0.5) + (x-0.5)*(x-0.5)) > 0.25) return 0;

  //sawed wire
  if(shape == "wavedcylinder2" && ((y-0.5)*(y-0.5) + (z-0.5)*(z-0.5)) > 0.25-(1-fabs(cos(2*M_PI*x)))*0.09) return 0;
  if(shape == "wavedcylinder4" && ((y-0.5)*(y-0.5) + (z-0.5)*(z-0.5)) > 0.25-(1-fabs(cos(4*M_PI*x)))*0.09) return 0;
  if(shape == "sawedcylinder2" && ((y-0.5)*(y-0.5) + (z-0.5)*(z-0.5)) > 0.25-(fabs(triwav(2*M_PI*x,0.5)))*0.09) return 0;
  if(shape == "sawedcylinder4" && ((y-0.5)*(y-0.5) + (z-0.5)*(z-0.5)) > 0.25-(fabs(triwav(4*M_PI*x,0.5)))*0.09) return 0;
  if(shape == "wavedprism2" && ((y-0.5)*(y-0.5)) > 0.25-(1-fabs(cos(2*M_PI*x)))*0.09) return 0;
  if(shape == "wavedprism4" && ((y-0.5)*(y-0.5)) > 0.25-(1-fabs(cos(4*M_PI*x)))*0.09) return 0;
  if(shape == "sawedprism2" && ((y-0.5)*(y-0.5)) > 0.25-(fabs(triwav(2*M_PI*x,0.5)))*0.09) return 0;
  if(shape == "sawedprism4" && ((y-0.5)*(y-0.5)) > 0.25-(fabs(triwav(4*M_PI*x,0.5)))*0.09) return 0;

  //winding wire
  if(shape == "sinwinding1" && (y > cos(2*M_PI*x+M_PI)*0.5/2+0.75 || y < cos(2*M_PI*x+M_PI)*0.5/2+0.25)) return 0;
  if(shape == "sinwinding2" && (y > cos(4*M_PI*x+M_PI)*0.5/2+0.75 || y < cos(4*M_PI*x+M_PI)*0.5/2+0.25)) return 0;
  if(shape == "triwinding1" && (y > triwav(2*M_PI*x-M_PI/2,0.5)*0.5/2+0.75 || y < triwav(2*M_PI*x-M_PI/2,0.5)*0.5/2+0.25)) return 0;
  if(shape == "triwinding2" && (y > triwav(4*M_PI*x-M_PI/2,0.5)*0.5/2+0.75 || y < triwav(4*M_PI*x-M_PI/2,0.5)*0.5/2+0.25)) return 0;
  if(shape == "sqrwinding1" && (y > sqrwav(2*M_PI*x-M_PI*0.26,0.74)*0.5/2+0.75 || y < sqrwav(2*M_PI*x-M_PI*0.74,0.26)*0.5/2+0.25)) return 0;
  if(shape == "sqrwinding2" && (y > sqrwav(4*M_PI*x-M_PI*0.26,0.74)*0.5/2+0.75 || y < sqrwav(4*M_PI*x-M_PI*0.74,0.26)*0.5/2+0.25)) return 0;

  //normalized wire
  if(shape == "cuboid_xn" && ((y-0.5)*(y-0.5) > 0.25/4.0 || (z-0.5)*(z-0.5) > 0.25/4.0)) return 0;
  if(shape == "cuboid_yn" && ((x-0.5)*(x-0.5) > 0.25/4.0 || (z-0.5)*(z-0.5) > 0.25/4.0)) return 0;
  if(shape == "cuboid_zn" && ((x-0.5)*(x-0.5) > 0.25/4.0 || (y-0.5)*(y-0.5) > 0.25/4.0)) return 0;
  if(shape == "cylinder_xn" && ((y-0.5)*(y-0.5) + (z-0.5)*(z-0.5)) > 0.25/M_PI) return 0;
  if(shape == "cylinder_yn" && ((x-0.5)*(x-0.5) + (z-0.5)*(z-0.5)) > 0.25/M_PI) return 0;
  if(shape == "cylinder_zn" && ((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5)) > 0.25/M_PI) return 0;
  if(shape == "prism_xn" && (z-0.3 > y || z-0.3 > (1-y) || z < 0.3)) return 0;
  if(shape == "prism_yn" && (z-0.3 > x || z-0.3 > (1-x) || z < 0.3)) return 0;
  if(shape == "prism_zn" && (x-0.3 > y || x-0.3 > (1-y) || x < 0.3)) return 0;
  return 1;
}

double sqrwav(double w, double d){
  double dpi = 2*M_PI*d;
  while(w>2*M_PI) w-=2*M_PI;
  while(w<0) w+=2*M_PI;
  if(w<dpi) return 1;
  else return -1;
}

double triwav(double w, double d){
  double dpi = 2*M_PI*d;
  while(w>2*M_PI) w-=2*M_PI;
  while(w<0) w+=2*M_PI;
  if(w<dpi/2) return 2*w/dpi;
  else if(w<(dpi/2+M_PI)) return 2*d-2*w/M_PI+1;
  else return w/(M_PI-dpi/2)-2/(1-d);
}


