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
int cluster(int argc, char **argv, optpara opts){
	string empty = "";    // initialize stringstream
	ostringstream os;
        string path  = "../";   // directory of mdmake
	
	//mdmake -c SiGe -r 0.5 -n "30 4 4" >> start.mdl
	//mdmake --cluster pxyz -n "" -c start.mdl >> goal.mdl
	  //mdmake --swap random -n "7 1 100" -c 0.mdl >> 1.mdl

        //Swap Condition
	int option = 1, swap_loop = 2, swap_step = 2;
	if(opts.cellnum.length() > 0) {
	  istringstream cell_iss(opts.cellnum.c_str());
	  cell_iss >> option >> swap_loop >> swap_step;
 	  if(option == 0 || swap_loop == 0 || swap_step == 0) {
 	    cerr << "Cell number is wrong" << endl;
 	    return 1;
 	  }
 	}
	
        //Periodic Boundary Condition   [1(periodic) or 0(noneriodic)}
        int p_x = 0, p_y = 0, p_z = 0;
	if(opts.periodic == "pxyz")p_x=1, p_y=1, p_z=1;
	else if(opts.periodic == "px")p_x=1;
	else if(opts.periodic == "py")p_y=1;
	else if(opts.periodic == "pz")p_z=1;
	else if(opts.periodic == "pxy")p_x=1, p_y=1;
	else if(opts.periodic == "pxz")p_x=1, p_z=1;
	else if(opts.periodic == "pyz")p_y=1, p_z=1;
	else if(opts.periodic == "none"){

	}
	else{
	    cerr << "\"" << opts.periodic << "\" is undefined!" << endl;
 	    return 1;
	}

	//0
        os << "mkdir data_" << opts.atomcomp << "_";
	system(os.str().c_str());
	os.str(empty);        

	os << path << "mdmake --convert mdl -c " << opts.atomcomp << " > 0.mdl";
	system(os.str().c_str());
	os.str(empty);
	os << path << "mdmake --convert mdlGROUP -n \"" << 2-p_x << " " << 2-p_y << " "  << 2-p_z << "\" -c 0.mdl > 0.grp";
        system(os.str().c_str());
	os.str(empty);
	os << path << "mdmake --convert xyz -c 0.mdl > 0.xyz";
	system(os.str().c_str());
	os.str(empty);
	os << "mv 0.grp data_" << opts.atomcomp << "_";
	system(os.str().c_str());
	os.str(empty);
	os << "mv 0.mdl.cdist data_" << opts.atomcomp << "_";
	system(os.str().c_str());
	os.str(empty);
	os << "mv 0.xyz data_" << opts.atomcomp << "_";
	system(os.str().c_str());
	os.str(empty);

	//loop(i-1>1)
        for(int i=1; i<=swap_loop; i++){
	  os << path << "mdmake --swap random -n \"" << p_x*1+p_y*2+p_z*4 << " 1 " << swap_step << "\" -c" << i-1 << ".mdl >"  << i << ".mdl";
	  system(os.str().c_str());
	  os.str(empty);
          if(i%1==0){
	    os << path << "mdmake --convert mdlGROUP -n \"" << 2-p_x << " " << 2-p_y << " "  << 2-p_z << "\" -c " << i << ".mdl > "  << i << ".grp";
	    system(os.str().c_str());
	    os.str(empty);
	    os << path << "mdmake --convert xyz -c " << i << ".mdl > " << i << ".xyz";
	    system(os.str().c_str());
	    os.str(empty);
	    os << "mv " << i << ".mdl.cdist data_" << opts.atomcomp << "_";
	    system(os.str().c_str());
	    os.str(empty);
	    os << "mv " << i << ".grp data_" << opts.atomcomp << "_";
	    system(os.str().c_str());
	    os.str(empty);
	    os << "mv " << i << ".xyz data_" << opts.atomcomp << "_";
	    system(os.str().c_str());
	    os.str(empty);
          }
	}

	//end
        os << "cp " << swap_loop << ".mdl cluster_" << swap_step*swap_loop << ".mdl";
	system(os.str().c_str());
	os.str(empty);
        for(int i=0; i<=swap_loop; i++){
	    os << "mv " << i << ".mdl data_" << opts.atomcomp << "_";
	    system(os.str().c_str());
	    os.str(empty);
        }

  return 0;
}

