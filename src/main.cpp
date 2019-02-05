/*#ifdef _WIN32
  #include <direct.h>
  #define GetCurrentDir _getcwd
  #define osParam 0	
#elif defined _WIN64
	#include <direct.h>
	#define GetCurrentDir _getcwd
	#define osParam 0	
#else
  #include <unistd.h>
  #define GetCurrentDir getcwd
  #define osParam 1
#endif*/

#include <iostream>
#include <time.h>
#include <fstream>
#include <stdio.h>
#include <string>
#include "vector3D.h"
#include "molecule2.h"
#include "monteTraj.h"
#include "diffEqSolver.h"
#include "myStructs.h"
#include "maxBoltz.h" //delete me
#include "fileParams.h"
#include <cmath>
#include <sstream>
#include <omp.h>
using namespace std;


//define external variables
ofstream fileStream_ext;
string omegaFilePath_ext;
string fileLocaiton_ext; //TYPO
string molFilePath_ext;
string molFileName_ext;
string atomFilePath_ext;
string dataFilePath_ext;
string xyzFilePath_ext;
string iterFilePath_ext;
string name_ext;
string subDir_ext;
string printMolName_ext;
double temp_ext;
double dt_ext;
double mass_ext;
double mu;
int trajNum_ext;
int iterations_ext;
int numThreads_ext;
long long int seed_ext;
int printRate_ext;
int maxClusterSize_ext;
int multipole_order_ext;
double cutOffRadius_ext;
double gMax_ext;
double gMin_ext;
double totalCharge_ext;
bool displayProg_ext;

double dispersion_radius_ext;
double multipole_radius_ext;

bool multipole_ext;
bool dispersionCutoff_ext;
bool charge_ext;
//bool orbitPi_ext;
bool projection_ext; //ellipsoid projection on/off
bool shorten_ext; //shorten trajectories
bool printTraj_ext; //print xyz files of trajectories
bool printTrajOnFail_ext; //print trajectory only if they fail = true, otherwise it prints all successful ones
bool printData_ext; //print starting data
bool printDataOnFail_ext; //print data only if they fail = true, otherwise it prints all successful ones
bool printMol_ext;
bool printToColsole_ext = false;

//const double avo_ext; //Avogadro constant
//const double calorie_ext; //number of calories in 1 jouel
//const double protonMassKg_ext; //proton mass in kg
//const double heliumMassU_ext; //atomic mass of helium n a.u.


//debugging only
bool disp = true;
int numThreads;
monteTraj *monte;
initParams params;
threadSharedParams shared;
ofstream iterData;

//double fixedVelocity_ext = 0; //debugging only
int main(int argc, char* argv[]) {

	fileParams files;
	char buffer[100];
	char printOut[1000] = "";

        //get current working directory
        char currentPath[FILENAME_MAX];
        if (!GetCurrentDir(currentPath, sizeof(currentPath)))
          std::cout << "Warning: Can't obtain current directory";


        if(!files.interpretCMD(argc, argv, currentPath))
        {
          cout << "ERROR: Check command line arguments.\nProgram will now exit...\n";
          return 0;
        }


        /*
        //import the parameters file and set output file names
        files.importParams(argv[1],currentPath, argv[2]);
		if(argc == 4)
			files.setFileInfo(argv[2], argv[3], currentPath);
		else
			files.setFileInfo(argv[2], "", currentPath);

        files.printInfo();
*/

        //display program start with current time
        time_t now = time(0); // current date/time based on current system
        char* dt = ctime(&now); // convert now to string form
        double startCPUTime = omp_get_wtime();
        double endCPUTime = 0;
        //import the molecule file and data
        molecule2 mol;



        if(files.getMolFileExt().compare("pdb") == 0)
          mol.importPDB();
        else
          mol.importMFJ(files.getNameChanged(), files.getChargeChanged());
		
		if (mol.numAtoms == 0)
		{
			cout << "ERROR: Cannot read molecule file.\nProgram will now exit...\n";
			return 0;
		}


        files.printInfo();
	mol.printMolStats();


        //initialize the program to start running trajectories
        //see details in the constructor
        monte = new monteTraj(mol, shared);

        monte->runTraj();


        //exit the program with the time
        now = time(0); // current date/time based on current system
        dt = ctime(&now); // convert now to string form
        endCPUTime = omp_get_wtime();
        cout << " Ending Program " << endl;
         cout << " Current Time: " << dt << endl;
         cout << " Total Time: " << endCPUTime - startCPUTime << endl;

		 fileStream_ext.open(omegaFilePath_ext.c_str(), std::ios_base::app);
		 fileStream_ext << " Ending Program " << endl;
		 fileStream_ext << " Current Time: " << dt << endl;
		 fileStream_ext << " Total Time: " << endCPUTime - startCPUTime << endl;
		 fileStream_ext.close();

         monte->~monteTraj();

                return 0;
}
