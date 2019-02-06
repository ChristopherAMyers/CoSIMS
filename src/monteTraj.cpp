/* 
* Copyright (C) 2019  Christopher A. Myers
* 
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
* 
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
* 
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>
*
* If you have found this code usefull, please cite the research paper
* associated with this package.
*/

#include "monteTraj.h"

 void *callFunction(void *args[]);

 unsigned int microS = 100000;
monteTraj::monteTraj(molecule2 &mol, threadSharedParams &share)
{
  char buffer[500];
  char printOut[1000] = "";

  heliumMass = 4.002602;
  mass_ext = mu = mol.totalMass*heliumMass / (mol.totalMass + heliumMass);

  //initialize crossSection and number of trajectories
   omega = new double[iterations_ext];
   omegaErr = new double[iterations_ext];
   iterations = iterations_ext;
   omegaAvg = 0; omegaStdDev = 0;
   weightedMean = 0; weightedErr = 0;
   maxErgyPctDev = 0;

   //seed the RNG with the current time of day
   /*
   if(seed_ext == 0)
     randGen.seed(1234567);
   else
     randGen.seed((int)seed_ext);
     */

   //point to the molecule data and shared memory between threads
   molPtr = &mol;
   shared = &share;


   //initialize system variables
   temp = temp_ext; //degrees kelvin
   beta = mu/(K*temp);


   distrPower = 5;
   velDistr.setProperties(distrPower, beta/2.0);
   gMin = velDistr.mean() - 2*velDistr.stdDev();
   gMax = velDistr.mean() + 4*velDistr.stdDev();
   velDistr.renormalize(gMin, gMax);
   maxVelDistr = velDistr(velDistr.mode());


   //initialize the maxBoltz distibution and velocity components of the gas atom
   maxB.setSpecials(mu, temp);
   //gMin = maxB.getMean() - 2*maxB.getStdDev();
   //gMax = maxB.getMean() + 4*maxB.getStdDev();
   //maxVelDistr = maxB.getMaxDist();
   missCount = 0;
   
   //determine impact parameter ellipsoid
   fileStream_ext.open(omegaFilePath_ext.c_str(), ::ios_base::app);
   sprintf(buffer, "\n Creating potential energy surface...");
   cout << buffer; fileStream_ext << buffer; fileStream_ext.close();
   impactB.setMol(mol);
   impactB.setBoundry(0.5*mass_ext*pow(velDistr.mean(2, beta/2.0) - 2*velDistr.stdDev(2, beta/2.0), 2));
   //impactB.setBoundry(0.5*mass_ext*gMin*gMin);

   sprintf(buffer, "done\n Forming ellipsoid to surface...");
   cout << buffer; fileStream_ext << buffer; fileStream_ext.close();
   impactB.findEllipsoidBound();

   sprintf(buffer, "done\n");
   cout << buffer; fileStream_ext << buffer; fileStream_ext.close();
   impactB.printBoundryResults();
   

   //holds the farthest distance form an atom of the molecule
   maxDist = 0;
   
   //set the atom the farthest distance away from origin
   findMaxDist();
   bMax = impactB.getMax();
   maxDist = bMax;
   
   gasAtom.setMol(mol);
   //gasAtom.setMaxPos(M_SQRT2*bMax);


   //set constants used for integration
   const1 = pow(beta, 3)/(64*M_PI);
   const1 = const1*8*M_PI*M_PI;
   const1 = const1/velDistr.normFactor();

   //const1 = pow((mu*M_PI_2)/(K*temp), 1.5)/4;
   //const1 = pow(beta, 1.5)*sqrt(M_PI_2)*M_PI/8;
   const2 = bMax*bMax;
   maxB.setNormalization(gMin, gMax);
   /*
   cout << "Normalization Constant: " << (velDistr.cdf(gMax) - velDistr.cdf(gMin)) << endl;
   cout << "mu:   " << mu << endl;
   cout << "gMax: " << gMax << endl;
   cout << "bMax: " << bMax << endl;
   cout << "gMin: " << gMin << endl;
   cout << "gMin2: " << velDistr.mean(2, beta / 2.0) - 2 * velDistr.stdDev(2, beta / 2.0) << endl;
   */
   //const2 = const2*maxB.normConst;


   sprintf(buffer, " Velocity distribution statistics: \n");
   strcat(printOut, buffer);
   sprintf(buffer, "\t Maximum velocity:  %6.3f ang/ps\n", gMax);
   strcat(printOut, buffer);
   sprintf(buffer, "\t Minimum velocity:  %6.3f ang/ps\n", gMin);
   strcat(printOut, buffer);
   sprintf(buffer, "\t Mean velocity:     %6.3f ang/ps\n", velDistr.mean());
   strcat(printOut, buffer);
   sprintf(buffer, "\t Mode velocity:     %6.3f ang/ps\n", velDistr.mode());
   strcat(printOut, buffer);
   sprintf(buffer, "\t Percent of velocity \n\t distribution used: %6.3f%%\n\n", 100*(velDistr.cdf(gMax) - velDistr.cdf(gMin)));
   strcat(printOut, buffer);


   //print the rotated molecule to file;
   if(printMol_ext)
   {
     outMolFile.open(printMolName_ext.c_str(), std::ios_base::out);
     printMolecule();
     outMolFile.close();
   }


   if (printData_ext)
   {
	   data.open(dataFilePath_ext.c_str(), ios::binary);
	   sprintf(buffer, "%8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s \n"
	               , "phi", "theta", "gamma", "g", "b","chi", "steps", "Px", "Py", "Pz", "Vx", "Vy", "Vz");
	   data << buffer;
	   data.close();
   }


   cout << printOut << endl;
   fileStream_ext.open(omegaFilePath_ext.c_str(), std::ios_base::app);
   fileStream_ext << printOut << endl;
   fileStream_ext.close();


   //calculate the number of random numbers to skip per thread.
   //this is to the random number generator for each thread
   //always has at a non-overlapping sequence of numbers between threads.
   areaRatio = 1.0/((M_PI*0.5)*(gMax - gMin)*velDistr(velDistr.mode())*2);
   randNumSkipSize = (int)max((int)ceil(2*1*(8*trajNum_ext*iterations)/(areaRatio*numThreads_ext)), 100000);
   restartGen = new unsigned int[numThreads_ext];


}

monteTraj::~monteTraj() {

  energyFile.close();
  xyzFile.close();
}

void monteTraj::findMaxDist()
{
        //This function finds the atom of the molecule that is the farthest
        //distance from the origin. This distance will become the farthest
        //distance any of the gas atoms can travel away from the origin
        float farDist = 0; //distance from farthest atom
        float curDist; //current vector magnitude distance of each atom

        //loops through each atom of molecule,
        //finds the atom the farthest away form origin
        for (int i = 0; i < molPtr->numAtoms; i ++)
        {

                curDist = molPtr->data[i].coords.mag();
                if (curDist > farDist) farDist = curDist;
        }

        //sets the farthest distance to 150% of the distance of the atom
        //farDist = 0.6;
        maxDist = farDist;
}



bool monteTraj::setUpTraj(int atomNum, int &trajNum, vector3D<double> &initPos, vector3D<double> &initVel
    , localThreadParams &localParams, mt19937_64 &generator, uniform_real_distribution<double>& distr)
{
        //This function will choose randomly an orientation for the gas atoms
        //by the use of three euler angles. The gas atoms are initially positioned
        //at (maxDist,0,0) in cartesian coords, and then rotated by the randomly
        //chosen angles

        vector3D<double> tempPos(0.0,0.0,0.0); //temp vector to generate random coordinates
        vector3D<double> tempVel(0.0,0.0,0.0); //temp vector to generate velocity
        vector3D<double> p(0.0,0.0,0.0); //temp vector used to send particle to the surface of the ellipsoid
        vector3D<double> v(0.0,0.0,0.0); //temp vector used to send particle to the surface of the ellipsoid
        vector3D<double> crossVec(0.0,0.0,0.0); //temp vector used to send particle to the surface of the ellipsoid
        double positParam = 0; // temp variable used to send particle to the surface of the ellipsoid
        double phi, theta, gamma; //three Euler angles
        double b; //impact parameter
        double g; //gas atom velocity
        double thetaCheck; //value to check theta distribution against
        double bCheck; //value to check impact parameter distribution against
        double gCheck; //value to check velocity against max-boltz disto
        bool accept = false; //used to check if pairs of numbers aggree with monte-carlo distributions
        double discrim = 0;




//              new way of generating points
////////////////////////////////////////////////////////////////////////////////////////////////////////////

          accept = false;
          do
          {


            {
			 
              localParams.ranNumCount += 8;
              //assign random velocity magnitude and a
              //velocity to check against max-boltz distro
              g = distr(generator)*(gMax - gMin) + gMin;
              gCheck = ((double)rand() / (double)RAND_MAX)*maxVelDistr;


              //assign three random angles
              phi = distr(generator)*2.0*M_PI;
              theta = distr(generator)*M_PI;
              thetaCheck = 0.5*distr(generator);
              gamma = distr(generator)*2.0*M_PI;

              //assign a random impact parameter
              b = distr(generator)*bMax;
              bCheck = distr(generator)*(2 / bMax);

              
            }
            //check whether or not the the values fall within the distribution
            //unless all three check out, we will pick another set of points
            //if(thetaCheck <= sin(theta)/2 && bCheck <= 2*b/(bMax*bMax) && gCheck <= maxB.distro(g))
			if (thetaCheck <= sin(theta) / 2 && bCheck <= 2 * b / (bMax*bMax) && gCheck <= velDistr(g))
                accept = true;


          }while(!accept);


          //create vector for the starting points of the particle
          tempPos.X = b*cos(gamma)*cos(phi)*cos(theta) - b*sin(gamma)*sin(phi) + maxDist*cos(phi)*sin(theta);
          tempPos.Y = b*cos(gamma)*sin(phi)*cos(theta) + b*sin(gamma)*cos(phi) + maxDist*sin(phi)*sin(theta);
          tempPos.Z = maxDist*cos(theta) - b*cos(gamma)*sin(theta);

          //create the velocity vector for the incoming particle
          tempVel.X = -g*cos(phi)*sin(theta);
          tempVel.Y = -g*sin(phi)*sin(theta);
          tempVel.Z = -g*cos(theta);

          //save starting components
          if(printData_ext)
          {
            localParams.b     = b;
            localParams.phi   = phi;
            localParams.theta = theta;
            localParams.gamma = gamma;
          }
          localParams.g = g;

          if(projection_ext)
          {
            //send vector to the surface of the ellipsoid
            v.setEqual(tempVel.X/impactB.a, tempVel.Y/impactB.b, tempVel.Z/impactB.c);
            p.setEqual(tempPos.X/impactB.a, tempPos.Y/impactB.b, tempPos.Z/impactB.c);
            crossVec = p.cross(v);
            discrim = -crossVec.dot(crossVec) + v.dot(v);

            if(discrim < 0)
            {

              localParams.missCount++;
              //localParams.trajCount++;

              /*localParams.currOmega = localParams.total / localParams.trajCount;
              localParams.currError = sqrt((localParams.total2/localParams.trajCount
                  - (localParams.total/localParams.trajCount)*(localParams.total/localParams.trajCount))/localParams.trajCount);
              */
              return false;
            }
            positParam = -(p.dot(v)/(v.dot(v))) - sqrt(discrim)/v.dot(v);
            tempPos = tempPos + tempVel*positParam;
          }


          localParams.initPos = tempPos;
          localParams.initVel = tempVel;


          ///////////////////////////////

          initPos = tempPos;
          initVel = tempVel;

          return true;
}


void monteTraj::runTraj()
{
  double total, total2;
  int numFailed = 0;
  double tmpError = 0.0;
  char buffer[200];
  char printOut[2000] = "";
  double startTime = 0, endTime = 0;
  double totalTime = 0;

        sprintf(buffer, " Starting Trajectory Calculations:\n");
        strcat(printOut, buffer);
        sprintf(buffer, " _________________________________\n");
        strcat(printOut, buffer);
        cout << printOut << endl;
        fileStream_ext.open(omegaFilePath_ext.c_str(), std::ios_base::app);
        fileStream_ext << printOut << endl;
        fileStream_ext.close();
        strcpy(printOut, "");

        for(int traj = 0; traj < iterations; traj ++)
        {

          total = 0;
          total2 = 0;
          currTotal = currTotal2 = 0.0;
          currTrajCount = 0;
		  currSuccTraj = 0;
          startTime = omp_get_wtime();

          omp_set_num_threads(numThreads_ext);
          #pragma omp parallel
          {

            runningNumThreads = omp_get_num_threads();
            #pragma omp master
            {
                    sprintf(buffer, "\n Computing CCS integral %i / %i using %i OpenMP thread(s)\n", traj + 1, iterations, runningNumThreads);
                    cout << buffer;
                    fileStream_ext.open(omegaFilePath_ext.c_str(), std::ios_base::app);
                    fileStream_ext << buffer;
                    fileStream_ext.close();
            }

	    int id = omp_get_thread_num();
            bool succTraj;
            bool hitEllip;
            double chi;
            vector3D<double> initPos, initVel;
            localThreadParams p;
            succTraj = false;
            p.trajCount = 0; p.prevTrajCount = 0;
            p.hitCount = 0;
            p.lowChiCount = 0;
            p.orbitCount = 0;
            p.ranNumCount = 0;
            p.total = 0; p.prevTotal = 0;
            p.total2 = 0; p.prevTotal2 = 0;
			p.maxErgyPctDev = 0;
            updateCond = (int)ceil(((double)printRate_ext / (double)runningNumThreads));
            //loop through the trajectory

            //initialize random number generator state
            mt19937_64 generator(seed_ext);
            if(traj == 0)
              generator.discard(randNumSkipSize*id);
            else
              generator.discard(randNumSkipSize*id + restartGen[id]);
            uniform_real_distribution<double> distr(0.0, 1.0);


            #pragma omp for schedule(static)
            for(int n = 0; n < trajNum_ext; n ++)
            {
                 //setup new trajectory
                 //hitEllip = setUpTraj(0, n, initPos, initVel, p, localState);
                 hitEllip = setUpTraj(0, n, initPos, initVel, p, generator, distr);
                 if (hitEllip)
                 {
                     p.hitCount++;
                     //compute trajectory of gas atom
                     succTraj = false;
                     succTraj = gasAtom.verletVel(dt_ext
                             , initPos, initVel
                             , impactB.a, impactB.b, impactB.c, chi, p);

                     //dispResults(0, p);
                 }
                 else
                   chi = 0;
                 //cout << chi << "  " << p.stepCount << endl;
		runThreads(n, succTraj, chi, p);

            }

            updateCounts(p);
            #pragma omp critical
            {

                restartGen[id] = p.ranNumCount;
                total += p.total;
                total2 += p.total2;
                maxErgyPctDev = max(maxErgyPctDev, p.maxErgyPctDev);
            }
          }
		  numFailed += trajNum_ext - currSuccTraj;
          omega[traj] = total/ currSuccTraj;
          endTime = omp_get_wtime();
          totalTime += endTime - startTime;

          //cout << "Number of completed Trajectories: " << trajNum_ext << endl;
          //cout << "Omega: " << omega[traj] << " +- " << sqrt((total2/trajNum_ext
          //    - (total/trajNum_ext)*(total/trajNum_ext))/trajNum_ext) << endl;
          //cout << "Completed Iteration " << traj << endl << endl;

          tmpError = sqrt((total2/ currSuccTraj
              - (total/ currSuccTraj)*(total/ currSuccTraj))/ currSuccTraj);
		  omegaErr[traj] = tmpError;

          sprintf(buffer, " Number of completed trajectories: %i\n", trajNum_ext);
          strcat(printOut, buffer);

          sprintf(buffer, " Omega: %11.2f +/- %-*.2f (%*.2f %%)\n", omega[traj]
              , (int)log10(tmpError) + 3, tmpError, (int)log10(100*tmpError/omega[traj]) + 3
              , 100*tmpError/omega[traj]);
          strcat(printOut, buffer);

          sprintf(buffer, " Total time: %5.3f; estimated time remaining: %5.3f\n", totalTime, totalTime*(iterations - traj - 1)/(traj + 1));
          strcat(printOut, buffer);



          cout << printOut;
          fileStream_ext.open(omegaFilePath_ext.c_str(), std::ios_base::app);
          fileStream_ext << printOut;
          fileStream_ext.close();
          strcpy(printOut, "");

          //fileStream_ext.open(omegaFilePath_ext.c_str(), std::ios_base::app);
          //fileStream_ext << omega[traj] << "\t" << shared->tempError << "\t" << 100 * shared->tempError / shared->tempOmega << endl;
          //fileStream_ext.close();
        }


        //calculate and print average cross-section and standard deviation endl;
        for(int i = 0; i < iterations; i ++)
          omegaAvg += omega[i];

        omegaAvg = omegaAvg/(double)iterations;

		for (int i = 0; i < iterations; i++)
		{
			omegaStdDev += (omega[i] - omegaAvg)*(omega[i] - omegaAvg);
			weightedErr += 1.0 / (omegaErr[i] * omegaErr[i]);
			weightedMean += omega[i] / (omegaErr[i] * omegaErr[i]);
		}
		weightedMean /= weightedErr;
        omegaStdDev = sqrt(omegaStdDev/(double)iterations);
		weightedErr = 1.0 / sqrt(weightedErr);

		
		strcpy(printOut, "");
		sprintf(buffer, "\n\n ___________________________________\n");
		strcat(printOut, buffer);
		sprintf(buffer, "   Mean Cross-section:   %9.3f\n", omegaAvg);
		strcat(printOut, buffer);
		sprintf(buffer, "   Standard Dev:         %9.3f\n", omegaStdDev);
		strcat(printOut, buffer);
		sprintf(buffer, "   StdDev Percent:       %9.3f\n\n", (omegaStdDev * 100 / omegaAvg));
		strcat(printOut, buffer);
		sprintf(buffer, "   Weighted CCS Mean:    %9.3f\n", weightedMean);
		strcat(printOut, buffer);
		sprintf(buffer, "   Weighted Error:       %9.3f\n", weightedErr);
		strcat(printOut, buffer);
		sprintf(buffer, "   Weighted Pct. Error:  %9.4f\n", (weightedErr * 100 / weightedMean));
		strcat(printOut, buffer);
		sprintf(buffer, "   Max Energy Deviation: %9.4f\n", maxErgyPctDev);
		strcat(printOut, buffer);
		//sprintf(buffer, "   No. of failed traj:   %9i\n", numFailed);
		//strcat(printOut, buffer);
		sprintf(buffer, " ___________________________________\n");
		strcat(printOut, buffer);

		cout << printOut;
		fileStream_ext.open(omegaFilePath_ext.c_str(), ::ios_base::app);
		fileStream_ext << printOut;
		fileStream_ext.close();

}


//void *monteTraj::runThreads(void *args)
void monteTraj::runThreads(int i, bool succTraj, double chi, localThreadParams &p)
{
         p.trajCount ++;
         if(succTraj)
         {
           //p.hitCount ++;
         }
         else
         {
           p.orbitCount++;
         }

         if(chi < M_PI/100.0)
           p.lowChiCount ++;

         //double dOmega = (p.g*p.g*p.g)*(1 - cos(abs(chi)));
		 //double dOmega = (1 - cos(abs(chi)));
		 double dOmega = pow(p.g, 5 - distrPower)*(1 - cos(abs(chi)));
         double tmpTotal = const1*const2*dOmega;

         p.total += tmpTotal;
         p.total2 += tmpTotal*tmpTotal;
         p.currOmega = p.total/p.trajCount;
         p.currError = sqrt((p.total2/p.trajCount - (p.total/p.trajCount)*(p.total/p.trajCount))/p.trajCount);

         //printf("Total: %i %8.5f\n", i, chi);

         //update global counters
         if (p.trajCount % updateCond == 0)
         {
               updateCounts(p);
               #pragma omp critical
               {
                 if ((currTrajCount % printRate_ext) <= prevCountMod && displayProg_ext)
                   dispResults(p);
               }


               /*
                 //each thread updates it's own local counters to the global one
                #pragma omp critical
                {
                        prevCountMod = currTrajCount % printRate_ext;
                        //currTrajCount += updateCond;
                        currTrajCount += p.trajCount - p.prevTrajCount;
                        //currSuccTraj += updateCond - p.orbitCount;
                        currSuccTraj = currTrajCount - p.orbitCount;

                        currTotal += (p.total - p.prevTotal);
                        currTotal2 += (p.total2 - p.prevTotal2);

                        p.orbitCount = 0;
                        p.prevTotal = p.total;
                        p.prevTotal2 = p.total2;

                        if ((currTrajCount % printRate_ext) <= prevCountMod && displayProg_ext)
                                dispResults(p);
                }
                p.prevTrajCount = p.trajCount;
                */
         }

         //check if any printing to file is performed
         if(printTraj_ext)
         {
           if(succTraj != printTrajOnFail_ext)
             printTraj(true);
           else
             printTraj(false);
         }

         if(printData_ext)
           if(succTraj != printDataOnFail_ext)
             printData(i, chi, p);

}

void monteTraj::updateCounts(localThreadParams &p)
{
          //each thread updates it's own local counters to the global one
         #pragma omp critical
         {
                 prevCountMod = currTrajCount % printRate_ext;
                 //currTrajCount += updateCond;
                 currTrajCount += p.trajCount - p.prevTrajCount;
                 //currSuccTraj += updateCond - p.orbitCount;
                 currSuccTraj = currTrajCount - p.orbitCount;

                 currTotal += (p.total - p.prevTotal);
                 currTotal2 += (p.total2 - p.prevTotal2);
         }
         p.orbitCount = 0;
         p.prevTotal = p.total;
         p.prevTotal2 = p.total2;
         p.prevTrajCount = p.trajCount;

}

void monteTraj::dispResults(localThreadParams &p)
{
	//double error = sqrt((currTotal2*currTrajCount - currTotal*currTotal)
		///(currTrajCount*(currTrajCount - 1)));

        fileStream_ext.open(omegaFilePath_ext.c_str(), std::ios_base::app);
        char buffer[1000];

		double numSucc = currSuccTraj;
        double error = sqrt((currTotal2/ numSucc
            - (currTotal/ numSucc)*(currTotal/ numSucc))/ numSucc);
	double errorPct = error * 100 * (double)numSucc / currTotal;


	sprintf(buffer, " progress: Traj %*i of %i (%4.1f %%); CCS = %#9.5g +- %9.5g (%4.2f %%)\n"
		, (int)log10(trajNum_ext), currTrajCount, trajNum_ext, double(100 * currTrajCount) / double(trajNum_ext)
		, (currTotal / (double)numSucc), error, errorPct);

	cout << buffer;
	fileStream_ext << buffer;
	fileStream_ext.close();

}


void monteTraj::printTraj(bool toFile)
{


  if(toFile)
  {
    cout << "Printing file: " << xyzStream.str() << endl;
    xyzFile << fixed << showpoint << setprecision(15) << xyzStream.str();
    //cout << fixed << showpoint << setprecision(15) << xyzStream.str();
  }
    xyzStream.str("");
}

void monteTraj::printData(int i, double chi, localThreadParams &p)
{

  data.open(dataFilePath_ext.c_str(),ios_base::app);

        char buffer[1000];
		sprintf(buffer, "%8.3g  %8.3g  %8.3g  %8.3g  %8.3g  %8.3g  %8i  %8.3g  %8.3g  %8.3g  %8.3g  %8.3g  %8.3g \n"
			, p.phi
			, p.theta, p.gamma
			, p.g, p.b
			, chi, p.stepCount, p.initPos.X
			, p.initPos.Y, p.initPos.Z
			, p.initVel.X, p.initVel.Y
			, p.initVel.Z);

        data << buffer;
        data.close();
}

void monteTraj::storeTraj()
{
	
        xyzStream << (numGas) << "\n";
        xyzStream << "no comments available \n";
        vector3D<double> tempVec(0,0,0);

		xyzStream << setprecision(15);

          xyzStream << fixed << showpoint << setprecision(15) << "He     " << gasAtom.getPosition().X << "       "
                << gasAtom.getPosition().Y << "       " << gasAtom.getPosition().Z;
        xyzStream << fixed << showpoint << setprecision(15) << "     " << gasAtom.getVelocity().X << "       "
                        << gasAtom.getVelocity().Y << "       " << gasAtom.getVelocity().Z <<"\n";

}

void monteTraj::printMolecule()
{
  outMolFile << (1*molPtr->numAtoms) << "\n";
  outMolFile << "Rotated Molecule \n";
  vector3D<double> tempVec(0,0,0);
  char buffer[1000];

  for (int i = 0; i < molPtr->numAtoms; i ++)
  {

    sprintf(buffer, "%s %9.4g %9.4g %9.4g\n", molPtr->data[i].atom.c_str(), molPtr->data[i].coords.X
        , molPtr->data[i].coords.Y, molPtr->data[i].coords.Z);
    outMolFile << buffer;
  }
}
