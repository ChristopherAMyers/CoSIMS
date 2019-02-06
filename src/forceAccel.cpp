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

#include "forceAccel.h"

namespace std
{

  forceAccel::forceAccel()
  {
    mass = 0;

    cutOffRadius = 12.0;
    resolution = 1.0;

    heliumMass = 4.002602;

    dipoleScalar = 28276.01363 / mass_ext;

    forceType = 6;
    mol = NULL;
	totalCenterReplaced = 0;

	avgClusterRad = 0;
  }

  forceAccel::~forceAccel()
  {
    // TODO Auto-generated destructor stub
  }

  void forceAccel::setMol(molecule2 &inMol)
  {
	  mol = &inMol;
  }

  void forceAccel::setClusters()
  {
	  char buffer[200];
	  char printOut[2000] = "";

	  for (int i = 0; i < mol->numAtoms; i++)
	  {
		  cluster.addCoord(mol->data[i].coords);
	  }

	  if (multipole_ext || dispersionCutoff_ext)
	  {



		  sprintf(buffer, " Cutt-off scheme requested via input variabl");
		  if (multipole_ext && !dispersionCutoff_ext)
			  strcat(buffer, "e MULTIPOLE: \n");
		  else if (!multipole_ext && dispersionCutoff_ext)
			  strcat(buffer, "e DISPERSION_CUTOFF: \n");
		  else
			  strcat(buffer, "es MULTIPOLE and DISPERSION_CUTOFF: \n");

		  cluster.formClusters();

		  for (int i = 0; i < (signed)cluster.getNumClusters(); i++)
		  {
			  clusterGroups.push_back(cluster.getCluster(i));
			  addMultipole(cluster.getCluster(i), cluster.getCenter(i));
			  //clusterPos.push_back(cluster.getCenter(i));
			  
			  clusterRadius.push_back(cluster.getMaxDist(i));
			  addClusterPos(cluster.getCenter(i), i);

			  //cout << "ClusterRad: " << clusterRadius[i] << "  " << avgClusterRad << endl;
			  avgClusterRad += clusterRadius[i];
			  avgClusterSize += clusterGroups[i].size();
		  }

		  //printf(buffer);
		  //fileStream_ext << buffer;
		  //fileStream_ext.close();

		  avgClusterSize /= (double)clusterGroups.size();
		  avgClusterRad /= (double)clusterGroups.size();

		  sprintf(buffer, " Cluster formation summary:\n");
		  strcat(printOut, buffer);

		  if (multipole_ext)
		  {
			  sprintf(buffer, "\t Number of centers replaced\n\t by center of charge:          %6i\n", totalCenterReplaced);
			  strcat(printOut, buffer);
		  }

		  sprintf(buffer, "\t Number of clusters formed:    %6i\n", (int)clusterGroups.size());
		  strcat(printOut, buffer);
		  sprintf(buffer, "\t Average cluster size:         %6.2f atoms\n", avgClusterSize);
		  strcat(printOut, buffer);
		  sprintf(buffer, "\t Average cluster radius:       %6.2f Ang.\n", avgClusterRad);
		  strcat(printOut, buffer);

		  if (multipole_ext)
		  {
			  sprintf(buffer, "\t Multipole cut-off radius:     %6.1f Ang.\n", multipole_radius_ext);
			  strcat(printOut, buffer);
		  }
		  if (dispersionCutoff_ext)
		  {
			  sprintf(buffer, "\t Dispersion cut-off radius:    %6.1f Ang.\n", dispersion_radius_ext);
			  strcat(printOut, buffer);
		  }

		  fileStream_ext.open(omegaFilePath_ext.c_str(), std::ios_base::app);
		  cout << printOut << endl;
		  fileStream_ext << printOut << endl;
		  fileStream_ext.close();
	  }


	  if (charge_ext)
		  if (multipole_ext)
			  if (dispersionCutoff_ext)
				  forceType = 1;
			  else
				  forceType = 2;
		  else
			  if (dispersionCutoff_ext)
				  forceType = 3;
			  else
				  forceType = 4;
	  else
		  if (dispersionCutoff_ext)
			  forceType = 5;
		  else
			  forceType = 6;

	  //printClusters();

  }

  void forceAccel::addMultipole(vector<int> cluster, vector3D<double> center)
  {
	  double charge = 0.0;
	  int index = 0;
	  double sumCharge = 0.0;
	  double pos[3] = { 0, 0, 0 };
	  double mag2 = 0;

	  const double arr[] = { 0, 0, 0 };
	  vector<double> sumDipole(arr, arr + sizeof(arr) / sizeof(arr[0]));
	  vector<vector<double> > sumQuadPol;
	  sumQuadPol.push_back(sumDipole);
	  sumQuadPol.push_back(sumDipole);
	  sumQuadPol.push_back(sumDipole);


	  for (int n = 0; n < cluster.size(); n++)
	  {
		  index = cluster[n];
		  charge = mol->data[index].charge;
		  sumCharge += charge;

		  pos[0] = mol->data[index].coords.at(0) - center.X;
		  pos[1] = mol->data[index].coords.at(1) - center.Y;
		  pos[2] = mol->data[index].coords.at(2) - center.Z;
		  mag2 = pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2];


		  for (int i = 0; i < 3; i++)
		  {
			  sumDipole[i] += charge * pos[i];

			  for (int j = 0; j < 3; j++)
			  {
				  sumQuadPol[i][j] += charge*(3 * pos[i] * pos[j] - mag2*delta(i, j));
			  }
		  }
	  }

	  totalCharge.push_back(sumCharge);
	  dipole.push_back(sumDipole);
	  quadPole.push_back(sumQuadPol);
  }

  void forceAccel::addClusterPos(vector3D<double> center, int index)
  {
	  vector3D<double> atomPos, centerOfCharge(0, 0, 0);
	  double percentDiff, charge;
	  for (int n = 0; n < (signed)clusterGroups[index].size(); n++)
	  {
		  atomPos = mol->data[clusterGroups[index][n]].coords;
		  charge = mol->data[clusterGroups[index][n]].charge;
		  centerOfCharge = centerOfCharge + atomPos*charge;
	  }
	  centerOfCharge = centerOfCharge / totalCharge[index];
	  percentDiff = 100*abs((centerOfCharge.mag() - center.mag()) / centerOfCharge.mag());
	  if (percentDiff <= 5)
	  {
		  centerReplace.push_back(true);
		  clusterPos.push_back(centerOfCharge);
		  totalCenterReplaced++;
	  }
	  else
	  {
		  centerReplace.push_back(false);
		  clusterPos.push_back(center);
	  }
  }

  double forceAccel::delta(int a, int b)
  {
	  if (a == b)
		  return 1.0;
	  else
		  return 0.0;
  }

 
  vector3D<double> forceAccel::multipoleCutoff(vector3D<double> &position)
  {
      double fieldX, fieldY, fieldZ;
      double sMatXY, sMatXZ, sMatYZ, aMatXY, aMatXZ, aMatYZ;
      double sumX, sumY, sumZ;
      double diagMatXX, diagMatYY, diagMatZZ;
      double xx, yy, zz, xy, xz, yz;
      double dx[3] = { 0.0, 0.0, 0.0 };
      double chargeR_5, chargeR_3;
      double const1, const2;
      double tMat2[3];
      double term1, term2, scalarComp;
      double r2, r_1, r_2, r_3, r_5, r_7, r_9, r6, charge;
      int index;
      int i, k;
      double multiPlusRad, dispPlusRad;

      fieldX = fieldY = fieldZ = 0;
      sMatXY = sMatXZ = sMatYZ = aMatXY = aMatXZ = aMatYZ = 0;
      diagMatXX = diagMatYY = diagMatZZ = 0;
      sumX = sumY = sumZ = 0;


      for(int n = 0; n < (signed)clusterGroups.size(); n++)
      {

          dx[0] = position.X - clusterPos[n].X;
          dx[1] = position.Y - clusterPos[n].Y;
          dx[2] = position.Z - clusterPos[n].Z;
          xx = dx[0]*dx[0]; yy = dx[1]*dx[1]; zz = dx[2]*dx[2];
          r2 = xx + yy + zz;

          multiPlusRad = multipole_radius_ext + clusterRadius[n];
          multiPlusRad *= multiPlusRad;
          dispPlusRad = dispersion_radius_ext + clusterRadius[n];
		  dispPlusRad *= dispPlusRad;
		  //cout << n << "  " << dispPlusRad << "  " << r2 << "  " << clusterRadius[n] << "  " << dispersion_radius_ext << endl;

          if (r2 > multiPlusRad)
          {

              r_2 = 1.0 / r2;
              r_1 = 1.0 / sqrt(r2);
              r_3 = r_2*r_1;
              r_5 = r_3*r_2;

              xy = dx[0]*dx[1]; xz = dx[0]*dx[2]; yz = dx[1]*dx[2];
              chargeR_5 = totalCharge[n]*r_5;
              chargeR_3 = chargeR_5*r2;

              if(multipole_order_ext >= 1 && !centerReplace[n])
              {
                r_7 = r_5*r_2;
                const1 = 3*(dipole[n][0]*dx[0] + dipole[n][1]*dx[1] + dipole[n][2]*dx[2]);
              }

              if (multipole_order_ext >= 2)
              {
                r_7 = r_5*r_2;
                r_9 = r_5*r_2*r_2;
                const2 = 0;
                tMat2[0] = (quadPole[n][0][0]*dx[0] + quadPole[n][1][0]*dx[1] + quadPole[n][2][0]*dx[2])*r_7;
                tMat2[1] = (quadPole[n][0][1]*dx[0] + quadPole[n][1][1]*dx[1] + quadPole[n][2][1]*dx[2])*r_7;
                tMat2[2] = (quadPole[n][0][2]*dx[0] + quadPole[n][1][2]*dx[1] + quadPole[n][2][2]*dx[2])*r_7;
                const2 = 5*(tMat2[0]*dx[0] + tMat2[1]*dx[1] + tMat2[2]*dx[2])*r_2*0.5;
              }

              fieldX += chargeR_3*dx[0];
              fieldY += chargeR_3*dx[1];
              fieldZ += chargeR_3*dx[2];

              diagMatXX += (r2 - 3*xx)*chargeR_5;
              diagMatYY += (r2 - 3*yy)*chargeR_5;
              diagMatZZ += (r2 - 3*zz)*chargeR_5;

              sMatXY -= 3*xy*chargeR_5;
              sMatXZ -= 3*xz*chargeR_5;
              sMatYZ -= 3*yz*chargeR_5;

              if(multipole_order_ext >= 1 && !centerReplace[n])
              {
                fieldX += (const1*dx[0] - r2*dipole[n][0])*r_5;
                fieldY += (const1*dx[1] - r2*dipole[n][1])*r_5;
                fieldZ += (const1*dx[2] - r2*dipole[n][2])*r_5;

                diagMatXX += (r2 - 5*xx)*r_7*const1;
                diagMatYY += (r2 - 5*yy)*r_7*const1;
                diagMatZZ += (r2 - 5*zz)*r_7*const1;

                sMatXY -= 5*xy*r_7*const1;
                sMatXZ -= 5*xz*r_7*const1;
                sMatYZ -= 5*yz*r_7*const1;

                aMatXY += 3*(dipole[n][1]*dx[0] - dipole[n][0]*dx[1])*r_5;
                aMatXZ += 3*(dipole[n][2]*dx[0] - dipole[n][0]*dx[2])*r_5;
                aMatYZ += 3*(dipole[n][2]*dx[1] - dipole[n][1]*dx[2])*r_5;
              }

              if(multipole_order_ext >= 2)
              {

                fieldX += (const2*dx[0] - tMat2[0])*r2;
                fieldY += (const2*dx[1] - tMat2[1])*r2;
                fieldZ += (const2*dx[2] - tMat2[2])*r2;

                diagMatXX += ((r2 - 7*xx)*const2 - quadPole[n][0][0]*r_5);
                diagMatYY += ((r2 - 7*yy)*const2 - quadPole[n][1][1]*r_5);
                diagMatZZ += ((r2 - 7*zz)*const2 - quadPole[n][2][2]*r_5);

                sMatXY -= (7*xy*const2 + quadPole[n][0][1]*r_5);
                sMatXZ -= (7*xz*const2 + quadPole[n][0][2]*r_5);
                sMatYZ -= (7*yz*const2 + quadPole[n][1][2]*r_5);

                aMatXY += 5*(tMat2[1]*dx[0] - tMat2[0]*dx[1]);
                aMatXZ += 5*(tMat2[2]*dx[0] - tMat2[0]*dx[2]);
                aMatYZ += 5*(tMat2[2]*dx[1] - tMat2[1]*dx[2]);

              }

              if (!(r2 > dispPlusRad))
              {
                for (k = 0; k < (signed)clusterGroups[n].size(); k++)
                 {
                         index = clusterGroups[n][k];

                         charge = mol->data[index].charge;
                         dx[0] = position.X - mol->data[index].coords.X;
                         dx[1] = position.Y - mol->data[index].coords.Y;
                         dx[2] = position.Z - mol->data[index].coords.Z;

                         r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

                         r6 = r2*r2*r2;
                         term2 = mol->data[n].sigAvg6 / r6;
                         term1 = 2 * term2*term2;
                         scalarComp = mol->data[n].eps24_mass*(term1 - term2) / r2;

                         sumX += scalarComp*dx[0];
                         sumY += scalarComp*dx[1];
                         sumZ += scalarComp*dx[2];
                 }
              }
          }

          else if(!(r2 > dispPlusRad))
          {

              for (k = 0; k < (signed)clusterGroups[n].size(); k++)
              {
                      index = clusterGroups[n][k];

                      charge = mol->data[index].charge;
                      dx[0] = position.X - mol->data[index].coords.X;
                      dx[1] = position.Y - mol->data[index].coords.Y;
                      dx[2] = position.Z - mol->data[index].coords.Z;

                      r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

                      r_5 = 1.0 / (r2*r2*sqrt(r2));
                      r_3 = r_5*r2;

                      r_3 *= charge;
                      r_5 *= charge;

                      fieldX += dx[0]*r_3;
                      fieldY += dx[1]*r_3;
                      fieldZ += dx[2]*r_3;

                      diagMatXX += (r2 - 3 * dx[0]*dx[0])*r_5;
                      diagMatYY += (r2 - 3 * dx[1]*dx[1])*r_5;
                      diagMatZZ += (r2 - 3 * dx[2]*dx[2])*r_5;

                      sMatXY -= 3 * dx[0]*dx[1]*r_5;
                      sMatXZ -= 3 * dx[0]*dx[2]*r_5;
                      sMatYZ -= 3 * dx[1]*dx[2]*r_5;

                      r6 = r2*r2*r2;
                      term2 = mol->data[index].sigAvg6 / r6;
                      term1 = 2 * term2*term2;
                      scalarComp = mol->data[index].eps24_mass*(term1 - term2) / r2;

                      sumX += scalarComp*dx[0];
                      sumY += scalarComp*dx[1];
                      sumZ += scalarComp*dx[2];
              }
          }

          else
          {
            for (k = 0; k < (signed)clusterGroups[n].size(); k++)
             {
                      index = clusterGroups[n][k];

                      charge = mol->data[index].charge;
                      dx[0] = position.X - mol->data[index].coords.X;
                      dx[1] = position.Y - mol->data[index].coords.Y;
                      dx[2] = position.Z - mol->data[index].coords.Z;

                      r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

                      r_5 = 1.0 / (r2*r2*sqrt(r2));
                      r_3 = r_5*r2;

                      r_3 *= charge;
                      r_5 *= charge;

                      fieldX += dx[0]*r_3;
                      fieldY += dx[1]*r_3;
                      fieldZ += dx[2]*r_3;

                      diagMatXX += (r2 - 3 * dx[0]*dx[0])*r_5;
                      diagMatYY += (r2 - 3 * dx[1]*dx[1])*r_5;
                      diagMatZZ += (r2 - 3 * dx[2]*dx[2])*r_5;

                      sMatXY -= 3 * dx[0]*dx[1]*r_5;
                      sMatXZ -= 3 * dx[0]*dx[2]*r_5;
                      sMatYZ -= 3 * dx[1]*dx[2]*r_5;
             }
          }
      }

      sumX +=         dipoleScalar * (diagMatXX*fieldX + (sMatXY + aMatXY)*fieldY + (sMatXZ + aMatXZ)*fieldZ);
      sumY += dipoleScalar * ((sMatXY - aMatXY)*fieldX +         diagMatYY*fieldY + (sMatYZ + aMatYZ)*fieldZ);
      sumZ += dipoleScalar * ((sMatXZ - aMatXZ)*fieldX + (sMatYZ - aMatYZ)*fieldY +         diagMatZZ*fieldZ);

      return vector3D<double>(sumX, sumY, sumZ);
  }

  vector3D<double> forceAccel::multipoleOnly(vector3D<double> &position)
  {
      double fieldX, fieldY, fieldZ;
      double sMatXY, sMatXZ, sMatYZ, aMatXY, aMatXZ, aMatYZ;
      double sumX, sumY, sumZ;
      double diagMatXX, diagMatYY, diagMatZZ;
      double xx, yy, zz, xy, xz, yz;
      double dx[3] = { 0.0, 0.0, 0.0 };
      double chargeR_5, chargeR_3;
      double const1, const2;
      double tMat2[3];
      double term1, term2, scalarComp;
      double r2, r_1, r_2, r_3, r_5, r_7, r_9, r6, charge;
      int index;
      int i, k;
      double multiPlusRad, dispPlusRad;

      fieldX = fieldY = fieldZ = 0;
      sMatXY = sMatXZ = sMatYZ = aMatXY = aMatXZ = aMatYZ = 0;
      diagMatXX = diagMatYY = diagMatZZ = 0;
      sumX = sumY = sumZ = 0;


      for(int n = 0; n < (signed)clusterGroups.size(); n++)
      {

          dx[0] = position.X - clusterPos[n].X;
          dx[1] = position.Y - clusterPos[n].Y;
          dx[2] = position.Z - clusterPos[n].Z;
          xx = dx[0]*dx[0]; yy = dx[1]*dx[1]; zz = dx[2]*dx[2];
          r2 = xx + yy + zz;

          multiPlusRad = multipole_radius_ext + clusterRadius[n];
          multiPlusRad *= multiPlusRad;
          dispPlusRad = dispersion_radius_ext + clusterRadius[n];
                  dispPlusRad *= dispPlusRad;
                  //cout << n << "  " << dispPlusRad << "  " << r2 << "  " << clusterRadius[n] << "  " << dispersion_radius_ext << endl;

          if (r2 > multiPlusRad)
          {

              r_2 = 1.0 / r2;
              r_1 = 1.0 / sqrt(r2);
              r_3 = r_2*r_1;
              r_5 = r_3*r_2;

              xy = dx[0]*dx[1]; xz = dx[0]*dx[2]; yz = dx[1]*dx[2];
              chargeR_5 = totalCharge[n]*r_5;
              chargeR_3 = chargeR_5*r2;

              if(multipole_order_ext >= 1 && !centerReplace[n])
              {

                r_7 = r_5*r_2;
                const1 = 3*(dipole[n][0]*dx[0] + dipole[n][1]*dx[1] + dipole[n][2]*dx[2]);
              }

              if (multipole_order_ext >= 2)
              {

                r_7 = r_5*r_2;
                r_9 = r_5*r_2*r_2;
                const2 = 0;
                tMat2[0] = (quadPole[n][0][0]*dx[0] + quadPole[n][1][0]*dx[1] + quadPole[n][2][0]*dx[2])*r_7;
                tMat2[1] = (quadPole[n][0][1]*dx[0] + quadPole[n][1][1]*dx[1] + quadPole[n][2][1]*dx[2])*r_7;
                tMat2[2] = (quadPole[n][0][2]*dx[0] + quadPole[n][1][2]*dx[1] + quadPole[n][2][2]*dx[2])*r_7;
                const2 = 5*(tMat2[0]*dx[0] + tMat2[1]*dx[1] + tMat2[2]*dx[2])*r_2*0.5;
              }
              fieldX += chargeR_3*dx[0];
              fieldY += chargeR_3*dx[1];
              fieldZ += chargeR_3*dx[2];

              diagMatXX += (r2 - 3*xx)*chargeR_5;
              diagMatYY += (r2 - 3*yy)*chargeR_5;
              diagMatZZ += (r2 - 3*zz)*chargeR_5;

              sMatXY -= 3*xy*chargeR_5;
              sMatXZ -= 3*xz*chargeR_5;
              sMatYZ -= 3*yz*chargeR_5;

              if(multipole_order_ext >= 1 && !centerReplace[n])
              {
                fieldX += (const1*dx[0] - r2*dipole[n][0])*r_5;
                fieldY += (const1*dx[1] - r2*dipole[n][1])*r_5;
                fieldZ += (const1*dx[2] - r2*dipole[n][2])*r_5;

                diagMatXX += (r2 - 5*xx)*r_7*const1;
                diagMatYY += (r2 - 5*yy)*r_7*const1;
                diagMatZZ += (r2 - 5*zz)*r_7*const1;

                sMatXY -= 5*xy*r_7*const1;
                sMatXZ -= 5*xz*r_7*const1;
                sMatYZ -= 5*yz*r_7*const1;

                aMatXY += 3*(dipole[n][1]*dx[0] - dipole[n][0]*dx[1])*r_5;
                aMatXZ += 3*(dipole[n][2]*dx[0] - dipole[n][0]*dx[2])*r_5;
                aMatYZ += 3*(dipole[n][2]*dx[1] - dipole[n][1]*dx[2])*r_5;
              }

              if(multipole_order_ext >= 2)
              {

                fieldX += (const2*dx[0] - tMat2[0])*r2;
                fieldY += (const2*dx[1] - tMat2[1])*r2;
                fieldZ += (const2*dx[2] - tMat2[2])*r2;

                diagMatXX += ((r2 - 7*xx)*const2 - quadPole[n][0][0]*r_5);
                diagMatYY += ((r2 - 7*yy)*const2 - quadPole[n][1][1]*r_5);
                diagMatZZ += ((r2 - 7*zz)*const2 - quadPole[n][2][2]*r_5);

                sMatXY -= (7*xy*const2 + quadPole[n][0][1]*r_5);
                sMatXZ -= (7*xz*const2 + quadPole[n][0][2]*r_5);
                sMatYZ -= (7*yz*const2 + quadPole[n][1][2]*r_5);

                aMatXY += 5*(tMat2[1]*dx[0] - tMat2[0]*dx[1]);
                aMatXZ += 5*(tMat2[2]*dx[0] - tMat2[0]*dx[2]);
                aMatYZ += 5*(tMat2[2]*dx[1] - tMat2[1]*dx[2]);

              }

              for (k = 0; k < (signed)clusterGroups[n].size(); k++)
               {
                       index = clusterGroups[n][k];

                       charge = mol->data[index].charge;
                       dx[0] = position.X - mol->data[index].coords.X;
                       dx[1] = position.Y - mol->data[index].coords.Y;
                       dx[2] = position.Z - mol->data[index].coords.Z;

                       r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

                       r6 = r2*r2*r2;
                       term2 = mol->data[index].sigAvg6 / r6;
                       term1 = 2 * term2*term2;
                       scalarComp = mol->data[index].eps24_mass*(term1 - term2) / r2;

                       sumX += scalarComp*dx[0];
                       sumY += scalarComp*dx[1];
                       sumZ += scalarComp*dx[2];
               }
          }

          else
          {

              for (k = 0; k < (signed)clusterGroups[n].size(); k++)
              {
                      index = clusterGroups[n][k];

                      charge = mol->data[index].charge;
                      dx[0] = position.X - mol->data[index].coords.X;
                      dx[1] = position.Y - mol->data[index].coords.Y;
                      dx[2] = position.Z - mol->data[index].coords.Z;

                      r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

                      r_5 = 1.0 / (r2*r2*sqrt(r2));
                      r_3 = r_5*r2;

                      r_3 *= charge;
                      r_5 *= charge;

                      fieldX += dx[0]*r_3;
                      fieldY += dx[1]*r_3;
                      fieldZ += dx[2]*r_3;

                      diagMatXX += (r2 - 3 * dx[0]*dx[0])*r_5;
                      diagMatYY += (r2 - 3 * dx[1]*dx[1])*r_5;
                      diagMatZZ += (r2 - 3 * dx[2]*dx[2])*r_5;

                      sMatXY -= 3 * dx[0]*dx[1]*r_5;
                      sMatXZ -= 3 * dx[0]*dx[2]*r_5;
                      sMatYZ -= 3 * dx[1]*dx[2]*r_5;

                      r6 = r2*r2*r2;
                      term2 = mol->data[index].sigAvg6 / r6;
                      term1 = 2 * term2*term2;
                      scalarComp = mol->data[index].eps24_mass*(term1 - term2) / r2;

                      sumX += scalarComp*dx[0];
                      sumY += scalarComp*dx[1];
                      sumZ += scalarComp*dx[2];
              }
          }
      }

      sumX +=         dipoleScalar * (diagMatXX*fieldX + (sMatXY + aMatXY)*fieldY + (sMatXZ + aMatXZ)*fieldZ);
      sumY += dipoleScalar * ((sMatXY - aMatXY)*fieldX +         diagMatYY*fieldY + (sMatYZ + aMatYZ)*fieldZ);
      sumZ += dipoleScalar * ((sMatXZ - aMatXZ)*fieldX + (sMatYZ - aMatYZ)*fieldY +         diagMatZZ*fieldZ);

      return vector3D<double>(sumX, sumY, sumZ);
  }

  vector3D<double> forceAccel::dispersionOnly(vector3D<double> &position)
  {
      double fieldX, fieldY, fieldZ;
      double sMatXY, sMatXZ, sMatYZ;
      double sumX, sumY, sumZ;
      double diagMatXX, diagMatYY, diagMatZZ;
      double xx, yy, zz, xy, xz, yz;
      double dx[3] = { 0.0, 0.0, 0.0 };
      double term1, term2, scalarComp;
      double r2, r_3, r_5, r6, charge;
      int index;
      int k;
      double dispPlusRad;

      fieldX = fieldY = fieldZ = 0;
      sMatXY = sMatXZ = sMatYZ = 0;
      diagMatXX = diagMatYY = diagMatZZ = 0;
      sumX = sumY = sumZ = 0;


      for(int n = 0; n < (signed)clusterGroups.size(); n++)
      {

          dx[0] = position.X - clusterPos[n].X;
          dx[1] = position.Y - clusterPos[n].Y;
          dx[2] = position.Z - clusterPos[n].Z;
          xx = dx[0]*dx[0]; yy = dx[1]*dx[1]; zz = dx[2]*dx[2];
          r2 = xx + yy + zz;

          dispPlusRad = dispersion_radius_ext + clusterRadius[n];
                  dispPlusRad *= dispPlusRad;

          if(!(r2 > dispPlusRad))
          {

              for (k = 0; k < (signed)clusterGroups[n].size(); k++)
              {
                      index = clusterGroups[n][k];

                      charge = mol->data[index].charge;
                      dx[0] = position.X - mol->data[index].coords.X;
                      dx[1] = position.Y - mol->data[index].coords.Y;
                      dx[2] = position.Z - mol->data[index].coords.Z;

                      r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

                      r_5 = 1.0 / (r2*r2*sqrt(r2));
                      r_3 = r_5*r2;

                      r_3 *= charge;
                      r_5 *= charge;

                      fieldX += dx[0]*r_3;
                      fieldY += dx[1]*r_3;
                      fieldZ += dx[2]*r_3;

                      diagMatXX += (r2 - 3 * dx[0]*dx[0])*r_5;
                      diagMatYY += (r2 - 3 * dx[1]*dx[1])*r_5;
                      diagMatZZ += (r2 - 3 * dx[2]*dx[2])*r_5;

                      sMatXY -= 3 * dx[0]*dx[1]*r_5;
                      sMatXZ -= 3 * dx[0]*dx[2]*r_5;
                      sMatYZ -= 3 * dx[1]*dx[2]*r_5;

                      r6 = r2*r2*r2;
                      term2 = mol->data[index].sigAvg6 / r6;
                      term1 = 2 * term2*term2;
                      scalarComp = mol->data[index].eps24_mass*(term1 - term2) / r2;

                      sumX += scalarComp*dx[0];
                      sumY += scalarComp*dx[1];
                      sumZ += scalarComp*dx[2];
              }
          }

          else
          {
            for (k = 0; k < (signed)clusterGroups[n].size(); k++)
             {
                      index = clusterGroups[n][k];

                      charge = mol->data[index].charge;
                      dx[0] = position.X - mol->data[index].coords.X;
                      dx[1] = position.Y - mol->data[index].coords.Y;
                      dx[2] = position.Z - mol->data[index].coords.Z;

                      r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

                      r_5 = 1.0 / (r2*r2*sqrt(r2));
                      r_3 = r_5*r2;

                      r_3 *= charge;
                      r_5 *= charge;

                      fieldX += dx[0]*r_3;
                      fieldY += dx[1]*r_3;
                      fieldZ += dx[2]*r_3;

                      diagMatXX += (r2 - 3 * dx[0]*dx[0])*r_5;
                      diagMatYY += (r2 - 3 * dx[1]*dx[1])*r_5;
                      diagMatZZ += (r2 - 3 * dx[2]*dx[2])*r_5;

                      sMatXY -= 3 * dx[0]*dx[1]*r_5;
                      sMatXZ -= 3 * dx[0]*dx[2]*r_5;
                      sMatYZ -= 3 * dx[1]*dx[2]*r_5;
             }
          }
      }

      sumX += dipoleScalar * (diagMatXX*fieldX + sMatXY*fieldY    + sMatXZ *fieldZ);
      sumY += dipoleScalar * (sMatXY*fieldX    + diagMatYY*fieldY + sMatYZ*fieldZ);
      sumZ += dipoleScalar * (sMatXZ*fieldX    + sMatYZ*fieldY    + diagMatZZ*fieldZ);

      return vector3D<double>(sumX, sumY, sumZ);
  }

  vector3D<double> forceAccel::chargeDipoleLJ(vector3D <double> &position)
    {
            double sumXX = 0.0, sumYY = 0.0, sumZZ = 0.0;
            double sumXY = 0.0, sumXZ = 0.0, sumYZ = 0.0;
            double xDiff, yDiff, zDiff, r2, r6, scalarComp;
            double dipX = 0.0, dipY = 0.0, dipZ = 0.0;
            double r, r3, r_3, r_5, r5;
            double xDiff2 = 0.0, yDiff2 = 0.0, zDiff2 = 0.0;

            double tempx = 0, tempy = 0, tempz = 0;
            double term1 = 0, term2 = 0;

            for (int n = 0; n < mol->numAtoms; n++)
            {

                    xDiff = position.X - mol->data[n].coords.X;
                    yDiff = position.Y - mol->data[n].coords.Y;
                    zDiff = position.Z - mol->data[n].coords.Z;
                    xDiff2 = xDiff*xDiff;
                    yDiff2 = yDiff*yDiff;
                    zDiff2 = zDiff*zDiff;


                    r2 = xDiff2 + yDiff2 + zDiff2;
                    r = sqrt(r2);
                    r3 = r2*r;
                    r_3 = 1.0 / r3;
                    r_5 = r_3 / r2;
                    r6 = r2*r2*r2;
                    r_3 *= mol->data[n].charge;
                    r_5 *= mol->data[n].charge;

                    term2 = mol->data[n].sigAvg6 / r6;
                    term1 = 2 * term2*term2;
                    scalarComp = mol->data[n].eps24_mass*(term1 - term2) / r2;

                    tempx += scalarComp*xDiff;
                    tempy += scalarComp*yDiff;
                    tempz += scalarComp*zDiff;


                    dipX += xDiff*r_3;
                    dipY += yDiff*r_3;
                    dipZ += zDiff*r_3;

                    sumXX += (3.0*xDiff2 - r2)*r_5;
                    sumYY += (3.0*yDiff2 - r2)*r_5;
                    sumZZ += (3.0*zDiff2 - r2)*r_5;
                    sumXY += 3*xDiff*yDiff*r_5;
                    sumXZ += 3*xDiff*zDiff*r_5;
                    sumYZ += 3*yDiff*zDiff*r_5;

            }


            tempx -= dipoleScalar * (sumXX*dipX + sumXY*dipY + sumXZ*dipZ);
            tempy -= dipoleScalar * (sumXY*dipX + sumYY*dipY + sumYZ*dipZ);
            tempz -= dipoleScalar * (sumXZ*dipX + sumYZ*dipY + sumZZ*dipZ);

            return vector3D<double>(tempx, tempy, tempz);
    }

  vector3D<double> forceAccel::dispersionNoCharge(vector3D<double> &position)
    {
        double sumX, sumY, sumZ;
        double xx, yy, zz, xy, xz, yz;
        double dx[3] = { 0.0, 0.0, 0.0 };
        double term1, term2, scalarComp;
        double r2, r6;
        int index;
        int k;
        double dispPlusRad;
        sumX = sumY = sumZ = 0;

        for(int n = 0; n < (signed)clusterGroups.size(); n++)
        {

            dx[0] = position.X - clusterPos[n].X;
            dx[1] = position.Y - clusterPos[n].Y;
            dx[2] = position.Z - clusterPos[n].Z;
            xx = dx[0]*dx[0]; yy = dx[1]*dx[1]; zz = dx[2]*dx[2];
            r2 = xx + yy + zz;

            dispPlusRad = dispersion_radius_ext + clusterRadius[n];
                    dispPlusRad *= dispPlusRad;

            if(!(r2 > dispPlusRad))
            {

                for (k = 0; k < (signed)clusterGroups[n].size(); k++)
                {


                        index = clusterGroups[n][k];



                        dx[0] = position.X - mol->data[index].coords.X;
                        dx[1] = position.Y - mol->data[index].coords.Y;
                        dx[2] = position.Z - mol->data[index].coords.Z;

                        r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];


                        r6 = r2*r2*r2;
                        term2 = mol->data[index].sigAvg6 / r6;
                        term1 = 2 * term2*term2;
                        scalarComp = mol->data[index].eps24_mass*(term1 - term2) / r2;

                        sumX += scalarComp*dx[0];
                        sumY += scalarComp*dx[1];
                        sumZ += scalarComp*dx[2];
                }
            }
            //cout << "  " << (signed)clusterGroups[n].size() << "  " ;
        }


        return vector3D<double>(sumX, sumY, sumZ);
    }

  vector3D<double> forceAccel::dispersionFull(vector3D <double> &position)
    {
            double xDiff, yDiff, zDiff, r2, r6, scalarComp;
            double xDiff2 = 0.0, yDiff2 = 0.0, zDiff2 = 0.0;
            double tempx = 0, tempy = 0, tempz = 0;
            double term1 = 0, term2 = 0;


            for (int n = 0; n < mol->numAtoms; n++)
            {
                    xDiff = position.X - mol->data[n].coords.X;
                    yDiff = position.Y - mol->data[n].coords.Y;
                    zDiff = position.Z - mol->data[n].coords.Z;
                    xDiff2 = xDiff*xDiff;
                    yDiff2 = yDiff*yDiff;
                    zDiff2 = zDiff*zDiff;
                    r2 = xDiff2 + yDiff2 + zDiff2;

                    r6 = r2*r2*r2;
                    term2 = mol->data[n].sigAvg6 / r6;
                    term1 = 2 * term2*term2;
                    scalarComp = mol->data[n].eps24_mass*(term1 - term2) / r2;

                    tempx += scalarComp*xDiff;
                    tempy += scalarComp*yDiff;
                    tempz += scalarComp*zDiff;
            }

            return vector3D<double>(tempx, tempy, tempz);
    }

 
  vector3D<double> forceAccel::getReducedVec(vector3D<double> &position)
  {

	  /*
	  if (!multipole_ext && !dispersionCutoff_ext)
		  return chargeDipoleLJ(position);
	  else
		  return multipoleCutoff(position);
	  */
	  //for (int n = 0; n < mol->numAtoms; n++)
		  //cout << mol->data[n].coords.X << "  " << mol->data[n].coords.Y << "  " << mol->data[n].coords.Z << endl;
	  //cin.get();

	  switch (forceType)
	  {
	  case 1:
	    return multipoleCutoff(position); //both multipole-expansion and LJ cutoff
	  case 2:
	    return multipoleOnly(position); //multipole-expansion and full LJ interactions
	  case 3:
	    return dispersionOnly(position); //exact electrostatics and LJ cutoff
	  case 4:
	    return chargeDipoleLJ(position); //exact electrostatics and full LJ interactions
	  case 5:
	    return dispersionNoCharge(position); //no electrostatics and LJ cutoff
	  case 6:
	    return dispersionFull(position); //no electrostatics and full LJ interactions
	  default:
		  return dispersionFull(position); //no electrostatics and full LJ interactions;
	  }
		  
  }
  
  vector3D<double> forceAccel::getVec(vector3D<double> &position)
  {
    return chargeDipoleLJ(position);
  }

  double forceAccel::getSecDerivMag(vector3D<double>& position)
  {
      //******************************************************
      //   This functions computes the second derivative
      //   of the magnitude of the force of the LJ potential
      //******************************************************
      double tempx = 0, tempy = 0, tempz = 0;
      double term1 = 0, term2 = 0, total = 0;
      double xDiff, yDiff, zDiff, r2, r6, r12;

      for (int n = 0; n < mol->numAtoms; n++)
      {

        xDiff = position.X - mol->data[n].coords.X;
        yDiff = position.Y - mol->data[n].coords.Y;
        zDiff = position.Z - mol->data[n].coords.Z;

        r2 = xDiff*xDiff + yDiff*yDiff + zDiff*zDiff;
        r6 = r2*r2*r2;
        r12 = r6*r6;

        term2 = mol->data[n].sigAvg6/r6;
        term1 = 26*term2*term2;
        total -=mol->data[n].eps24*(term1 - 6*term2)/r2;
      }
      return total;
  }

  double forceAccel::getPotential(vector3D<double> &position)
  {
            double total = 0;
            double term1 = 0, term2 = 0;
            double sumX = 0.0;
            double sumY = 0.0;
            double sumZ = 0.0;
            double xDiff, yDiff, zDiff, r2, r6, r_3;

            for (int n = 0; n < mol->numAtoms; n++)
              {

                xDiff = position.X - mol->data[n].coords.X;
                yDiff = position.Y - mol->data[n].coords.Y;
                zDiff = position.Z - mol->data[n].coords.Z;

                r2 = xDiff*xDiff + yDiff*yDiff + zDiff*zDiff;
                r6 = r2*r2*r2;

                term2 = mol->data[n].sigAvg6/r6;
                term1 = term2*term2;

                total += mol->data[n].eps4*(term1 - term2);

                if(charge_ext)
                {
                  r_3 = mol->data[n].charge/(r2*sqrt(r2));
                  sumX += xDiff*r_3;
                  sumY += yDiff*r_3;
                  sumZ += zDiff*r_3;
                }


              }

            if(charge_ext)
            {
              total -= dipoleScalar*0.5*(sumX*sumX + sumY*sumY + sumZ*sumZ);
            }

  return total;

  }

  void forceAccel::printClusters()
  {
          vector3D<double> testPoint(-18.0, 0.0, 0.0);
	  char buffer[100];
	  char type[4][10] = {"C", "N", "O", "F"};
	  string atomType;
	  int index = 0;
	  double x, y, z;
	  ofstream clusterFile;
	  clusterFile.open("cluster.xyz");
	  clusterFile << (mol->numAtoms + 1) << endl;
	  clusterFile << "Clustered Atoms\n";
	  clusterFile << "He999  " << testPoint.X << "  " << testPoint.Y << "  " << testPoint.Z << endl;

	  //cout << endl;
	  for (int n = 0; n < (signed)clusterGroups.size(); n++)
	  {
	          if((testPoint - clusterPos[n]).mag() <= 25)
	            atomType = "Ne" + ::to_string(n);
	          else
	            atomType = type[n % 4] + ::to_string(n);
	          //cout << "ClusterPos " << n << ": " << clusterPos[n].X << "  " << clusterPos[n].Y << "  " << clusterPos[n].Z << endl;
		  for (int i = 0; i < (signed)clusterGroups[n].size(); i++)
		  {
			  index = clusterGroups[n][i];
			  x = mol->data[index].coords.X;
			  y = mol->data[index].coords.Y;
			  z = mol->data[index].coords.Z;
			  sprintf(buffer, "%4s  %8.4f  %8.4f  %8.4f\n", atomType.c_str(), x, y, z);
			  clusterFile << buffer;
		  }
	  }
	  clusterFile.close();
  }

} /* namespace std */
