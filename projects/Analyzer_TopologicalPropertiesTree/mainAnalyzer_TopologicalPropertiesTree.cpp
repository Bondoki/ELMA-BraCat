/*--------------------------------------------------------------------------------
    ooo      L   attice-based  |
  o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
oo---0---oo  A   lgorithm and  |
 o/./|\.\o   D   evelopment    | Copyright (C) 2017,2021 by
  o/.|.\o    E   nvironment    | LeMonADE Principal Developers (Ron Dockhorn)
    ooo                        |
----------------------------------------------------------------------------------

This file is part of LeMonADE.

LeMonADE is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

LeMonADE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with LeMonADE.  If not, see <http://www.gnu.org/licenses/>.

--------------------------------------------------------------------------------*/

#include <cstring>
#include <sstream>

#include <LeMonADE/utility/RandomNumberGenerators.h>
#include <LeMonADE/core/ConfigureSystem.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureAttributes.h>
#include <LeMonADE/utility/TaskManager.h>
#include <LeMonADE/updater/UpdaterReadBfmFile.h>

#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
#include <LeMonADE/feature/FeatureFixedMonomers.h>


#include "Analyzer_TopologicalPropertiesTree.h"
//#include "UpdaterCreatorChainWalking.h"

int main(int argc, char* argv[])
{
  try{
	std::string infile;
	std::string outfile;
	uint32_t numMonoPerObjects=0;
	uint32_t numObjects=0;
	uint32_t startFrame=0;
	uint32_t diffFrame=0;
	uint32_t numFiles=0;
	uint32_t idxMiddleMonomer=0;
	uint32_t idxTerminalMonomer=0;
	std::string dstDir ="./";
	
	outfile="outfile.bfm";
	
	double probability = 1.0;

	if(!(argc==11)|| (argc==2 && strcmp(argv[1],"--help")==0 ))
	{
		std::string errormessage;
		errormessage="usage: ./AnalyzerMSDGeneralObjects input_filename NrOfMonoPerObject NrOfObjects startAtFrame DiffFrames NrOfFiles DstDir IdxMiddleMonomer IdxTerminalMonomer\n";
		errormessage+="\nSimple Analyzer for calculating the MSD of walker ON the structure\n";
		throw std::runtime_error(errormessage);
		
	}
	else
	{
		infile=argv[1];
		numMonoPerObjects=atoi(argv[2]);
		numObjects=atoi(argv[3]);
		startFrame=atoi(argv[4]);
		diffFrame=atoi(argv[5]);
		numFiles=atoi(argv[6]);
		dstDir=argv[7];
		idxMiddleMonomer=atoi(argv[8]);
		idxTerminalMonomer=atoi(argv[9]);
		probability=atof(argv[10]);

		if(startFrame < 1)
		{
			std::string errormessage;
			errormessage="StartFrame should be > 0 \n";
			throw std::runtime_error(errormessage);
		}

		if(diffFrame < 1)
		{
			std::string errormessage;
			errormessage="diffFrame should be > 0 \n";
			throw std::runtime_error(errormessage);
		}
	}
	
	//seed the globally available random number generators
	RandomNumberGenerators rng;
	rng.seedAll();
	
	// FeatureExcludedVolume<> is equivalent to FeatureExcludedVolume<FeatureLattice<bool> >
	//typedef LOKI_TYPELIST_4(FeatureMoleculesIO, FeatureFixedMonomers,FeatureAttributes,FeatureExcludedVolume<>) Features;
	typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureAttributes<>,FeatureExcludedVolumeSc<>) Features;

	typedef ConfigureSystem<VectorInt3,Features, 4> Config;
	typedef Ingredients<Config> Ing;
	Ing myIngredients;

	myIngredients.setName(infile);
	
	
	TaskManager taskmanager;
	taskmanager.addUpdater(new UpdaterReadBfmFile<Ing>(infile,myIngredients,UpdaterReadBfmFile<Ing>::READ_LAST_CONFIG_SAVE),0);

	taskmanager.addAnalyzer(new Analyzer_TopologicalPropertiesTree<Ing>(myIngredients, dstDir));
	

	taskmanager.initialize();
	taskmanager.run();
	taskmanager.cleanup();

	
	/*
	//double probability = 0.1;
	//int number_of_monomers=1024;
	//int box=512;

	//set box size
	myIngredients.setBoxX(128);
	myIngredients.setBoxY(128);
	myIngredients.setBoxZ(128);

		//set periodicity
	myIngredients.setPeriodicX(true);
	myIngredients.setPeriodicY(true);
	myIngredients.setPeriodicZ(true);

	//add Bondset
	myIngredients.modifyBondset().clear();
	myIngredients.modifyBondset().addBFMclassicBondset();

	myIngredients.modifyMolecules().clear();
	myIngredients.synchronize(myIngredients);

	myIngredients.modifyMolecules().addMonomer(0, 0, 0);
	myIngredients.modifyMolecules().addMonomer(-2, -2, 0);
	myIngredients.modifyMolecules().addMonomer(-2, 0, 0);
	myIngredients.modifyMolecules().addMonomer(-4, -2, 0);
	myIngredients.modifyMolecules().addMonomer(0, -2, 0);
	myIngredients.modifyMolecules().addMonomer(-6, -2, 0);
	myIngredients.modifyMolecules().addMonomer(-2, -4, 0);
	myIngredients.modifyMolecules().addMonomer(2, -1, 0);

	myIngredients.modifyMolecules().connect( 0, 2);
	//myIngredients.modifyMolecules().connect( 0, 7);
	myIngredients.modifyMolecules().connect( 1, 2);
	myIngredients.modifyMolecules().connect( 1, 3);
	myIngredients.modifyMolecules().connect( 1, 4);
	myIngredients.modifyMolecules().connect( 1, 6);
	myIngredients.modifyMolecules().connect( 3, 5);
	myIngredients.modifyMolecules().connect( 4, 7);

	myIngredients.synchronize(myIngredients);

	Analyzer_TopologicalPropertiesTree<Ing> A(myIngredients, numMonoPerObjects, numObjects, startFrame, diffFrame, dstDir, idxMiddleMonomer, idxTerminalMonomer);

	A.initialize();
	A.execute();
	A.cleanup();
	*/
	/*
	UpdaterCreatorChainWalking<Ing> U(myIngredients, number_of_monomers, probability, box, box, box);
	Analyzer_TopologicalPropertiesTree<Ing> A(myIngredients, numMonoPerObjects, numObjects, startFrame, diffFrame, dstDir, idxMiddleMonomer, idxTerminalMonomer);

	A.initialize();

	for(int f = 0; f <= 1; f++)
	{
		U.initialize();
		U.execute();
		A.execute();

	}

	U.cleanup();
	A.cleanup();
	*/
	/*
	UpdaterReadBfmFile<Ing> U_RBFM(infile,myIngredients,UpdaterReadBfmFile<Ing>::READ_STEPWISE);
	AnalyzerCalculateMSDGeneralObjects<Ing> A_MSD(myIngredients, numMonoPerObjects, numObjects, startFrame, diffFrame, dstDir, idxMiddleMonomer, idxTerminalMonomer);

	U_RBFM.initialize();
	A_MSD.initialize();

	while(U_RBFM.execute() == true)
	{
		A_MSD.execute();
	}

	for(int files = 2; files <= numFiles; files++)
	{
		std::string filenameGeneral=infile;
		// delete the .bfm in the name
		filenameGeneral.erase (infile.length()-4, infile.length());
		std::stringstream ss;
		ss << files;

		//new filename
		std::string filename = filenameGeneral + "_0" + ss.str() + ".bfm";
		std::cout << "read-in file " << filename << std::endl;

		// work-around for clearing the bond-set for multiple file read-in
		myIngredients.modifyBondset().clear();
		UpdaterReadBfmFile<Ing> U_RBFM(filename,myIngredients,UpdaterReadBfmFile<Ing>::READ_STEPWISE);
		U_RBFM.initialize();

		while(U_RBFM.execute() == true)
			{
				A_MSD.execute();
			}

	}

	U_RBFM.cleanup();
	A_MSD.cleanup();

	*/
	}
	catch(std::exception& err){std::cerr<<err.what();}
	return 0;
  
}

