/*--------------------------------------------------------------------------------
    ooo      L   attice-based  | LeMonADE: An Open Source Implementation of the
  o\.|./o    e   xtensible     |           Bond-Fluctuation-Model for Polymers
 o\.\|/./o   Mon te-Carlo      |           
oo---0---oo  A   lgorithm and  | ELMA-BraCat: Extension for branched catalysis 
 o/./|\.\o   D   evelopment    | Copyright (C) 2018,2021,2022 by
  o/.|.\o    E   nvironment    | LeMonADE Principal Developers
    ooo      ELMA-BraCat       | Ron Dockhorn
----------------------------------------------------------------------------------

This file is part of LeMonADE and ELMA-BraCat extension.

LeMonADE and ELMA-BraCat extension is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

LeMonADE  and ELMA-BraCat extension is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with LeMonADE and ELMA-BraCat extension. If not, see <http://www.gnu.org/licenses/>.

--------------------------------------------------------------------------------*/

#include <cstring>
#include <stdlib.h> //for atoi
#include <unistd.h> //for getopt

#include <LeMonADE/utility/RandomNumberGenerators.h>
#include <LeMonADE/core/ConfigureSystem.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureAttributes.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
#include <LeMonADE/feature/FeatureFixedMonomers.h>
#include <LeMonADE/utility/TaskManager.h>
#include <LeMonADE/updater/UpdaterReadBfmFile.h>
#include <LeMonADE/updater/UpdaterSimpleSimulator.h>
#include <LeMonADE/utility/DepthIteratorPredicates.h>
#include <LeMonADE/analyzer/AnalyzerWriteBfmFileSubGroup.h>

#include "UpdaterCreatorBarabasiAlbert.h"
//#include "Analyzer_EigenvaluesRouseMatrix.h"


int main(int argc, char* argv[])
{
  try{
	// FeatureExcludedVolume<> is equivalent to FeatureExcludedVolume<FeatureLattice<bool> >
        // FeatureFeatureAttributes<> is equivalent to FeatureAttributes<int32_t>
	typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureFixedMonomers, FeatureAttributes<>) Features;
	
	typedef ConfigureSystem<VectorInt3,Features, 2000> Config;
	typedef Ingredients<Config> Ing;
	Ing myIngredients;

	std::string filename="BarabasiAlbert.bfm";
        uint32_t numbermonomers = 128;
        
	int option_char(0);
	int box(128);
        
        int numSamples=1;

	//read in options by getopt
	while ((option_char = getopt (argc, argv, "f:n:b:s:h"))  != EOF){
		switch (option_char)
		{
		case 'f':
			filename=optarg;
			break;
		case 'n':
			numbermonomers = atof(optarg);
			break;
		case 'b':
			box = atoi(optarg);
			break;
                case 's':
			numSamples = atoi(optarg);
			break;        
		case 'h': std::cout << "This program creates a BarabasiAlbert branched structure." << std::endl;
                          std::cout << "It only adds one node with exactly one edge to the graph. The creation of functionality of nodes is limited to 2000." << std::endl;
                          //std::cout << "It evaluates the ideal Gaussian radius of gyration" << std::endl;
                          //std::cout << "Option 'p' represents the walking probability proportional to the walking rate~walking steps/time between reaction event" << std::endl;
                          //std::cout << "Option 't' represents the walking probability over tertiary bond e.g. t=1.0 no restriction and t=0.0 full barrier." << std::endl;
                
		default:
			std::cerr << "Usage: " << argv[0] << " [-f filename] [-n numbermonomers(=128)] [-s number_of_samples(=1)] [-b box size] \n";
			return 0;
		}
	}
	std::cout << "write out options: filename=" << filename <<  " maxmonomers=" << numbermonomers << " in box " << box <<std::endl;
        std::cout << "Number of samples for ensemble average: " << numSamples <<std::endl;

	//seed the globally available random number generators
    RandomNumberGenerators randomNumbers;
    randomNumbers.seedAll();

    myIngredients.setName(filename);


    UpdaterCreatorBarabasiAlbert<Ing> U_CW(myIngredients, numbermonomers, box, box, box);

    AnalyzerWriteBfmFile<Ing> A_Write(filename,myIngredients);

        
        U_CW.initialize();
        U_CW.execute();

        A_Write.initialize();
        A_Write.execute();

        
        U_CW.cleanup();
        A_Write.cleanup();



	//TaskManager taskmanager;
	//taskmanager.addUpdater(new UpdaterReadBfmFile<Ing>(infile,myIngredients,UpdaterReadBfmFile<Ing>::READ_LAST_CONFIG_SAVE),0);
	//here you can choose to use MoveLocalBcc instead. Careful though: no real tests made yet
	//(other than for latticeOccupation, valid bonds, frozen monomers...)
	//taskmanager.addUpdater(new UpdaterSimpleSimulator<Ing,MoveLocalSc>(myIngredients,save_interval));

	//taskmanager.addAnalyzer(new AnalyzerWriteBfmFile<Ing>(outfile,myIngredients));
	
	//taskmanager.initialize();
	//taskmanager.run(max_mcs/save_interval);
	//taskmanager.cleanup();
	
	}
	catch(std::exception& err){std::cerr<<err.what();}
	return 0;
  
}

