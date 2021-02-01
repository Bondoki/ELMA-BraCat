

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

#include "UpdaterCreatorChainWalkingTertiaryBondWalking_CreationProperties.h"


int main(int argc, char* argv[])
{
  try{
	// FeatureExcludedVolume<> is equivalent to FeatureExcludedVolume<FeatureLattice<bool> >
        // FeatureFeatureAttributes<> is equivalent to FeatureAttributes<int32_t>
	typedef LOKI_TYPELIST_4(FeatureMoleculesIO, FeatureFixedMonomers, FeatureAttributes<>,FeatureExcludedVolumeSc<>) Features;
	
	typedef ConfigureSystem<VectorInt3,Features, 3> Config;
	typedef Ingredients<Config> Ing;
	Ing myIngredients;

	std::string filename="test.bfm";
	uint32_t number_of_monomers=100;
	double probability = 1.0;

	double probabilityTertiaryWalk = 1.0;

	int option_char(0);
	int box(128);
        
        int numSamples=250;

	//read in options by getopt
	while ((option_char = getopt (argc, argv, "f:n:p:t:b:s:h"))  != EOF){
		switch (option_char)
		{
		case 'f':
			filename=optarg;
			break;
		case 'n':
			number_of_monomers = atoi(optarg);
			break;
		case 't':
			probabilityTertiaryWalk = atof(optarg);
			break;
		case 'p':
			probability = atof(optarg);
			break;
		case 'b':
			box = atoi(optarg);
			break;
                case 's':
			numSamples = atoi(optarg);
			break;        
		case 'h': std::cout << "This program creates a branched structure on an evolving creation process." << std::endl;
                          std::cout << "It uses full excluded volume condition for creation." << std::endl;
                          std::cout << "It evaluates the average distance between creation events and number of steps." << std::endl;
                          std::cout << "Option 'p' represents the walking probability proportional to the walking rate~walking steps/time between reaction event" << std::endl;
                          std::cout << "Option 't' represents the walking probability over tertiary bond e.g. t=1.0 no restriction and t=0.0 full barrier." << std::endl;
                
		default:
			std::cerr << "Usage: " << argv[0] << " [-f filename] [-n number_of_monomers] [-p probability(=1.0)] [-t probabilityTertiaryWalk(=1.0)] [-s number_of_samples(=250)] [-b box size] \n";
			return 0;
		}
	}
	std::cout << "write out options: filename=" << filename <<  " number_of_monomers=" << number_of_monomers <<  " probability" << probability  << " probabilityTertiaryWalk:" << probabilityTertiaryWalk << " in box " << box <<std::endl;
        std::cout << "Number of samples for ensemble average: " << numSamples <<std::endl;

	//seed the globally available random number generators
    RandomNumberGenerators randomNumbers;
    randomNumbers.seedAll();

    myIngredients.setName(filename);


    UpdaterCreatorChainWalkingTertiaryBondWalking_CreationProperties<Ing> U_CW_CreationProperties(myIngredients, number_of_monomers, probability, probabilityTertiaryWalk, box, box, box);

        for(int i = 0; i < numSamples; i++)
        {
        U_CW_CreationProperties.initialize();
        U_CW_CreationProperties.execute();

        }
        U_CW_CreationProperties.cleanup();
        


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

