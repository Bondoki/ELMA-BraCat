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

#include <LeMonADE/utility/RandomNumberGenerators.h>
#include <LeMonADE/core/ConfigureSystem.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureAttributes.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
#include <LeMonADE/feature/FeatureFixedMonomers.h>
#include <LeMonADE/feature/FeatureNNInteractionSc.h>
#include <LeMonADE/utility/TaskManager.h>
#include <LeMonADE/updater/UpdaterReadBfmFile.h>
#include <LeMonADE/updater/UpdaterSimpleSimulator.h>

#include "catchorg/clara/clara.hpp"

int main(int argc, char* argv[])
{
  try{
	std::string infile  = "input.bfm";
	std::string outfile = "outfile.bfm";
	uint32_t max_mcs=100;
	uint32_t save_interval=100;
	
    bool showHelp = false;
    
    auto parser
    = clara::Opt( infile, "input (=input.bfm)" )
        ["-i"]["--infile"]
        ("BFM-file to load.")
        .required()
    | clara::Opt( outfile, "output (=outfile.bfm)" )
        ["-o"]["--outfile"]
        ("BFM-file to save.")
        .required()
    | clara::Opt( [&max_mcs](int const m)
        {
         if (m <= 0)
         {
            return clara::ParserResult::runtimeError("Simulation time must be greater than 0");
         }
         else
         {
            max_mcs = m;
            return clara::ParserResult::ok(clara::ParseResultType::Matched);
         }
        }, "max MCS(=100)" )
        ["-m"]["--max-mcs"]
        ("(required) specifies the total Monte-Carlo steps to simulate.")
        .required()
    | clara::Opt( [&save_interval](int const s)
        {
         if (s < 0)
         {
            return clara::ParserResult::runtimeError("Save intervall must be greater than 0");
         }
         else
         {
            save_interval = s;
            return clara::ParserResult::ok(clara::ParseResultType::Matched);
         }
      }, "save MCS(=100)")
        ["-s"]["--save-mcs"]
        ("(required) Save after every <integer> Monte-Carlo steps to the output file." )
        .required()
    | clara::Help( showHelp );
        
    auto result = parser.parse( clara::Args( argc, argv ) );
    if( !result ) {
    std::cerr << "Error in command line: " << result.errorMessage() << std::endl;
    exit(1);
    }
    else if(showHelp == true)
    {
        std::cout << "Simulator for the ScBFM with Ex.Vol and BondCheck" << std::endl
                  << "maximum number of connections per monomer is 6" << std::endl
                  << "Features used: FeatureBondset, FeatureAttributes, FeatureExcludedVolumeSc<FeatureLattice<bool> >" << std::endl
		          << "Updaters used: ReadFullBFMFile, SimpleSimulator" << std::endl
		          << "Analyzers used: WriteBfmFile" << std::endl;
        
        parser.writeToStream(std::cout);
        exit(0);
    }
    else
    {
        std::cout << "infile:        " << infile << std::endl
                  << "outfile:       " << outfile << std::endl
                  << "max_mcs:       " << max_mcs << std::endl
                  << "save_interval: " << save_interval << std::endl;
    }
       
	//seed the globally available random number generators
	RandomNumberGenerators rng;
	rng.seedAll();
	
	// FeatureExcludedVolume<> is equivalent to FeatureExcludedVolume<FeatureLattice<bool> >
	typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureAttributes< >,FeatureExcludedVolumeSc<>) Features;
	
	typedef ConfigureSystem<VectorInt3,Features, 6> Config;
	typedef Ingredients<Config> Ing;
	Ing myIngredients;

	TaskManager taskmanager;
	taskmanager.addUpdater(new UpdaterReadBfmFile<Ing>(infile,myIngredients,UpdaterReadBfmFile<Ing>::READ_LAST_CONFIG_SAVE),0);
	//here you can choose to use MoveLocalBcc instead. Careful though: no real tests made yet
	//(other than for latticeOccupation, valid bonds, frozen monomers...)
	taskmanager.addUpdater(new UpdaterSimpleSimulator<Ing,MoveLocalSc>(myIngredients,save_interval));

	taskmanager.addAnalyzer(new AnalyzerWriteBfmFile<Ing>(outfile,myIngredients));
	
	taskmanager.initialize();
	taskmanager.run(max_mcs/save_interval);
	taskmanager.cleanup();
	
	}
	catch(std::exception& err){std::cerr<<err.what();}
	return 0;
  
}

