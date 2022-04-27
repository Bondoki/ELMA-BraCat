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

#ifndef LEMONADE_UPDATER_CREATOR_SLOWGROWTH_H
#define LEMONADE_UPDATER_CREATOR_SLOWGROWTH_H
/**
 * @file
 *
 * @class UpdaterCreatorSlowGrowth
 *
 * @brief create Updater for a brached molecule using the slow growth algorithm
 * 
 * @details A single branched molecule is created by successivels attachment of monomers 
 * to the structure with restriction of functionality f <= 3 and excluded volume.
 * For more details, please see: https://doi.org/10.1039/C4SM00711E
 *
 * @tparam IngredientsType
 *
 **/

#include <LeMonADE/updater/UpdaterAbstractCreate.h>
#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/updater/moves/MoveAddMonomerSc.h>

#include <LeMonADE/utility/FastBondset.h>

#include <map>
#include <vector>
#include <map>
#include <utility>
#include <functional>
#include <queue>
#include <vector>

#include "StatisticMoment.h"
#include "Histogram1D.h"

template<class IngredientsType>
class UpdaterCreatorSlowGrowth: public UpdaterAbstractCreate<IngredientsType>
{
	typedef UpdaterAbstractCreate<IngredientsType> BaseClass;

public:
	UpdaterCreatorSlowGrowth(IngredientsType& ingredients_, uint32_t _maxmonomers, uint32_t _box_x=128, uint32_t _box_y=128, uint32_t _box_z=128);

	virtual void initialize();
	virtual bool execute();
	virtual void cleanup();

	bool addMonomerToParentWithProbability(uint32_t parent_id, int32_t type);

private:
	using BaseClass::ingredients;

	using BaseClass::addMonomerToParent;
	using BaseClass::addSingleMonomer;
	using BaseClass::addMonomerAtPosition;
	using BaseClass::addMonomerInsideConnectedPair;
	using BaseClass::linearizeSystem;
	using BaseClass::moveSystem;
	using BaseClass::randomBondvector;

	//! linear chain length
	uint32_t number_of_monomers;

	//! simulation box sizes
	uint32_t boxX,boxY,boxZ;

	
	//! Maximum number of monomers inside the SG structure
	uint32_t maxmonomers;

	// RNG
	RandomNumberGenerators rng;

	// (small) bondset for creation of bonds
	FastBondset smallBondSet;
        
        //! provides a histogram of functionality 
        Histogram1D* HG_FunctionalityNodes; 
        
        int nrOfSamples;
        
};

/** 
 * @brief Constructor handling the new systems paramters
 *
 * @param ingredients_ a reference to the IngredientsType - mainly the system
 * @param L_ linear chain length
 * @param box_ size of the box
 */
template < class IngredientsType >
UpdaterCreatorSlowGrowth<IngredientsType>::UpdaterCreatorSlowGrowth(IngredientsType& ingredients_, uint32_t _maxmonomers, uint32_t _box_x, uint32_t _box_y, uint32_t _box_z):
BaseClass(ingredients_), maxmonomers(_maxmonomers), boxX(_box_x), boxY(_box_y), boxZ(_box_z)
{

    HG_FunctionalityNodes = new Histogram1D(0 - 0.5, 2000 - 0.5, 2000);
    nrOfSamples = 0;
}

/**
 * The initialize function handles the new systems information.
 *
 * @tparam IngredientsType Features used in the system. See Ingredients.
 */
template < class IngredientsType >
void UpdaterCreatorSlowGrowth<IngredientsType>::initialize(){
	std::cout << "initializeUpdaterCreatorChainWalkingTertiaryBondWalking" << std::endl;

	//set box size
	ingredients.setBoxX(boxX);
	ingredients.setBoxY(boxY);
	ingredients.setBoxZ(boxZ);

	//set periodicity
	ingredients.setPeriodicX(true);
	ingredients.setPeriodicY(true);
	ingredients.setPeriodicZ(true);

	//add Bondset
	// supress std::cout output of addBondset
	std::streambuf *old = std::cout.rdbuf(); // <-- save
	std::stringstream ss;
	std::cout.rdbuf (ss.rdbuf());       // <-- redirect
	ingredients.modifyBondset().clear();
	ingredients.modifyBondset().addBFMclassicBondset();
	std::cout.rdbuf (old);

	ingredients.modifyMolecules().resize(0);
	ingredients.modifyMolecules().clear();
	ingredients.modifyMolecules().addMonomer(boxX/2, boxY/2, boxZ/2);
	ingredients.modifyMolecules()[0].setAttributeTag(1);
	ingredients.modifyMolecules()[0].setMovableTag(false);
        
        ingredients.modifyMolecules().addMonomer(boxX/2+2, boxY/2, boxZ/2);
	ingredients.modifyMolecules()[1].setAttributeTag(2);
	
        ingredients.modifyMolecules().addMonomer(boxX/2+4, boxY/2, boxZ/2);
	ingredients.modifyMolecules()[2].setAttributeTag(1);
	
	ingredients.modifyMolecules().connect( 0, 1);
	ingredients.modifyMolecules().connect( 1, 2);
        
        //ingredients.modifyMolecules().addMonomer(boxX/2+2, boxY/2, boxZ/2);
	//ingredients.modifyMolecules()[1].setAttributeTag(2);
	//idxWalker = ingredients.getMolecules().size()-1;

	//ingredients.modifyMolecules().connect( 0, 1);

	ingredients.synchronize(ingredients);

	//create and fill small bondset
	smallBondSet.clear();

	smallBondSet.addBond(2,0,0, 0);
	smallBondSet.addBond(0,0,2, 1);
	smallBondSet.addBond(0,2,0, 2);
	smallBondSet.addBond(-2,0,0, 3);
	smallBondSet.addBond(0,0,-2, 4);
	smallBondSet.addBond(0,-2,0, 5);

	

	smallBondSet.updateLookupTable();

	// execute();
}

/**
 * Execution of the system creation
 *
 * @tparam IngredientsType Features used in the system. See Ingredients.
 */
template < class IngredientsType >
bool UpdaterCreatorSlowGrowth<IngredientsType>::execute(){
	std::cout << "executeUpdaterCreatorSlowGrowth" << std::endl;

	std::cout << "small set of bonds used for creation of monomers:" << std::endl;

	for(FastBondset::iterator it=smallBondSet.begin();it!=smallBondSet.end();++it)
	{
		std::cout << it->first << " -> " << it->second << std::endl;
	}

	std::cout << "Number of monomers: " << ingredients.getMolecules().size() << std::endl;

	for (int i = 0; i < ingredients.getMolecules().size(); i++)
	{
		std::cout << "Monomer " << i << " at : " << ingredients.getMolecules()[i] << " with tag: " <<  ingredients.getMolecules()[i].getAttributeTag() << std::endl;
	}
	
	uint32_t counter =0 ;
	// create the slow growth structure with Chain Walking

        nrOfSamples++;

        ingredients.modifyMolecules()[0].setMovableTag(true);

        
	uint32_t numMonomers = 3;
	do
	{
		// std::cout << "Number of monomers in the system: "  << ingredients.getMolecules().size() << std::endl;

		bool addingSuccessful=false;

		
		
                    
                // randomly select a monomer to add a monomer
                // idxWalker = position of monomer as index
                uint32_t idxWalker=rng.r250_rand32()%ingredients.getMolecules().size();
		// try to connect to new structure if functionality is smaller than 3
		if(ingredients.getMolecules().getNumLinks(idxWalker) < 3)
		{
			// std::cout << "Try to add monomer from: " << idxWalker;

			addingSuccessful=addMonomerToParentWithProbability(idxWalker, (ingredients.getMolecules()[idxWalker].getAttributeTag()%2)+1);

			if(addingSuccessful == true)
			{
				std::cout << "Number of monomers in the system: "  << ingredients.getMolecules().size() << std::endl;
				std::cout << "   -> equilibrated: " << (1000*counter) << " attempted monomer moves = " << (1000*counter/ingredients.getMolecules().size()) << " MCS"<< std::endl;
				counter =0;

				numMonomers++;
				// std::cout << " to monomer: " << (ingredients.getMolecules().size()-1) << " with tag " << (ingredients.getMolecules()[(ingredients.getMolecules().size()-1)].getAttributeTag()) << std::endl;
			}

		}

		
		
		// std::cout << " to " << idxWalker << std::endl;

		//move the system 1000 monomer sweeps
		//if(addingSuccessful == true)
		{
			//  std::cout << "move 1000 sweeps" << std::endl;
			MoveLocalSc move;

			for(int32_t m=0;m<1000;m++){
				move.init(ingredients);
				if(move.check(ingredients)==true){
					move.apply(ingredients);
				}
			}

			counter++;

		}
	}while (numMonomers != maxmonomers);

	
	//std::cout << "Creation done with " << ingredients.getMolecules().size() << " monomers and Re2e of " << (ingredients.getMolecules()[idxEnd]-ingredients.getMolecules()[0]) << std::endl;

        
        // create the histogram of occurence of functionality in structure    
        for (int idx = 0; idx < ingredients.getMolecules().size(); idx++)
        {
            HG_FunctionalityNodes->AddValue(ingredients.getMolecules().getNumLinks(idx));
        }
        
        
	// free the first monomer
	ingredients.modifyMolecules()[0].setMovableTag(true);

	std::cout << "Finial sync...";
	ingredients.synchronize();
	std::cout << "done." << std::endl;

	//std::cout << "Finial linearization...";
	//linearizeSystem();
	//ingredients.synchronize();
	//std::cout << "done." << std::endl;

        return true;
}


template<class IngredientsType>
bool UpdaterCreatorSlowGrowth<IngredientsType>::addMonomerToParentWithProbability(uint32_t parent_id, int32_t type){

	// collect all bond-vectors leading to an allowed adding of monomers
	std::vector<VectorInt3> allowedBondsForAdding;
	allowedBondsForAdding.clear();

	// test all bond-vectors in the small set and collect them
	for(FastBondset::iterator it=smallBondSet.begin();it!=smallBondSet.end();++it)
	{
		// set properties of add Monomer Move
		MoveAddMonomerSc<> addmove;
		addmove.init(ingredients);
		addmove.setTag(type);

		VectorInt3 bondvector(it->second);

		// set position of new monomer
		addmove.setPosition(ingredients.getMolecules()[parent_id]+bondvector);

		// check new position (excluded volume) etc
		if(addmove.check(ingredients)==true){
			//adding monomer is possible -> collect the vector
			allowedBondsForAdding.push_back(bondvector);
			//std::cout << "Adding at idx " << (allowedBondsForAdding.size()-1) << " the bond " << it->second  << " is possible" << std::endl;
		}

	}

	//std::cout << "For monomer with idx: " << parent_id  << " are " << allowedBondsForAdding.size() << " are possible." << std::endl;

	// Adding monomer with probability into the vicinity or return
	if(allowedBondsForAdding.size() != 0)
	{
		MoveAddMonomerSc<> addmove;
		addmove.init(ingredients);
		addmove.setTag(type);

		// chose bond-vector out of the allowed set
		uint32_t randomBondIdx = rng.r250_rand32() % allowedBondsForAdding.size();
		//std::cout << "Used random bond at idx: " << randomBondIdx << " == " << allowedBondsForAdding.at(randomBondIdx) << std::endl;

		VectorInt3 bondvector(allowedBondsForAdding.at(randomBondIdx));

		// set position of new monomer
		addmove.setPosition(ingredients.getMolecules()[parent_id]+bondvector);

		// check new position (excluded volume) -> should never be false
		if(addmove.check(ingredients)==true){

			/*
			// test the adding against the probability of adding
			if(rng.r250_drand() <= probability)
			{
				addmove.apply(ingredients);
				ingredients.modifyMolecules().connect( parent_id, (ingredients.getMolecules().size()-1) );
				return true;
			}
			else return false;
			*/


			addmove.apply(ingredients);
			ingredients.modifyMolecules().connect( parent_id, (ingredients.getMolecules().size()-1) );
			return true;

		}
		else throw std::runtime_error("Adding of monomer with incorrect conditions");
	}
	else
	{ 	// nothing to add
		return false;
	}

}

/**
 * Standard clean up.
 *
 * @tparam IngredientsType Features used in the system. See Ingredients.
 */
template < class IngredientsType >
void UpdaterCreatorSlowGrowth<IngredientsType>::cleanup(){

    // print results into a file
    // get the filename and path
    std::string filenameGeneral = ingredients.getName();
    // delete the .bfm in the name
    filenameGeneral.erase(ingredients.getName().length() - 4, ingredients.getName().length());

    /********************************/
    // construct a list
    std::vector < std::vector<double> > tmpResultsHG_FunctionalityNodes;

    // we have 6 columns and row
    int columns = 4;
    int rows = 1;

    // we have columns
    tmpResultsHG_FunctionalityNodes.resize(columns);

    // we have rows
    for (int i = 0; i < columns; i++)
        tmpResultsHG_FunctionalityNodes[i].resize(0);

    // fill the histogram
    for (int bin = 0; bin < HG_FunctionalityNodes->GetNrBins(); bin++) {
        if (HG_FunctionalityNodes->GetNrInBin(bin) != 0) {
            tmpResultsHG_FunctionalityNodes[0].push_back(HG_FunctionalityNodes->GetRangeInBin(bin));
            tmpResultsHG_FunctionalityNodes[1].push_back(HG_FunctionalityNodes->GetNrInBin(bin));
            tmpResultsHG_FunctionalityNodes[2].push_back(HG_FunctionalityNodes->GetNrInBinNormiert(bin) / HG_FunctionalityNodes->GetIntervallThickness());
            tmpResultsHG_FunctionalityNodes[3].push_back(HG_FunctionalityNodes->GetCumulativeNrInBinNormiert(bin) / HG_FunctionalityNodes->GetIntervallThickness());
            
        }
    }
    
    std::stringstream commentHG_FunctionalityNodes;
    commentHG_FunctionalityNodes << " File produced by analyzer UpdaterCreatorSlowGrowth" << std::endl
            << " Statistics of hyperstar structure providing the occurence of functionality f" << std::endl
            << " with " << ingredients.getMolecules().size() << " monomers" << std::endl
            << std::endl
            << " Total counts in all bins: " << HG_FunctionalityNodes->GetNrCounts() << std::endl
            << " Total number of samples: " << nrOfSamples << std::endl
            << " Total number of bins " << HG_FunctionalityNodes->GetNrBins() << std::endl
            << " Interval [" << HG_FunctionalityNodes->GetRangeInBin(0) << " ; " << HG_FunctionalityNodes->GetRangeInBin(HG_FunctionalityNodes->GetNrBins()) << "]" << std::endl
            << " Interval thickness dI = " << HG_FunctionalityNodes->GetIntervallThickness() << std::endl
            //<< " Attachment probability p: " << probabilityForInsertion << std::endl
            //<< " Walking Rate w: " << (1.0/probabilityForInsertion) << std::endl
            << std::endl
            << " f          ... functionality of node" << std::endl
            << " <c(f)>     ... total counts of functionality (ensemble average)" << std::endl
            << " <cnorm(f)> ... average normalized counts by all counts of functionality (ensemble sum)" << std::endl
            << " <snorm(f)> ... cumulative average normalized counts by all counts of functionality = sum_d <cnorm(d)> " << std::endl
            << "f <c(f)> <cnorm(f)> <snorm(f)>";

    //new filename
    std::string filenameHG_FunctionalityNodes = filenameGeneral + "_CreationPropertyHGFunctionality.dat";

    ResultFormattingTools::writeResultFile(filenameHG_FunctionalityNodes, this->ingredients, tmpResultsHG_FunctionalityNodes, commentHG_FunctionalityNodes.str());
    
    // delete all allocated memory
       
    delete HG_FunctionalityNodes; 
}


#endif /* LEMONADE_UPDATER_CREATOR_SLOWGROWTH_H */
