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

#ifndef LEMONADE_UPDATER_CREATOR_BGRW_TertiaryBondWalking_H
#define LEMONADE_UPDATER_CREATOR_BGRW_TertiaryBondWalking_H
/**
 * @file
 *
 * @class UpdaterCreatorBGRWTertiaryBondWalking
 *
 * @brief create Updater for a BGRW structures and related creation properties
 * 
 * @details A single BGRW structure is created in a cubix box of arbitrary size
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

#include "StatisticMoment.h"
#include "Histogram1D.h"

template<class IngredientsType>
class UpdaterCreatorBGRWTertiaryBondWalking: public UpdaterAbstractCreate<IngredientsType>
{
	typedef UpdaterAbstractCreate<IngredientsType> BaseClass;

public:
	UpdaterCreatorBGRWTertiaryBondWalking(IngredientsType& ingredients_, uint32_t _number_of_monomers, double _probability, double _probabilityTertiaryWalk, uint32_t _box_x=128, uint32_t _box_y=128, uint32_t _box_z=128);

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

	//! idxWalker = position of monomer as index
	uint32_t idxWalker;

	//! probability of connecting and adding monomer
	double probabilityForInsertion;

	//! probability to walk on/over tertiary monomer
	double probabilityTertiaryWalk;

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
UpdaterCreatorBGRWTertiaryBondWalking<IngredientsType>::UpdaterCreatorBGRWTertiaryBondWalking(IngredientsType& ingredients_, uint32_t _number_of_monomers, double _probabilityForInsertion, double _probabilityTertiaryWalk, uint32_t _box_x, uint32_t _box_y, uint32_t _box_z):
BaseClass(ingredients_), number_of_monomers(_number_of_monomers), probabilityForInsertion(_probabilityForInsertion), probabilityTertiaryWalk(_probabilityTertiaryWalk), boxX(_box_x), boxY(_box_y), boxZ(_box_z)
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
void UpdaterCreatorBGRWTertiaryBondWalking<IngredientsType>::initialize(){
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
	idxWalker = ingredients.getMolecules().size()-1;

	ingredients.modifyMolecules().connect( 0, 1);

	ingredients.synchronize(ingredients);

	//create and fill small bondset
	smallBondSet.clear();

	smallBondSet.addBond(2,0,0, 0);
	smallBondSet.addBond(0,0,2, 1);
	smallBondSet.addBond(0,2,0, 2);
	smallBondSet.addBond(-2,0,0, 3);
	smallBondSet.addBond(0,0,-2, 4);
	smallBondSet.addBond(0,-2,0, 5);

	/*  smallBondSet.addBond(2,1,0, 6);
   	smallBondSet.addBond(1,0,2, 7);
   	smallBondSet.addBond(0,2,1, 8);
   	smallBondSet.addBond(-2,1,0, 9);
   	smallBondSet.addBond(1,0,-2, 10);
   	smallBondSet.addBond(0,-2,1, 11);

   	smallBondSet.addBond(2,-1,0, 12);
   	smallBondSet.addBond(-1,0,2, 13);
   	smallBondSet.addBond(0,2,-1, 14);
   	smallBondSet.addBond(-2,-1,0, 15);
   	smallBondSet.addBond(-1,0,-2, 16);
   	smallBondSet.addBond(0,-2,-1, 17);

   	smallBondSet.addBond(1,2,0, 18);
   	smallBondSet.addBond(2,0,1, 19);
   	smallBondSet.addBond(0,1,2, 20);
   	smallBondSet.addBond(-1,2,0, 21);
   	smallBondSet.addBond(2,0,-1, 22);
   	smallBondSet.addBond(0,-1,2, 23);

   	smallBondSet.addBond(1,-2,0, 24);
   	smallBondSet.addBond(-2,0,1, 25);
   	smallBondSet.addBond(0,1,-2, 26);
   	smallBondSet.addBond(-1,-2,0, 27);
   	smallBondSet.addBond(-2,0,-1, 28);
   	smallBondSet.addBond(0,-1,-2, 29);

   	smallBondSet.addBond(2,1,1, 30);
   	smallBondSet.addBond(1,1,2, 31);
   	smallBondSet.addBond(1,2,1, 32);
   	smallBondSet.addBond(-2,1,1, 33);
   	smallBondSet.addBond(1,1,-2, 34);
   	smallBondSet.addBond(1,-2,1, 35);

   	smallBondSet.addBond(2,-1,1, 36);
   	smallBondSet.addBond(-1,1,2, 37);
   	smallBondSet.addBond(1,2,-1, 38);
   	smallBondSet.addBond(-2,-1,1, 39);
   	smallBondSet.addBond(-1,1,-2, 40);
   	smallBondSet.addBond(1,-2,-1, 41);

   	smallBondSet.addBond(2,1,-1, 42);
   	smallBondSet.addBond(1,-1,2, 43);
   	smallBondSet.addBond(-1,2,1, 44);
   	smallBondSet.addBond(-2,1,-1, 45);
   	smallBondSet.addBond(1,-1,-2, 46);
   	smallBondSet.addBond(-1,-2,1, 47);

   	smallBondSet.addBond(2,-1,-1, 48);
   	smallBondSet.addBond(-1,-1,2, 49);
   	smallBondSet.addBond(-1,2,-1, 50);
   	smallBondSet.addBond(-2,-1,-1, 51);
   	smallBondSet.addBond(-1,-1,-2, 52);
   	smallBondSet.addBond(-1,-2,-1, 53);
	 */

	smallBondSet.updateLookupTable();

	// execute();
}

/**
 * Execution of the system creation
 *
 * @tparam IngredientsType Features used in the system. See Ingredients.
 */
template < class IngredientsType >
bool UpdaterCreatorBGRWTertiaryBondWalking<IngredientsType>::execute(){
	std::cout << "executeUpdaterCreatorChainWalkingTertiaryBondWalking" << std::endl;

	std::cout << "Creating structure with " << number_of_monomers << " with pI" << probabilityForInsertion << " probability for insertion and pT" << probabilityTertiaryWalk << " probability for tertiary C-walk " <<std::endl;

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

	std::cout << "Walker at monomer index: " << idxWalker << std::endl << std::endl;

	uint32_t counter =0 ;
	// create the slow growth structure with Chain Walking

	uint32_t numMonomers = 2;
        
        nrOfSamples++;
        
	do
	{
		// std::cout << "Number of monomers in the system: "  << ingredients.getMolecules().size() << std::endl;

		bool addingSuccessful=false;

		double rngInsertion = rng.r250_drand();

		// try to insert monomer at given position with pI
		if( rngInsertion < probabilityForInsertion)
		{

		// try to connect to new structure if functionality is smaller than 3
		// if(ingredients.getMolecules().getNumLinks(idxWalker) < 3)
		{
			// std::cout << "Try to add monomer from: " << idxWalker;

			addingSuccessful=addMonomerToParentWithProbability(idxWalker, (ingredients.getMolecules()[idxWalker].getAttributeTag()%2)+1);

			if(addingSuccessful == true)
			{
				std::cout << "Number of monomers in the system: "  << ingredients.getMolecules().size() << std::endl;
				std::cout << "   -> equilibrated: " << (1000*counter) << " attempted monomer moves = " << (1000*counter/ingredients.getMolecules().size()) << " MCS"<< std::endl;
				counter =0;

				numMonomers++;
				// move Walker to new monomer
				idxWalker = ingredients.getMolecules().size()-1;
				// std::cout << " to monomer: " << (ingredients.getMolecules().size()-1) << " with tag " << (ingredients.getMolecules()[(ingredients.getMolecules().size()-1)].getAttributeTag()) << std::endl;
			}

		}

		}
		else // walking with probability pWalk = 1-PI
		{

		// move the walker on the structure
		uint32_t possible_directions = ingredients.getMolecules().getNumLinks(idxWalker);


		// RW-node has neighbors also select one between [0; f)
		uint32_t rand_dir = rng.r250_rand32()%possible_directions;

		if(rand_dir < 0)
			throw std::runtime_error("RNG corrupted");

		// check for tertiary C-barrier
		if( ingredients.getMolecules().getNumLinks(ingredients.getMolecules().getNeighborIdx(idxWalker, rand_dir)) < 3 )
		{
			// no barrier at all
			idxWalker = ingredients.getMolecules().getNeighborIdx(idxWalker, rand_dir);
		}
		else // check against tertiary C-barrier
		{
			if(rng.r250_drand() < probabilityTertiaryWalk)
			{
				idxWalker = ingredients.getMolecules().getNeighborIdx(idxWalker, rand_dir);
			}
		}

	}
		// std::cout << " to " << idxWalker << std::endl;

		//move the system 1000 monomer sweeps
		//if(addingSuccessful == true)
		{
			//  std::cout << "move 1000 sweeps" << std::endl;
                        // no need to relax structure in BGRW
                        /*
			MoveLocalSc move;

			for(int32_t m=0;m<1000;m++){
				move.init(ingredients);
				if(move.check(ingredients)==true){
					move.apply(ingredients);
				}
			}
                        */

			counter++;

		}
	}while (numMonomers != number_of_monomers);

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

	std::cout << "Finial linearization...";
	linearizeSystem();
	ingredients.synchronize();
	std::cout << "done." << std::endl;

        return true;
}


template<class IngredientsType>
bool UpdaterCreatorBGRWTertiaryBondWalking<IngredientsType>::addMonomerToParentWithProbability(uint32_t parent_id, int32_t type){

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
void UpdaterCreatorBGRWTertiaryBondWalking<IngredientsType>::cleanup(){

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
    commentHG_FunctionalityNodes << " File produced by analyzer UpdaterCreatorBGRWTertiaryBondWalking" << std::endl
            << " Statistics of BGRW structure providing the occurence of functionality f" << std::endl
            << " with " << ingredients.getMolecules().size() << " monomers" << std::endl
            << std::endl
            << " Total counts in all bins: " << HG_FunctionalityNodes->GetNrCounts() << std::endl
            << " Total number of samples: " << nrOfSamples << std::endl
            << " Total number of bins " << HG_FunctionalityNodes->GetNrBins() << std::endl
            << " Interval [" << HG_FunctionalityNodes->GetRangeInBin(0) << " ; " << HG_FunctionalityNodes->GetRangeInBin(HG_FunctionalityNodes->GetNrBins()) << "]" << std::endl
            << " Interval thickness dI = " << HG_FunctionalityNodes->GetIntervallThickness() << std::endl
            << " Attachment probability p: " << probabilityForInsertion << std::endl
            << " Walking Rate w: " << (1.0/probabilityForInsertion) << std::endl
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


#endif /* LEMONADE_UPDATER_CREATOR_BGRW_TertiaryBondWalking_H */
