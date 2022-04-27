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

#ifndef LEMONADE_UPDATER_CREATOR_CHAIN_WALKING_TertiaryBondWalking_CreationProperties_H
#define LEMONADE_UPDATER_CREATOR_CHAIN_WALKING_TertiaryBondWalking_CreationProperties_H
/**
 * @file
 *
 * @class UpdaterCreatorChainWalkingTertiaryBondWalking_CreationProperties
 *
 * @brief create Updater for a CW structures and related creation properties
 * 
 * @details A single CW structure is created in a cubix box of arbitrary size
 *
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
#include <iostream>
#include <string>
#include <utility>
#include <functional>
#include <queue>

#include "StatisticMoment.h"
#include "Histogram1D.h"

template<class IngredientsType>
class UpdaterCreatorChainWalkingTertiaryBondWalking_CreationProperties: public UpdaterAbstractCreate<IngredientsType>
{
	typedef UpdaterAbstractCreate<IngredientsType> BaseClass;

public:
	UpdaterCreatorChainWalkingTertiaryBondWalking_CreationProperties(IngredientsType& ingredients_, uint32_t _number_of_monomers, double _probability, double _probabilityTertiaryWalk, uint32_t _box_x=128, uint32_t _box_y=128, uint32_t _box_z=128);

	virtual void initialize();
	virtual bool execute();
	virtual void cleanup();

	bool addMonomerToParentWithProbability(uint32_t parent_id, int32_t type);

        typedef uint32_t NodeIdx;
        typedef std::map<NodeIdx, int> Nodelist;
        typedef std::map<NodeIdx, Nodelist> Graph;
        //typedef std::pair<NodeIdx, Nodelist> Node;
        typedef std::pair<int, NodeIdx> Edge; //! (tentative distance, idxNode)
        typedef std::vector<NodeIdx> NodeVector;

        void dijkstra(Graph &graph, NodeIdx source, Nodelist &distance);
        int calculateTopologicalDistance(NodeIdx, NodeIdx);
        
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
        
        
        int nrOfSamples;
        
        StatisticMoment WalksBetweenAttemptedInsertion;
        StatisticMoment WalksBetweenPerformedInsertion;
        
        
        Histogram1D* HG_WalksBetweenAttemptedInsertion;
        Histogram1D* HG_WalksBetweenPerformedInsertion;
        
        StatisticMoment TopologicalDistanceBetweenAttemptedInsertion;
        StatisticMoment TopologicalDistanceBetweenPerformedInsertion;
        //! idxWalker = position of monomer as index
	uint32_t idxWalkerBeforePerformedInsertion;
        uint32_t idxWalkerBeforeAttemptedInsertion;
        
        Histogram1D* HG_TopologicalDistanceBetweenAttemptedInsertion;
        Histogram1D* HG_TopologicalDistanceBetweenPerformedInsertion;
};

/** 
 * @brief Constructor handling the new systems paramters
 *
 * @param ingredients_ a reference to the IngredientsType - mainly the system
 * @param L_ linear chain length
 * @param box_ size of the box
 */
template < class IngredientsType >
UpdaterCreatorChainWalkingTertiaryBondWalking_CreationProperties<IngredientsType>::UpdaterCreatorChainWalkingTertiaryBondWalking_CreationProperties(IngredientsType& ingredients_, uint32_t _number_of_monomers, double _probabilityForInsertion, double _probabilityTertiaryWalk, uint32_t _box_x, uint32_t _box_y, uint32_t _box_z):
BaseClass(ingredients_), number_of_monomers(_number_of_monomers), probabilityForInsertion(_probabilityForInsertion), probabilityTertiaryWalk(_probabilityTertiaryWalk), boxX(_box_x), boxY(_box_y), boxZ(_box_z)
{
  WalksBetweenAttemptedInsertion.clear();
  WalksBetweenPerformedInsertion.clear();
  
  HG_WalksBetweenAttemptedInsertion = new Histogram1D(0 - 0.5, 20000 - 0.5, 20000);
  HG_WalksBetweenPerformedInsertion = new Histogram1D(0 - 0.5, 20000 - 0.5, 20000);
  
  nrOfSamples = 0;
  
  TopologicalDistanceBetweenPerformedInsertion.clear();
  TopologicalDistanceBetweenAttemptedInsertion.clear();
  
  HG_TopologicalDistanceBetweenAttemptedInsertion = new Histogram1D(0 - 0.5, 20000 - 0.5, 20000);
  HG_TopologicalDistanceBetweenPerformedInsertion = new Histogram1D(0 - 0.5, 20000 - 0.5, 20000);
}

/**
 * The initialize function handles the new systems information.
 *
 * @tparam IngredientsType Features used in the system. See Ingredients.
 */
template < class IngredientsType >
void UpdaterCreatorChainWalkingTertiaryBondWalking_CreationProperties<IngredientsType>::initialize(){
	std::cout << "initializeUpdaterCreatorChainWalkingTertiaryBondWalking_CreationProperties" << std::endl;

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
        idxWalkerBeforePerformedInsertion = ingredients.getMolecules().size()-1;
        idxWalkerBeforeAttemptedInsertion = ingredients.getMolecules().size()-1;

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
bool UpdaterCreatorChainWalkingTertiaryBondWalking_CreationProperties<IngredientsType>::execute(){
	std::cout << "executeUpdaterCreatorChainWalkingTertiaryBondWalking_CreationProperties" << std::endl;

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
        
        uint32_t WalksBetweenPerformedInsertion_counter = 0;
        uint32_t WalksBetweenAttemptedInsertion_counter = 0;
	do
	{
		// std::cout << "Number of monomers in the system: "  << ingredients.getMolecules().size() << std::endl;

		bool addingSuccessful=false;

		double rngInsertion = rng.r250_drand();

		// try to insert monomer at given position with pI
		if( rngInsertion < probabilityForInsertion)
		{
                  WalksBetweenAttemptedInsertion.AddValue(WalksBetweenAttemptedInsertion_counter);
                  
                  // add value to the histogram
                  HG_WalksBetweenAttemptedInsertion->AddValue(WalksBetweenAttemptedInsertion_counter);
                  
                  WalksBetweenAttemptedInsertion_counter = 0;
                  
                  // +1 as virtually treated to be attached at new position
                  int td = calculateTopologicalDistance(idxWalkerBeforeAttemptedInsertion, idxWalker)+1;
                  TopologicalDistanceBetweenAttemptedInsertion.AddValue(td);
                  
                  HG_TopologicalDistanceBetweenAttemptedInsertion->AddValue(td);
                  
                  idxWalkerBeforeAttemptedInsertion=idxWalker;
                                

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
                                
                               
                                WalksBetweenPerformedInsertion.AddValue(WalksBetweenPerformedInsertion_counter);
                                // add value to the histogram
                                HG_WalksBetweenPerformedInsertion->AddValue(WalksBetweenPerformedInsertion_counter);

                                WalksBetweenPerformedInsertion_counter=0;

				numMonomers++;
				// move Walker to new monomer
				idxWalker = ingredients.getMolecules().size()-1;
				// std::cout << " to monomer: " << (ingredients.getMolecules().size()-1) << " with tag " << (ingredients.getMolecules()[(ingredients.getMolecules().size()-1)].getAttributeTag()) << std::endl;
			        //std::cout << "Calc Perf" << idxWalkerBeforePerformedInsertion << " to " <<idxWalker  << std::endl; 
                                //std::cout << "Calc Perf" << calculateTopologicalDistance(idxWalkerBeforePerformedInsertion, idxWalker) << std::endl;
                                int td = calculateTopologicalDistance(idxWalkerBeforePerformedInsertion, idxWalker);
                                TopologicalDistanceBetweenPerformedInsertion.AddValue(td);
                                HG_TopologicalDistanceBetweenPerformedInsertion->AddValue(td);
                                idxWalkerBeforePerformedInsertion=idxWalker;
                                
                                // correct the position in case the attachment was successful
                                idxWalkerBeforeAttemptedInsertion=idxWalker;
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
		
		WalksBetweenPerformedInsertion_counter++;
                WalksBetweenAttemptedInsertion_counter++;

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
	}while (numMonomers != number_of_monomers);

	//std::cout << "Creation done with " << ingredients.getMolecules().size() << " monomers and Re2e of " << (ingredients.getMolecules()[idxEnd]-ingredients.getMolecules()[0]) << std::endl;

	// free the first monomer
	ingredients.modifyMolecules()[0].setMovableTag(true);

	std::cout << "Finial sync...";
	ingredients.synchronize();
	std::cout << "done." << std::endl;

	std::cout << "Finial linearization...";
	linearizeSystem();
	ingredients.synchronize();
	std::cout << "done." << std::endl << std::endl;
        
        std::cout << "MeanWalksBetweenInsertionEventAndRealInsertion: " << WalksBetweenPerformedInsertion.ReturnM1() << " " << WalksBetweenAttemptedInsertion.ReturnM1() << std::endl;

        return true;
}

template<class IngredientsType>
int UpdaterCreatorChainWalkingTertiaryBondWalking_CreationProperties<IngredientsType>::calculateTopologicalDistance(NodeIdx StartNode, NodeIdx EndNode)
{
  Graph treeGraph;

  Nodelist degreeNode; // Degree of node

    for (int k = 0; k < ingredients.getMolecules().size(); k++) {
        for (int l = 0; l < ingredients.getMolecules().getNumLinks(k); l++) {
            // connect from idxAA to idxBB with value ZZ
            // all (symmetric) connections
            treeGraph[k][ingredients.getMolecules().getNeighborIdx(k, l)] = 1;

        }
        // fill with degree of nodes
        degreeNode[k] = ingredients.getMolecules().getNumLinks(k);
    }
    
    Nodelist topo_distance;
    dijkstra(treeGraph, StartNode, topo_distance);
    
    std::cout << "startNode -> node => distance" << std::endl;
    
    Nodelist::iterator it;
    for (it = topo_distance.begin(); it != topo_distance.end(); ++it) {
      std::cout << StartNode << " -> " << it->first << "  => " << it->second << std::endl;
    }
    
    int distance = 0;
    Nodelist::iterator itEN = topo_distance.find(EndNode);
    if (itEN != topo_distance.end())
      distance = itEN->second;
    else
      {
        std::stringstream errormessage;
        errormessage << "UpdaterCreatorChainWalkingTertiaryBondWalking_CreationProperties::execute()...Topological distance is undefined\n"
                << "This should NEVER happen as we have only a single tree\n";

        throw std::runtime_error(errormessage.str());
      }
    
    return distance;
}


template<class IngredientsType>
void UpdaterCreatorChainWalkingTertiaryBondWalking_CreationProperties<IngredientsType>::dijkstra(Graph &graph, NodeIdx source, Nodelist &distance) {

    distance.clear(); //! clear all tentative distance information

    //! list of all nodes to be visit and addressed already with least distance on top
    std::priority_queue<Edge, std::vector<Edge>, std::greater<Edge> > queueNode;

    // starting node with tentative distance=0, starting node = source
    queueNode.push(Edge(0, source));

    while (!queueNode.empty()) {

        //get the element with least tentative distance
        Edge tmped = queueNode.top();

        // access the node index
        NodeIdx tmpnl = tmped.second;

        // removes the top element
        queueNode.pop();

        // if the node never visited before
        if (distance.count(tmpnl) == 0) {

            // tentative distance to the recent node
            int dist = tmped.first;

            // set the tentative distance to the node
            distance[tmpnl] = dist;

            // new subgraph of all neighbors
            Nodelist tempgraph = graph[tmpnl];

            Nodelist::iterator it;
            for (it = tempgraph.begin(); it != tempgraph.end(); ++it) {
                int distint = it->second;
                NodeIdx distlabel = it->first;
                queueNode.push(Edge(dist + distint, distlabel));
            }

        }
    }

}

template<class IngredientsType>
bool UpdaterCreatorChainWalkingTertiaryBondWalking_CreationProperties<IngredientsType>::addMonomerToParentWithProbability(uint32_t parent_id, int32_t type){

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
void UpdaterCreatorChainWalkingTertiaryBondWalking_CreationProperties<IngredientsType>::cleanup(){
  
    std::cout << "MeanWalksBetweenInsertionAttemptedAndPerformedInsertion: " << WalksBetweenPerformedInsertion.ReturnM1() << " " << WalksBetweenAttemptedInsertion.ReturnM1() << std::endl;
    
    // print results into a file
    // get the filename and path
    std::string filenameGeneral = ingredients.getName();
    // delete the .bfm in the name
    filenameGeneral.erase(ingredients.getName().length() - 4, ingredients.getName().length());

    // construct a list
    std::vector < std::vector<double> > tmpResultsWalking;

    // we have 3 columns and row
    uint32_t columns = 5;
    uint32_t rows = 1;

    // we have columns
    tmpResultsWalking.resize(columns);

    // we have rows
    for (int i = 0; i < columns; i++)
        tmpResultsWalking[i].resize(0);

    tmpResultsWalking[0].push_back(probabilityForInsertion);
    tmpResultsWalking[1].push_back(WalksBetweenAttemptedInsertion.ReturnM1());
    tmpResultsWalking[2].push_back(WalksBetweenPerformedInsertion.ReturnM1());
    tmpResultsWalking[3].push_back(WalksBetweenAttemptedInsertion.ReturnN());
    tmpResultsWalking[4].push_back(WalksBetweenPerformedInsertion.ReturnN());
    
    std::stringstream comment;
    comment << " File produced by analyzer UpdaterCreatorChainWalkingTertiaryBondWalking_CreationProperties" << std::endl
            << " Statistics of creation process involving the average number of steps (walks) along the structure" << std::endl
            << " between an attempted (and maybe successful) adding of monomer Sattempt" << std::endl
            << " between an performed adding of monomer Sperf" << std::endl
            << " with " << ingredients.getMolecules().size() << " monomers" << std::endl
            << " Total number of samples: " << nrOfSamples << std::endl
            << std::endl
            << " p          ... attachment probability" << std::endl
            << " Sattempted ... average steps (walks) between attempted attachment" << std::endl
            << " Sperformed ... average steps (walks) between performed attachment" << std::endl
            << " samplesSat ... total number of samples for Sattempted" << std::endl
            << " samplesSpe ... total number of samples for Sperformed" << std::endl
            << "p Sattempted Sperformed samplesSat samplesSpe";

    //new filename
    std::string filenameWalking = filenameGeneral + "_CreationPropertyWalks.dat";

    ResultFormattingTools::writeResultFile(filenameWalking, this->ingredients, tmpResultsWalking, comment.str());

    /*******************/
    
    // construct a list
    std::vector < std::vector<double> > tmpResultsHG_WalkingAttempted;

    // we have 6 columns and row
    columns = 4;
    rows = 1;

    // we have columns
    tmpResultsHG_WalkingAttempted.resize(columns);

    // we have rows
    for (int i = 0; i < columns; i++)
        tmpResultsHG_WalkingAttempted[i].resize(0);

    // fill the histogram
    for (int bin = 0; bin < HG_WalksBetweenAttemptedInsertion->GetNrBins(); bin++) {
        if (HG_WalksBetweenAttemptedInsertion->GetNrInBin(bin) != 0) {
            tmpResultsHG_WalkingAttempted[0].push_back(HG_WalksBetweenAttemptedInsertion->GetRangeInBin(bin));
            tmpResultsHG_WalkingAttempted[1].push_back(HG_WalksBetweenAttemptedInsertion->GetNrInBin(bin));
            tmpResultsHG_WalkingAttempted[2].push_back(HG_WalksBetweenAttemptedInsertion->GetNrInBinNormiert(bin) / HG_WalksBetweenAttemptedInsertion->GetIntervallThickness());
            tmpResultsHG_WalkingAttempted[3].push_back(HG_WalksBetweenAttemptedInsertion->GetCumulativeNrInBinNormiert(bin) / HG_WalksBetweenAttemptedInsertion->GetIntervallThickness());
            
        }
    }
    
    std::stringstream commentHG_WalkingAttempted;
    commentHG_WalkingAttempted << " File produced by analyzer UpdaterCreatorChainWalkingTertiaryBondWalking_CreationProperties" << std::endl
            << " Statistics of creation process involving the average number of steps (walks) along the structure" << std::endl
            << " between an attempted (and maybe successful) adding of monomer Sattempted" << std::endl
            << " between an performed adding of monomer Sperformed" << std::endl
            << " with " << ingredients.getMolecules().size() << " monomers" << std::endl
            << " Total counts in all bins: " << HG_WalksBetweenAttemptedInsertion->GetNrCounts() << std::endl
            << " Total number of samples: " << nrOfSamples << std::endl
            << " Total number of bins " << HG_WalksBetweenAttemptedInsertion->GetNrBins() << std::endl
            << " Interval [" << HG_WalksBetweenAttemptedInsertion->GetRangeInBin(0) << " ; " << HG_WalksBetweenAttemptedInsertion->GetRangeInBin(HG_WalksBetweenAttemptedInsertion->GetNrBins()) << "]" << std::endl
            << " Interval thickness dI = " << HG_WalksBetweenAttemptedInsertion->GetIntervallThickness() << std::endl
            << " Attachment probability p: " << probabilityForInsertion << std::endl
            << " Average walk length between attempted attachment Sattempted: " << WalksBetweenAttemptedInsertion.ReturnM1() << std::endl
            << " Average walk length between successful attachment Sperformed: " << WalksBetweenPerformedInsertion.ReturnM1()  << std::endl
            << std::endl
            << " d ... steps (walk length) between performed attachment" << std::endl
            << " <c(d)> ... total counts at walk length (ensemble average)" << std::endl
            << " <cnorm(d)> ... average normalized counts by all counts at walk length (ensemble sum)" << std::endl
            << " <snorm(d)> ... cumulative average normalized counts by all counts at walk length d = sum_d <cnorm(d)> " << std::endl
            << "d <c(d)> <cnorm(d)> <snorm(d)>";

    //new filename
    std::string filenameHG_WalkingAttempted = filenameGeneral + "_CreationPropertyWalksHGWalkingAttempted.dat";

    ResultFormattingTools::writeResultFile(filenameHG_WalkingAttempted, this->ingredients, tmpResultsHG_WalkingAttempted, commentHG_WalkingAttempted.str());
    
    
    /********************************/
    
    // construct a list
    std::vector < std::vector<double> > tmpResultsHG_WalkingPerformed;

    // we have 6 columns and row
    columns = 4;
    rows = 1;

    // we have columns
    tmpResultsHG_WalkingPerformed.resize(columns);

    // we have rows
    for (int i = 0; i < columns; i++)
        tmpResultsHG_WalkingPerformed[i].resize(0);

    // fill the histogram
    for (int bin = 0; bin < HG_WalksBetweenPerformedInsertion->GetNrBins(); bin++) {
        if (HG_WalksBetweenPerformedInsertion->GetNrInBin(bin) != 0) {
            tmpResultsHG_WalkingPerformed[0].push_back(HG_WalksBetweenPerformedInsertion->GetRangeInBin(bin));
            tmpResultsHG_WalkingPerformed[1].push_back(HG_WalksBetweenPerformedInsertion->GetNrInBin(bin));
            //tmpResultsHG_WalkingPerformed[1].push_back(HG_WalksBetweenPerformedInsertion->GetNrInBin(bin) / (1.0 * nrOfSamples));
            tmpResultsHG_WalkingPerformed[2].push_back(HG_WalksBetweenPerformedInsertion->GetNrInBinNormiert(bin) / HG_WalksBetweenPerformedInsertion->GetIntervallThickness());
            //tmpResultsHG_WalkingPerformed[3].push_back(HG_WalksBetweenPerformedInsertion->GetCumulativeNrInBin(bin) / (1.0 * nrOfSamples));
            tmpResultsHG_WalkingPerformed[3].push_back(HG_WalksBetweenPerformedInsertion->GetCumulativeNrInBinNormiert(bin) / HG_WalksBetweenPerformedInsertion->GetIntervallThickness());
            //tmpResultsHG_WalkingPerformed[5].push_back(HG_WalksBetweenPerformedInsertion->GetNrInBin(bin));

        }
    }
    
    std::stringstream commentHG_WalkingPerformed;
    commentHG_WalkingPerformed << " File produced by analyzer UpdaterCreatorChainWalkingTertiaryBondWalking_CreationProperties" << std::endl
            << " Statistics of creation process involving the average number of steps (walks) along the structure" << std::endl
            << " between an attempted (and maybe successful) adding of monomer Sattempted" << std::endl
            << " between an performed adding of monomer Sperformed" << std::endl
            << " with " << ingredients.getMolecules().size() << " monomers" << std::endl
            << " Total counts in all bins: " << HG_WalksBetweenPerformedInsertion->GetNrCounts() << std::endl
            << " Total number of samples: " << nrOfSamples << std::endl
            << " Total number of bins " << HG_WalksBetweenPerformedInsertion->GetNrBins() << std::endl
            << " Interval [" << HG_WalksBetweenPerformedInsertion->GetRangeInBin(0) << " ; " << HG_WalksBetweenPerformedInsertion->GetRangeInBin(HG_WalksBetweenPerformedInsertion->GetNrBins()) << "]" << std::endl
            << " Interval thickness dI = " << HG_WalksBetweenPerformedInsertion->GetIntervallThickness() << std::endl
            << " Attachment probability p: " << probabilityForInsertion << std::endl
            << " Average walk length between attempted attachment Sattempted: " << WalksBetweenAttemptedInsertion.ReturnM1() << std::endl
            << " Average walk length between successful attachment Sperformed: " << WalksBetweenPerformedInsertion.ReturnM1()  << std::endl
            << std::endl
            << " d ... steps (walk length) between performed attachment" << std::endl
            << " <c(d)> ... total counts at walk length (ensemble average)" << std::endl
            << " <cnorm(d)> ... average normalized counts by all counts at walk length (ensemble sum)" << std::endl
            << " <snorm(d)> ... cumulative average normalized counts by all counts at walk length d = sum_d <cnorm(d)> " << std::endl
            << "d <c(d)> <cnorm(d)> <snorm(d)>";

    //new filename
    std::string filenameHG_WalkingPerformed = filenameGeneral + "_CreationPropertyWalksHGWalkingPerformed.dat";

    ResultFormattingTools::writeResultFile(filenameHG_WalkingPerformed, this->ingredients, tmpResultsHG_WalkingPerformed, commentHG_WalkingPerformed.str());
    
    
    
    /****************************/
     // construct a list
    std::vector < std::vector<double> > tmpResultsTopologicalDistance;

    // we have 3 columns and row
    columns = 5;
    rows = 1;

    // we have columns
    tmpResultsTopologicalDistance.resize(columns);

    // we have rows
    for (int i = 0; i < columns; i++)
        tmpResultsTopologicalDistance[i].resize(0);

    tmpResultsTopologicalDistance[0].push_back(probabilityForInsertion);
    tmpResultsTopologicalDistance[1].push_back(TopologicalDistanceBetweenAttemptedInsertion.ReturnM1());
    tmpResultsTopologicalDistance[2].push_back(TopologicalDistanceBetweenPerformedInsertion.ReturnM1());
    tmpResultsTopologicalDistance[3].push_back(TopologicalDistanceBetweenAttemptedInsertion.ReturnN());
    tmpResultsTopologicalDistance[4].push_back(TopologicalDistanceBetweenPerformedInsertion.ReturnN());
    
    std::stringstream commentTopologicalDistance;
    commentTopologicalDistance << " File produced by analyzer UpdaterCreatorChainWalkingTertiaryBondWalking_CreationProperties" << std::endl
            << " Statistics of creation process involving the average Topological distance d along the structure" << std::endl
            << " between an attempted (and maybe successful) adding of monomer Dattempted" << std::endl
            << " between an performed adding of monomer Dperformed" << std::endl
            << " with " << ingredients.getMolecules().size() << " monomers" << std::endl
            << " Total number of samples: " << nrOfSamples << std::endl
            << std::endl
            << " p          ... attachment probability" << std::endl
            << " Dattempted ... average topological distance along the structure between attempted attachment" << std::endl
            << " Dperformed ... average topological distance along the structure between performed attachment" << std::endl
            << " SamplesDat ... total number of samples for Dattempted" << std::endl
            << " SamplesDpe ... total number of samples for Dperformed" << std::endl
            << "p Dattempted Dperformed samplesDat samplesDpe";

    //new filename
    std::string filenameTopologicalDistance = filenameGeneral + "_CreationPropertyTopologicalDistance.dat";

    ResultFormattingTools::writeResultFile(filenameTopologicalDistance, this->ingredients, tmpResultsTopologicalDistance, commentTopologicalDistance.str());

    /*******************/
    
    // construct a list
    std::vector < std::vector<double> > tmpResultsHG_TopologicalDistanceBetweenAttemptedInsertion;

    // we have 6 columns and row
    columns = 4;
    rows = 1;

    // we have columns
    tmpResultsHG_TopologicalDistanceBetweenAttemptedInsertion.resize(columns);

    // we have rows
    for (int i = 0; i < columns; i++)
        tmpResultsHG_TopologicalDistanceBetweenAttemptedInsertion[i].resize(0);

    // fill the histogram
    for (int bin = 0; bin < HG_TopologicalDistanceBetweenAttemptedInsertion->GetNrBins(); bin++) {
        if (HG_TopologicalDistanceBetweenAttemptedInsertion->GetNrInBin(bin) != 0) {
            tmpResultsHG_TopologicalDistanceBetweenAttemptedInsertion[0].push_back(HG_TopologicalDistanceBetweenAttemptedInsertion->GetRangeInBin(bin));
            tmpResultsHG_TopologicalDistanceBetweenAttemptedInsertion[1].push_back(HG_TopologicalDistanceBetweenAttemptedInsertion->GetNrInBin(bin));
            tmpResultsHG_TopologicalDistanceBetweenAttemptedInsertion[2].push_back(HG_TopologicalDistanceBetweenAttemptedInsertion->GetNrInBinNormiert(bin) / HG_TopologicalDistanceBetweenAttemptedInsertion->GetIntervallThickness());
            tmpResultsHG_TopologicalDistanceBetweenAttemptedInsertion[3].push_back(HG_TopologicalDistanceBetweenAttemptedInsertion->GetCumulativeNrInBinNormiert(bin) / HG_TopologicalDistanceBetweenAttemptedInsertion->GetIntervallThickness());
            
        }
    }
    
    std::stringstream commentHG_TopologicalDistanceBetweenAttemptedInsertion;
    commentHG_TopologicalDistanceBetweenAttemptedInsertion << " File produced by analyzer UpdaterCreatorChainWalkingTertiaryBondWalking_CreationProperties" << std::endl
            << " Statistics of creation process involving the average topological distance d along the structure" << std::endl
            << " between an attempted (and maybe successful) adding of monomer Dattempted" << std::endl
            << " between an performed adding of monomer Dperformed" << std::endl
            << " with " << ingredients.getMolecules().size() << " monomers" << std::endl
            << std::endl
            << " Total counts in all bins: " << HG_TopologicalDistanceBetweenAttemptedInsertion->GetNrCounts() << std::endl
            << " Total number of samples: " << nrOfSamples << std::endl
            << " Total number of bins " << HG_TopologicalDistanceBetweenAttemptedInsertion->GetNrBins() << std::endl
            << " Interval [" << HG_TopologicalDistanceBetweenAttemptedInsertion->GetRangeInBin(0) << " ; " << HG_TopologicalDistanceBetweenAttemptedInsertion->GetRangeInBin(HG_TopologicalDistanceBetweenAttemptedInsertion->GetNrBins()) << "]" << std::endl
            << " Interval thickness dI = " << HG_TopologicalDistanceBetweenAttemptedInsertion->GetIntervallThickness() << std::endl
            << " Attachment probability p: " << probabilityForInsertion << std::endl
            << " Average topological distance between attempted attachment Dattempted: " << TopologicalDistanceBetweenAttemptedInsertion.ReturnM1() << std::endl
            << " Average topological distance between successful attachment Dperformed: " << TopologicalDistanceBetweenPerformedInsertion.ReturnM1()  << std::endl
            << std::endl
            << " d ... topological distance between performed attachment" << std::endl
            << " <c(d)> ... total counts at topological distance (ensemble average)" << std::endl
            << " <cnorm(d)> ... average normalized counts by all counts at topological distance (ensemble sum)" << std::endl
            << " <snorm(d)> ... cumulative average normalized counts by all counts at topological distance d = sum_d <cnorm(d)> " << std::endl
            << "d <c(d)> <cnorm(d)> <snorm(d)>";

    //new filename
    std::string filenameHG_TopologicalDistanceBetweenAttemptedInsertion = filenameGeneral + "_CreationPropertyTopologicalDistanceHGAttempted.dat";

    ResultFormattingTools::writeResultFile(filenameHG_TopologicalDistanceBetweenAttemptedInsertion, this->ingredients, tmpResultsHG_TopologicalDistanceBetweenAttemptedInsertion, commentHG_TopologicalDistanceBetweenAttemptedInsertion.str());
    
    
    /********************************/
    
    // construct a list
    std::vector < std::vector<double> > tmpResultsHG_TopologicalDistanceBetweenPerformedInsertion;

    // we have 6 columns and row
    columns = 4;
    rows = 1;

    // we have columns
    tmpResultsHG_TopologicalDistanceBetweenPerformedInsertion.resize(columns);

    // we have rows
    for (int i = 0; i < columns; i++)
        tmpResultsHG_TopologicalDistanceBetweenPerformedInsertion[i].resize(0);

    // fill the histogram
    for (int bin = 0; bin < HG_TopologicalDistanceBetweenPerformedInsertion->GetNrBins(); bin++) {
        if (HG_TopologicalDistanceBetweenPerformedInsertion->GetNrInBin(bin) != 0) {
            tmpResultsHG_TopologicalDistanceBetweenPerformedInsertion[0].push_back(HG_TopologicalDistanceBetweenPerformedInsertion->GetRangeInBin(bin));
            tmpResultsHG_TopologicalDistanceBetweenPerformedInsertion[1].push_back(HG_TopologicalDistanceBetweenPerformedInsertion->GetNrInBin(bin));
            tmpResultsHG_TopologicalDistanceBetweenPerformedInsertion[2].push_back(HG_TopologicalDistanceBetweenPerformedInsertion->GetNrInBinNormiert(bin) / HG_TopologicalDistanceBetweenPerformedInsertion->GetIntervallThickness());
            tmpResultsHG_TopologicalDistanceBetweenPerformedInsertion[3].push_back(HG_TopologicalDistanceBetweenPerformedInsertion->GetCumulativeNrInBinNormiert(bin) / HG_TopologicalDistanceBetweenPerformedInsertion->GetIntervallThickness());
            
        }
    }
    
    std::stringstream commentHG_TopologicalDistanceBetweenPerformedInsertion;
    commentHG_TopologicalDistanceBetweenPerformedInsertion << " File produced by analyzer UpdaterCreatorChainWalkingTertiaryBondWalking_CreationProperties" << std::endl
            << " Statistics of creation process involving the average topological distance d along the structure" << std::endl
            << " between an attempted (and maybe successful) adding of monomer Dattempted" << std::endl
            << " between an performed adding of monomer Dperformed" << std::endl
            << " with " << ingredients.getMolecules().size() << " monomers" << std::endl
            << std::endl
            << " Total counts in all bins: " << HG_TopologicalDistanceBetweenPerformedInsertion->GetNrCounts() << std::endl
            << " Total number of samples: " << nrOfSamples << std::endl
            << " Total number of bins " << HG_TopologicalDistanceBetweenPerformedInsertion->GetNrBins() << std::endl
            << " Interval [" << HG_TopologicalDistanceBetweenPerformedInsertion->GetRangeInBin(0) << " ; " << HG_TopologicalDistanceBetweenPerformedInsertion->GetRangeInBin(HG_TopologicalDistanceBetweenPerformedInsertion->GetNrBins()) << "]" << std::endl
            << " Interval thickness dI = " << HG_TopologicalDistanceBetweenPerformedInsertion->GetIntervallThickness() << std::endl
            << " Attachment probability p: " << probabilityForInsertion << std::endl
            << " Average topological distance between attempted attachment Dattempted: " << TopologicalDistanceBetweenAttemptedInsertion.ReturnM1() << std::endl
            << " Average topological distance between successful attachment Dperformed: " << TopologicalDistanceBetweenPerformedInsertion.ReturnM1()  << std::endl
            << std::endl
            << " d ... topological distance between performed attachment" << std::endl
            << " <c(d)> ... total counts at topological distance (ensemble average)" << std::endl
            << " <cnorm(d)> ... average normalized counts by all counts at topological distance (ensemble sum)" << std::endl
            << " <snorm(d)> ... cumulative average normalized counts by all counts at topological distance d = sum_d <cnorm(d)> " << std::endl
            << "d <c(d)> <cnorm(d)> <snorm(d)>";

    //new filename
    std::string filenameHG_TopologicalDistanceBetweenPerformedInsertion = filenameGeneral + "_CreationPropertyTopologicalDistanceHGPerformed.dat";

    ResultFormattingTools::writeResultFile(filenameHG_TopologicalDistanceBetweenPerformedInsertion, this->ingredients, tmpResultsHG_TopologicalDistanceBetweenPerformedInsertion, commentHG_TopologicalDistanceBetweenPerformedInsertion.str());
    
    
    /********************************/
}


#endif /* LEMONADE_UPDATER_CREATOR_CHAIN_WALKING_TertiaryBondWalking_H */
