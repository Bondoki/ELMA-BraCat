/*--------------------------------------------------------------------------------
    ooo      L   attice-based  | LeMonADE: An Open Source Implementation of the
  o\.|./o    e   xtensible     |           Bond-Fluctuation-Model for Polymers
 o\.\|/./o   Mon te-Carlo      |           
oo---0---oo  A   lgorithm and  | ELMA-BraCat: Extension for branched catalysis 
 o/./|\.\o   D   evelopment    | Copyright (C) 2013-2015, 2018,2021,2022 by
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

#ifndef LEMONADE_UPDATER_CREATE_DENDRIMER_FGS_H
#define LEMONADE_UPDATER_CREATE_DENDRIMER_FGS_H
/**
 * @file
 *
 * @class UpdaterCreateDendrimerGeneral_FGS
 *
 * @brief General Updater for a single dendrimer with functionality F, generation G and spacer length S
 * 
 * @details A single dendrimer with arbitrary generation, spacer length and functionality 
 * is created in a cubix box of arbitrary size togehter with linear chains of arbitrary 
 * length with arbitrary concentration
 *
 * @tparam IngredientsType
 *
 **/

#include <LeMonADE/updater/UpdaterAbstractCreate.h>
#include <LeMonADE/utility/Vector3D.h>


#include "StatisticMoment.h"
#include "Histogram1D.h"

template<class IngredientsType>
class UpdaterCreateDendrimerGeneral_FGS: public UpdaterAbstractCreate<IngredientsType>
{
  typedef UpdaterAbstractCreate<IngredientsType> BaseClass;
  
public:
  UpdaterCreateDendrimerGeneral_FGS(IngredientsType& ingredients_, uint32_t generation_, uint32_t spacerlength_, uint32_t functionality_, uint32_t box_x, uint32_t box_y, uint32_t box_z);
  
  virtual void initialize();
  virtual bool execute();
  virtual void cleanup();
  
private:
  using BaseClass::ingredients;
  
  using BaseClass::addMonomerToParent;
  using BaseClass::addSingleMonomer;
  using BaseClass::addMonomerAtPosition;
  using BaseClass::addMonomerInsideConnectedPair;
  using BaseClass::linearizeSystem;
  using BaseClass::moveSystem;

  //! Generation of Dendrimer
  uint32_t generation;
  
  //! simulation box sizes
  uint32_t boxX,boxY,boxZ;
  
  //! Spacer length
  uint32_t spacerlength;
  
  //! Functionality
  uint32_t functionality;

  //! provides a histogram of functionality 
  Histogram1D* HG_FunctionalityNodes; 
};

/** 
* @brief Constructor handling the new systems paramters
*
* @param ingredients_ a reference to the IngredientsType - mainly the system
* @param L_ linear chain length
* @param box_ size of the box
*/
template < class IngredientsType >
UpdaterCreateDendrimerGeneral_FGS<IngredientsType>::UpdaterCreateDendrimerGeneral_FGS(IngredientsType& ingredients_, uint32_t generation_, uint32_t spacerlength_, uint32_t functionality_, uint32_t box_x, uint32_t box_y, uint32_t box_z):
BaseClass(ingredients_), generation(generation_), spacerlength(spacerlength_), functionality(functionality_), boxX(box_x), boxY(box_y), boxZ(box_z)
{
    HG_FunctionalityNodes = new Histogram1D(0 - 0.5, 2000 - 0.5, 2000);
        
}

/**
* The initialize function handles the new systems information.
*
* @tparam IngredientsType Features used in the system. See Ingredients.
*/
template < class IngredientsType >
void UpdaterCreateDendrimerGeneral_FGS<IngredientsType>::initialize(){
  std::cout << "initialize UpdaterCreateDendrimerGeneral_FGS" << std::endl;
  
  //set box size
  ingredients.setBoxX(boxX);
  ingredients.setBoxY(boxY);
  ingredients.setBoxZ(boxZ);
  
  //set periodicity
  ingredients.setPeriodicX(false);//true);
  ingredients.setPeriodicY(false);//true);
  ingredients.setPeriodicZ(false);
  
  ingredients.modifyMolecules().setAge(0);
  //add Bondset
  // supress std::cout output of addBondset
  std::streambuf *old = std::cout.rdbuf(); // <-- save
  std::stringstream ss;
  std::cout.rdbuf (ss.rdbuf());       // <-- redirect
  ingredients.modifyBondset().addBFMclassicBondset();
  std::cout.rdbuf (old);
  
  ingredients.synchronize(ingredients);
  
  

 // execute();
}

/**
* Execution of the system creation
*
* @tparam IngredientsType Features used in the system. See Ingredients.
*/
template < class IngredientsType >
bool UpdaterCreateDendrimerGeneral_FGS<IngredientsType>::execute(){
  std::cout << "execute UpdaterCreateDendrimerGeneral_FGS" << std::endl;
  
  //check if system was already created
  if(ingredients.getMolecules().size()!=0)
    return true;
  
  // create first monomer of the chain
  addMonomerAtPosition(VectorInt3(boxX/2,boxY/2,boxZ/2),1) ? : throw std::runtime_error("UpdaterCreateDendrimerGeneral_FGS: addSingleMonomer is not able to place a monomer!");
  //ingredients.modifyMolecules()[ingredients.getMolecules().size()-1].setMovableTag(false);
  ingredients.modifyMolecules()[ingredients.getMolecules().size()-1].setAttributeTag(1);

  int idxCenter = ingredients.getMolecules().size()-1;

  ingredients.synchronize();

  std::cout << "g = "<<generation <<", s = "<< spacerlength<< ", functionality" << functionality << ", box = "<< boxX << ", " << boxY << ", " << boxZ<< std::endl;

    //look for size of generation, if smaller than 1 set to one and generate smallest dendrimer
  if(generation<1){throw std::runtime_error("ERROR: number of generations smaller than 1");}

  std::vector<int32_t> branches;


  //generate generation one ad hoc

    //add all branches to center
    for(uint32_t branch = 0; branch < functionality; branch++)
    {
    	for(uint32_t i=0;i<spacerlength;i++){

    		// add at center
    		if(i==0){
    			// tag between 1 and 2
    			int idxParent = idxCenter;
    			int tag = (ingredients.getMolecules()[idxParent].getAttributeTag()%2)+1;

    			while( !addMonomerToParent(idxParent, tag))
    			{
    				moveSystem(1);
    			}
    		}
    		else // add at spacer
    		{
    			// tag between 1 and 2
    			int idxParent = ingredients.getMolecules().size()-1;
    			int tag = (ingredients.getMolecules()[idxParent].getAttributeTag()%2)+1;

    			while( !addMonomerToParent(idxParent, tag))
    			{
    				moveSystem(1);
    			}
    		}

    	}
    	branches.push_back(ingredients.getMolecules().size()-1);
    }

  //set up system, occupy lattice
  ingredients.synchronize();


  //build up system for g>1 dendrimers FC=F=functionality
  if(generation>1){
	  //loop over number of generations
	  std::vector<int32_t>dummy;
	  for(int32_t g=1;g<generation;g++){
		  //loop over all entries in branches vector
		  for(int32_t i=0;i<branches.size();i++){
			  //cout << branches.at(i)<<std::endl;

			  //loop over all remaining sidebranches
			  for(int32_t functionality=0; functionality < (this->functionality-1); functionality++)
			  {
			  //loop over spacerlength the first time (first new branch)
			  for(int32_t s=0;s<spacerlength;s++){
				  //get index of parent monomer
				  int32_t idxParent;
				  if(s==0){
					  idxParent=branches.at(i);
				  }else{
					  idxParent=ingredients.getMolecules().size()-1;
				  }
				  int tag = (ingredients.getMolecules()[idxParent].getAttributeTag()%2)+1;

				  while( !addMonomerToParent(idxParent, tag))
				  {
				  	moveSystem(1);
				  }
				  std::cout << "added monomer :" <<  (ingredients.getMolecules().size()-1)  << std::endl;
			  }
			  //adding new branching point
			  dummy.push_back(ingredients.getMolecules().size()-1);

			  }
		  }
		  branches=dummy;
		  dummy.clear();
	  }
  }

/*
  for(int z = 0; z < ingredients.getMolecules().size(); z++)
  {
	  ingredients.modifyMolecules()[z].setAttributeTag(0);
  }

  for(int z = 0; z < ingredients.getMolecules().size(); z++)
    {
	  if(ingredients.getMolecules().getNumLinks(z) == 1)
  	  ingredients.modifyMolecules()[z].setAttributeTag(2);
    }
*/
/*
  // adding cross-linker to the system
  for(int cl=0; cl < crosslinker; cl++)
  {
	  // Cross-linker are B-Types
	  while( !addSingleMonomer(1))
	  {
		  moveSystem(1);
	  }
	  std::cout << "added cross-linker :" <<  (ingredients.getMolecules().size()-1)  << std::endl;

  }
*/

  std::cout << "Creation done with " << ingredients.getMolecules().size() << " monomers" << std::endl;

  std::cout << "Finial sync...";
  ingredients.synchronize();
  std::cout << "done." << std::endl;

  std::cout << "Finial linearization...";

  linearizeSystem();


  std::cout << "done." << std::endl;

  // create the histogram of occurence of functionality in structure    
        for (int idx = 0; idx < ingredients.getMolecules().size(); idx++)
        {
            HG_FunctionalityNodes->AddValue(ingredients.getMolecules().getNumLinks(idx));
        }
  
  
   std::cout << "fill histogram done" << std::endl;
   return true;
}

/**
* Standard clean up.
*
* @tparam IngredientsType Features used in the system. See Ingredients.
*/
template < class IngredientsType >
void UpdaterCreateDendrimerGeneral_FGS<IngredientsType>::cleanup(){
  
    std::cout << "file out" << std::endl;
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
    commentHG_FunctionalityNodes << " File produced by analyzer UpdaterCreateDendrimerGeneral_FGS" << std::endl
            << " Statistics of dendrimer structure providing the occurence of functionality f" << std::endl
            << " with " << ingredients.getMolecules().size() << " monomers" << std::endl
            << std::endl
            << " Total counts in all bins: " << HG_FunctionalityNodes->GetNrCounts() << std::endl
            << " Total number of samples: " << 1 << std::endl
            << " Total number of bins " << HG_FunctionalityNodes->GetNrBins() << std::endl
            << " Interval [" << HG_FunctionalityNodes->GetRangeInBin(0) << " ; " << HG_FunctionalityNodes->GetRangeInBin(HG_FunctionalityNodes->GetNrBins()) << "]" << std::endl
            << " Interval thickness dI = " << HG_FunctionalityNodes->GetIntervallThickness() << std::endl
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


#endif /* LEMONADE_UPDATER_CREATE_SINGLE_DENDRIMER_IN_LINEAR_MELT_H */
