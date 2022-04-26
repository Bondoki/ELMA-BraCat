/*--------------------------------------------------------------------------------
    ooo      L   attice-based  |
  o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
oo---0---oo  A   lgorithm and  |
 o/./|\.\o   D   evelopment    | Copyright (C) 2018,2022 by
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

#ifndef Analyzer_TopologicalPropertiesTree_H
#define Analyzer_TopologicalPropertiesTree_H

#include <vector>
#include <string>
#include <utility>      // std::pair
#include <map>
#include <vector>
#include <algorithm>    // std::find_if
#include <iostream>
#include <functional>
#include <queue>

#include <LeMonADE/utility/Vector3D.h>

#include "StatisticMoment.h"
#include "HistogramGeneralStatistik1D.h"



template<class IngredientsType>
class Analyzer_TopologicalPropertiesTree:public AbstractAnalyzer
{
public:
	Analyzer_TopologicalPropertiesTree(const IngredientsType& ing, std::string dstDir_);


	virtual ~Analyzer_TopologicalPropertiesTree(){

	};

	//typedef typename IngredientsType::molecules_type molecules_type;
	const typename IngredientsType::molecules_type& molecules ;

	const IngredientsType& getIngredients() const {return ingredients;}

	virtual void initialize();
	virtual bool execute();
	virtual void cleanup();

	typedef uint32_t NodeIdx;
	typedef std::map<NodeIdx, int> Nodelist;
	typedef std::map<NodeIdx, Nodelist> Graph;
	  //typedef std::pair<NodeIdx, Nodelist> Node;
	  typedef std::pair<int, NodeIdx> Edge; //! (tentative distance, idxNode)
	  typedef std::vector<NodeIdx> NodeVector;

	  typedef std::map<NodeIdx, NodeVector> NodePath;

	  void dijkstra(Graph &graph, NodeIdx source, Nodelist &distance);

	  void dijkstraWithPath(Graph &graph, NodeIdx source, Nodelist &distance, NodePath &distanceNodes);

	  void getPath(NodeIdx startNode, NodeIdx endNode, NodePath predecessor, NodeVector &path);

	  void findCentersOfTree(Graph &graph,  Nodelist degreeNodes, NodeVector &centerNodes);

private:

	const IngredientsType& ingredients;

	std::map<uint32_t, uint64_t> frameToTime;
	std::map<uint64_t, std::vector<VectorInt3> > timeToPosChainEnd; //end monomer
	std::map<uint64_t, std::vector<VectorInt3> > timeToPosChainStart; // start monomer
	std::map<uint64_t, std::vector<VectorInt3> > timeToPosChainHalf; // middle monomer

	std::map<uint64_t, std::vector<VectorDouble3> > timeToPosCOMChain; // COM Chain
	std::map<uint64_t, VectorDouble3 > timeToPosCOMBox; // COM Box


	uint32_t counterFrames;

	uint32_t numMonoPerObjects;
	uint32_t numObjects;
	uint32_t idxMiddleMonomer;
	uint32_t idxTerminalMonomer;

	uint32_t startFrame;
	uint32_t diffFrame;


	std::string filename;
	std::string dstdir;

	 HistogramGeneralStatistik1D HG_Statistik_TopoDistance_Walking;
	// RNG
	RandomNumberGenerators rng;

	// graph properties

	NodeVector centerNodes;// center of the tree (either one or two)
	Nodelist degreeNode; // Degree of node == functionality
	Nodelist topo_distance; // Topological distance from start node to desired node
};




/////////////////////////////////////////////////////////////////////////////

template<class IngredientsType>
Analyzer_TopologicalPropertiesTree<IngredientsType>::Analyzer_TopologicalPropertiesTree(const IngredientsType& ing, std::string dstDir_)
:ingredients(ing), molecules(ing.getMolecules()), dstdir(dstDir_)
 {
	counterFrames = 1;

 }

template<class IngredientsType>
void Analyzer_TopologicalPropertiesTree<IngredientsType>::initialize()
{
	std::cout << "init histogram" <<std::endl;


	HG_Statistik_TopoDistance_Walking.reset(0.0, 10000, 10000);




	execute();

}

template<class IngredientsType>
bool Analyzer_TopologicalPropertiesTree<IngredientsType>::execute()
{
	Graph treeGraph;


		//Nodelist  degreeNode; // Degree of node

		for (int k= 0; k < ingredients.getMolecules().size(); k++)
		{
			for (int l = 0; l < ingredients.getMolecules().getNumLinks(k); l++)
			{
				// connect from idxAA to idxBB with value ZZ
				// all (symmetric) connections
				treeGraph[k][ingredients.getMolecules().getNeighborIdx(k, l)]=1;

			}
			// fill with degree of nodes
			degreeNode[k]=ingredients.getMolecules().getNumLinks(k);
		}


		// calculate eccentricity of graph for all nodes
		std::cout << "Calculate eccentricity" << std::endl;
		std::map<NodeIdx, int> eccentricity;
		std::map<NodeIdx, Nodelist> distance;

		for (int k= 0; k < ingredients.getMolecules().size(); k++)
		{
			NodeIdx startNode = k;
			dijkstra(treeGraph, startNode, topo_distance);

			int maxTopologicalDistance = 0;

			for(Nodelist::iterator it=topo_distance.begin(); it!=topo_distance.end(); ++it)
			{
				//std::cout<< startNode << " -> " << it->first << "  => "<<it->second<<std::endl;
				distance[startNode][it->first] = it->second;

				if(maxTopologicalDistance < it->second)
					maxTopologicalDistance=it->second;
			}

			eccentricity[startNode]=maxTopologicalDistance;
			std::cout<< startNode << " -> eccentricity: " << eccentricity[startNode] <<std::endl;
		}

		// calculate the radius of the graph = min eccentricity
		int radiusGraph = std::numeric_limits<int>::max();

		for (std::map<NodeIdx, int>::iterator it_ecc=eccentricity.begin(); it_ecc!=eccentricity.end(); ++it_ecc)
		{
			if(it_ecc->second < radiusGraph)
				radiusGraph = it_ecc->second;
		}
		std::cout << "Radius Graph = " << radiusGraph <<std::endl;

		// calculate center nodes: eccentricity[center] == radiusGraph
		NodeVector centerGraph;
		for (std::map<NodeIdx, int>::iterator it_ecc=eccentricity.begin(); it_ecc!=eccentricity.end(); ++it_ecc)
		{
			if(it_ecc->second == radiusGraph)
				centerGraph.push_back(it_ecc->first);
		}

		std::cout << "Center Graph = ( ";
		for(NodeVector::iterator it_center=centerGraph.begin(); it_center!=centerGraph.end(); ++it_center){
			std::cout << *it_center << " : ";
		}
		std::cout << " )" << std::endl;

		// calculate the diameter of the graph = max eccentricity
		int diameterGraph = std::numeric_limits<int>::min();

		for (std::map<NodeIdx, int>::iterator it_ecc=eccentricity.begin(); it_ecc!=eccentricity.end(); ++it_ecc)
		{
			if(diameterGraph < it_ecc->second)
				diameterGraph = it_ecc->second;
		}
		std::cout << "Diameter Graph = " << diameterGraph <<std::endl;

		// calculate periphery nodes: eccentricity[center] == diameterGraph
		NodeVector peripheryGraph;
		for (std::map<NodeIdx, int>::iterator it_ecc=eccentricity.begin(); it_ecc!=eccentricity.end(); ++it_ecc)
		{
			if(it_ecc->second == diameterGraph)
				peripheryGraph.push_back(it_ecc->first);
		}

		std::cout << "Periphery Graph = ( ";
		for(NodeVector::iterator it_periphery=peripheryGraph.begin(); it_periphery!=peripheryGraph.end(); ++it_periphery){
			std::cout << *it_periphery << " : ";
		}
		std::cout << " )" << std::endl;

		// print path between nodes
		NodeIdx startNode = peripheryGraph[0];
		Nodelist topologicalDistance;
		NodePath PathToNode;

		dijkstraWithPath(treeGraph, startNode, topologicalDistance, PathToNode);

		std::cout << "startNode -> node => distance" <<std::endl;

		for(Nodelist::iterator it=topologicalDistance.begin(); it!=topologicalDistance.end(); ++it)
		{
			std::cout<< "From " << startNode << " -> " << it->first << "   with distance => "<<it->second<<std::endl;

			NodeIdx startRecursive = it->first;
			NodeVector PathFromNode;

			do
			{
				PathFromNode.push_back(startRecursive);
				startRecursive = PathToNode[ startRecursive ].at(0);

			}while(startRecursive != startNode);

			// add destiny node
			PathFromNode.push_back(startNode);

			// print path
			std::cout << "Path: ( ";

			for(NodeVector::iterator it_nv=PathFromNode.begin(); it_nv!=PathFromNode.end(); ++it_nv)
			{
				std::cout << *it_nv << ":" ;
			}
			std::cout<<" )" << std::endl;
		}

		// backbone or longest path on the tree
		// start from periphery and search for one longest distance
		//
		// every tree has at least one leave with highest distance (eccentricity)
		NodeIdx peripheryNode = peripheryGraph[0];
		Nodelist peripheryDistance = distance[peripheryNode];

		int maxLongestPathDistance = 0;
		NodeIdx maxLongestPathDistanceNode = 0;

		std::cout << std::endl << "Backbone (longest path) " << std::endl;

		for(Nodelist::iterator it=peripheryDistance.begin(); it!=peripheryDistance.end(); ++it)
		{
			std::cout<< peripheryNode << " -> " << it->first << "  => "<<it->second<<std::endl;

			if(maxLongestPathDistance < it->second)
			{
				maxLongestPathDistance=it->second;
				maxLongestPathDistanceNode=it->first;
			}
		}

		std::cout << "Backbone start node " <<  peripheryNode << " -> end node " << maxLongestPathDistanceNode << std::endl;



		// print path along backbone
		Nodelist topologicalDistanceBackbone;
		NodePath PathPredecessor;
		NodeVector PathOnBackbone;

		dijkstraWithPath(treeGraph, peripheryNode, topologicalDistanceBackbone, PathPredecessor);

		std::cout << "Backbone path: peripheryNode -> node => distance" <<std::endl;

		getPath(peripheryNode, maxLongestPathDistanceNode, PathPredecessor, PathOnBackbone);

		std::cout << std::endl;
		for(NodeVector::iterator it_nv=PathOnBackbone.begin(); it_nv!=PathOnBackbone.end(); ++it_nv)
		{
		std::cout << "!setColor:" << *it_nv << "-"<< *it_nv <<"=(1,1,1);" << std::endl;
		}
		std::cout<< std::endl;

		//finding the center of the tree (either one or two)
		//NodeVector centerNodes;
		findCentersOfTree(treeGraph, degreeNode, centerNodes);

		//calculate the topological distance of center to all nodes and fill histogram
		NodeVector::iterator itv;
		for(itv=centerNodes.begin(); itv!=centerNodes.end(); ++itv){



			NodeIdx startNode = (*itv);
			std::cout << std::endl << "startNode: " << startNode << std::endl;

			//Nodelist topo_distance;
			dijkstra(treeGraph, startNode, topo_distance);

			std::cout << "startNode -> node => distance" <<std::endl;

			Nodelist::iterator it;
			for(it=topo_distance.begin(); it!=topo_distance.end(); ++it){

				std::cout<< startNode << " -> " << it->first << "  => "<<it->second<<std::endl;


				//HG_Statistik_TopoDistance_zMean (TopoDistance, Z)
			//	HG_Statistik_TopoDistance_zMean.addValue(it->second, ingredients.getMolecules()[it->first].getZ());
				//HGStatistikRg2_XY.addValue(COM_Z, Rg2_x+Rg2_y);
			//	HG_TopoDistance.AddValue(it->second);
				/*
			       	    // add value to the histogram
			       	    HG_TopoDistance->AddValue(it->second);

		                    // add value to the histogram for number links (functionality)
		                    for(int o=0; o < degreeNode[it->first]; o++)
		                    {
		                    HG_FunctionalityTopoDistance->AddValue(it->second);
		                    }
				 */
			}
		}

		{
		//calculate the topological distance of center to all nodes and fill histogram
		NodeVector::iterator itv;
		for(itv=centerNodes.begin(); itv!=centerNodes.end(); ++itv){



			NodeIdx startNode = (*itv);
			std::cout << std::endl << "startNode: " << startNode << std::endl;

			NodePath PathToNode;
			//Nodelist topo_distance;
			dijkstraWithPath(treeGraph, startNode, topo_distance, PathToNode);

			std::cout << "startNode -> node => distance" <<std::endl;

			Nodelist::iterator it;
			for(it=topo_distance.begin(); it!=topo_distance.end(); ++it){

				std::cout<< "From " << startNode << " -> " << it->first << "   with distance => "<<it->second<<std::endl;

				NodeIdx startRecursive = it->first;
				NodeVector PathFromNode;

				do
				{
					PathFromNode.push_back(startRecursive);
					startRecursive = PathToNode[ startRecursive ].at(0);

				}while(startRecursive != startNode);

				// add destiny node
				PathFromNode.push_back(startNode);

				// print path
				std::cout << "Path: ( ";

				for(NodeVector::iterator it_nv=PathFromNode.begin(); it_nv!=PathFromNode.end(); ++it_nv)
				{
					std::cout << *it_nv << ":" ;
				}
				std::cout<<" )" << std::endl;


				NodeVector nv = PathToNode[it->first ];
				std::cout<< "size nv: " << nv.size()<<std::endl;
				NodeVector::iterator it_nv;
							for(it_nv=nv.begin(); it_nv!=nv.end(); ++it_nv){

								std::cout<< startNode << " -> " << *it_nv<<std::endl;
							}
				std::cout<<std::endl;

				//HG_Statistik_TopoDistance_zMean (TopoDistance, Z)
			//	HG_Statistik_TopoDistance_zMean.addValue(it->second, ingredients.getMolecules()[it->first].getZ());
				//HGStatistikRg2_XY.addValue(COM_Z, Rg2_x+Rg2_y);
			//	HG_TopoDistance.AddValue(it->second);
				/*
			       	    // add value to the histogram
			       	    HG_TopoDistance->AddValue(it->second);

		                    // add value to the histogram for number links (functionality)
		                    for(int o=0; o < degreeNode[it->first]; o++)
		                    {
		                    HG_FunctionalityTopoDistance->AddValue(it->second);
		                    }
				 */
			}
		}
		}

/*
	NodeIdx startNode = *(--centerNodes.end());
	NodeIdx rwOnNode = startNode;

	std::cout << "Start random walk at node: startNode == " << startNode <<std::endl;
	std::cout << "with functionality : f == " << degreeNode[startNode] <<std::endl;

	std::cout << std::endl << "Place random walker and run simulation" << std::endl;

	for(int rw_steps = 0; rw_steps < 10000; rw_steps++)
	{
		//std::cout<< " RW Step: " << rw_steps << " reside on node " << rwOnNode << " with distance -> " << topo_distance[rwOnNode] << std::endl;

		// calculating MSD on structure
		HG_Statistik_TopoDistance_Walking.addValue(rw_steps, topo_distance[rwOnNode]*topo_distance[rwOnNode]);

		// move the walker
		uint32_t possible_directions = degreeNode[rwOnNode];

		// RW-node has neighbors also select one between [0; f)
		uint32_t rand_dir = rng.r250_rand32()%possible_directions;

		rwOnNode = ingredients.getMolecules().getNeighborIdx(rwOnNode, rand_dir);

	}
	*/


	counterFrames++;
	std::cout << " run : " << counterFrames << std::endl;

	return true;
}


template<class IngredientsType>
void Analyzer_TopologicalPropertiesTree<IngredientsType>::findCentersOfTree(Graph &graph, Nodelist degreeNodes, NodeVector &centerNodes){

	std::queue <NodeIdx> queueLeaves; // Queue for algorithm holding (reduced) leaves
	std::queue <NodeIdx> queueNextLeaves; // Queue for algorithm holding next (reduced) leaves
	std::queue <NodeIdx> tmpOldLeaves; // Queue for algorithm holding (reduced) leaves

	// Fill and start from leaves

	Nodelist::iterator it;
	for(it=degreeNodes.begin(); it!=degreeNodes.end(); ++it)
	{
	        int degreeOfNode = it->second;
	        NodeIdx labelNode = it->first;

	        if (degreeOfNode == 1)
	        {
	        	queueLeaves.push(labelNode);
	        }
	}

	NodeIdx Seperator = -1;

	queueLeaves.push(Seperator);

	//  while ( (queueLeaves.front()!=Seperator)||(queueLeaves.size() > 2) ) {
	while ( (queueLeaves.size() > 0) ) {
		  NodeIdx node = queueLeaves.front();
		  queueLeaves.pop();

		  if( node == Seperator)
		  {
			  if(queueNextLeaves.empty())
				  break;

			 // std::cout << "next iteration: " << node << std::endl;
			  tmpOldLeaves = queueNextLeaves;
			  queueLeaves = queueNextLeaves;

			  if(queueNextLeaves.size() >0)
				  queueLeaves.push(Seperator);

			  queueNextLeaves=std::queue <NodeIdx>();

			  continue;
		  }


	 	      int dN_old = degreeNodes[node];

	 	      degreeNodes[node]--;

	 	     // new subgraph of all neighbors
	 	     Nodelist tempgraph=graph[node];

	 	     // connected nodes

	 	    // std::cout << "connected nodes from: " << node << std::endl;

	 	     Nodelist::iterator it;
	 	     for(it=tempgraph.begin(); it!=tempgraph.end(); ++it){

	 	    //	std::cout << it->first << " with degree: " << degreeNodes[it->first];

	 	    	 if(degreeNodes[it->first] > 0)
	 	    	 {
	 	    		degreeNodes[it->first]--;

	 	    	//	std::cout << " reduced to degree: " << degreeNodes[it->first];
	 	    	 }

	 	    	 if(degreeNodes[it->first] == 1)
	 	    	 {
	 	    		//queueLeaves.push(it->first);
	 	    		queueNextLeaves.push(it->first);

	 	    	//	std::cout << " added to queue the idx: " << (it->first);
	 	    	 }

	 	    	//std::cout << std::endl;
	 	     }


	 	  }

	std::cout << "center queue contains: ";
		  while (!tmpOldLeaves.empty())
		  {
			  centerNodes.push_back(tmpOldLeaves.front());
			  std::cout << ' ' <<  tmpOldLeaves.front();
			  tmpOldLeaves.pop();
		  }
		  std::cout << '\n';
}

template<class IngredientsType>
void Analyzer_TopologicalPropertiesTree<IngredientsType>::getPath(NodeIdx startNode, NodeIdx endNode, NodePath predecessor, NodeVector &path){

	NodeIdx startRecursive = endNode;

	// add start node
	do
	{
		path.push_back(startRecursive);
		startRecursive = predecessor[ startRecursive ].at(0);

	}while(startRecursive != startNode);

	// add destiny node
	path.push_back(startNode);

	// print path
	std::cout << "Path: ( ";

	for(NodeVector::iterator it_nv=path.begin(); it_nv!=path.end(); ++it_nv)
	{
		std::cout << *it_nv << ":" ;
	}
	std::cout<<" )" << std::endl;
}

template<class IngredientsType>
void Analyzer_TopologicalPropertiesTree<IngredientsType>::dijkstra(Graph &graph, NodeIdx source, Nodelist &distance){

	distance.clear(); //! clear all tentative distance information

	//! list of all nodes to be visit and addressed already with least distance on top
  std::priority_queue<Edge, std::vector<Edge>, std::greater<Edge> > queueNode;

  // starting node with tentative distance=0, starting node = source
  queueNode.push( Edge(0, source) );

  while(!queueNode.empty()){

	  //get the element with least tentative distance
    Edge tmped=queueNode.top();

    // access the node index
    NodeIdx tmpnl=tmped.second;

    // removes the top element
    queueNode.pop();

    // if the node never visited before
    if(distance.count(tmpnl)==0){

    	// tentative distance to the recent node
      int dist=tmped.first;

      // set the tentative distance to the node
      distance[tmpnl]=dist;

      // new subgraph of all neighbors
      Nodelist tempgraph=graph[tmpnl];

      Nodelist::iterator it;
      for(it=tempgraph.begin(); it!=tempgraph.end(); ++it){
        int distint=it->second;
        NodeIdx distlabel=it->first;
        queueNode.push(Edge(dist+distint, distlabel));
      }

    }
  }

}

template<class IngredientsType>
void Analyzer_TopologicalPropertiesTree<IngredientsType>::dijkstraWithPath(Graph &graph, NodeIdx source, Nodelist &distance, NodePath &path){

	distance.clear(); //! clear all tentative distance information

	typedef std::pair<int, std::pair<NodeIdx, NodeIdx>> EdgeWithPrevious; //! (tentative distance, idxNode)
	//! list of all nodes to be visit and addressed already with least distance on top
  std::priority_queue<EdgeWithPrevious, std::vector<EdgeWithPrevious>, std::greater<EdgeWithPrevious> > queueNode;

  // starting node with tentative distance=0, starting node = source
  queueNode.push( EdgeWithPrevious(0, std::make_pair(source, source)) );

  while(!queueNode.empty()){

	  //get the element with least tentative distance
	  EdgeWithPrevious tmped=queueNode.top();

    // access the node index
    NodeIdx tmpnl=tmped.second.first;

    NodeIdx tmpnl_prev=tmped.second.second;

    // removes the top element
    queueNode.pop();

    // if the node never visited before
    if(distance.count(tmpnl)==0){

    	// tentative distance to the recent node
      int dist=tmped.first;

      // set the tentative distance to the node
      distance[tmpnl]=dist;

      //push path to node
      path[tmpnl].push_back(tmpnl_prev);

      // new subgraph of all neighbors
      Nodelist tempgraph=graph[tmpnl];

      Nodelist::iterator it;
      for(it=tempgraph.begin(); it!=tempgraph.end(); ++it){
        int distint=it->second;
        NodeIdx distlabel=it->first;
        queueNode.push(EdgeWithPrevious(dist+distint, std::make_pair(distlabel, tmpnl)));


      }

    }
  }

}


struct PathSeparator
{
    bool operator()( char ch ) const
    {
        return ch == '\\' || ch == '/';
    }
};

template<class IngredientsType>
void Analyzer_TopologicalPropertiesTree<IngredientsType>::cleanup()
{
/*	std::cout << "Calculation of the MSD" << std::endl;

	std::map<uint32_t, uint64_t>::reverse_iterator it=frameToTime.rbegin();

	std::cout << "Max Frame is " << frameToTime.size() << " at " << it->first << " -> time " << uint64_t(it->second) << std::endl;

	std::map<uint64_t, std::vector<VectorInt3> >::reverse_iterator itTime=timeToPosChainStart.rbegin();
	std::cout << "Max Time at " << timeToPosChainStart.size() << " at time " << itTime->first << std::endl;

	uint32_t maxframe = frameToTime.size();
	//std::cout << "Max Frame is " << frameToTime.size() << '\n'
	//calculate g1, g2, g3 by time reversal


	std::cout << "mapping contains:\n";
			  for ( std::map<uint32_t, uint64_t>::iterator it=frameToTime.begin(); it!=frameToTime.end(); ++it)
			  {
			    std::cout << "Frame: " << it->first << " => Time: " << it->second << '\n';
			  }

	std::cout << std::endl << "Calculation of the MSD" << std::endl;

	std::map<uint64_t, StatisticMoment > g1; // MSD middle monomer
	std::map<uint64_t, StatisticMoment > g2; // MSD middle monomer - COM chain
	std::map<uint64_t, StatisticMoment > g3; // MSD COM chain
	std::map<uint64_t, StatisticMoment > g4; // MSD start and end monomer

	std::cout << "Diff Frame: ";
	// statistical overcounting of "independent" configurations
	for(int32_t dframe = 0; dframe <= int32_t(maxframe-startFrame); dframe+=diffFrame)
	{
		std::cout << " " << dframe << " ";
		// frame numbering start at 1 to incl. maxframe
		for(uint32_t frame = 0; frame < maxframe; frame++)
		{
			// time reversal
			if( (maxframe-frame-dframe) >= startFrame)
			{
				// find time to frame
				std::map<uint32_t, uint64_t>::iterator itFrameToTime;

				itFrameToTime = frameToTime.find(maxframe-frame-dframe);

				if (itFrameToTime != frameToTime.end())
				{
					uint64_t timeForEval = itFrameToTime->second;
					//std::cout << itFrameToTime->first << " => " << itFrameToTime->second << '\n';

					std::map<uint32_t, uint64_t>::iterator itFrameToTimeRef;
					itFrameToTimeRef = frameToTime.find(maxframe-dframe);

					if (itFrameToTimeRef != frameToTime.end())
					{
						uint64_t timeForEvalRef = itFrameToTimeRef->second;

						//std::cout << "Eval at: " << timeForEval << " => Ref " << timeForEvalRef << '\n';

						for (int nrObjects= 0; nrObjects < numObjects; nrObjects++)
						{
							//if(nrChains < 2)
							//	std::cout << timeToPosChainStart[timeForEval].at(nrChains) << std::endl;

							// Corrected by the COM-Diffusion of all particles
							VectorDouble3 diffRcom = (timeToPosCOMChain[timeForEval].at(nrObjects)-timeToPosCOMBox[timeForEval])-(timeToPosCOMChain[timeForEvalRef].at(nrObjects)-timeToPosCOMBox[timeForEvalRef]);

							VectorDouble3 diffRN2 = (timeToPosChainHalf[timeForEval].at(nrObjects)-timeToPosCOMBox[timeForEval])-(timeToPosChainHalf[timeForEvalRef].at(nrObjects)-timeToPosCOMBox[timeForEvalRef]);

							VectorDouble3 diffRFirst = (timeToPosChainStart[timeForEval].at(nrObjects)-timeToPosCOMBox[timeForEval])-(timeToPosChainStart[timeForEvalRef].at(nrObjects)-timeToPosCOMBox[timeForEvalRef]);

							VectorDouble3 diffRLast = (timeToPosChainEnd[timeForEval].at(nrObjects)-timeToPosCOMBox[timeForEval])-(timeToPosChainEnd[timeForEvalRef].at(nrObjects)-timeToPosCOMBox[timeForEvalRef]);

							VectorDouble3 diffRN2MRcom = (timeToPosChainHalf[timeForEval].at(nrObjects)-timeToPosCOMChain[timeForEval].at(nrObjects))-(timeToPosChainHalf[timeForEvalRef].at(nrObjects)-timeToPosCOMChain[timeForEvalRef].at(nrObjects));


							uint64_t differenceTime = timeForEvalRef - timeForEval;

							g1[differenceTime].AddValue(diffRN2*diffRN2);

							g2[differenceTime].AddValue(diffRN2MRcom*diffRN2MRcom);

							g3[differenceTime].AddValue(diffRcom*diffRcom);

							//g4[differenceTime].AddValue(diffRFirst*diffRFirst);
							g4[differenceTime].AddValue(diffRLast*diffRLast);

						}
					}
				}
			}
		}
	}


	std::cout  << std::endl << "Size g1: " << g1.size() << std::endl;

	// find the filename without path and extensions
	std::string filenameGeneral= std::string( std::find_if( ingredients.getName().rbegin(), ingredients.getName().rend(), PathSeparator() ).base(), ingredients.getName().end() );

	std::string::size_type const p(filenameGeneral.find_last_of('.'));
	filenameGeneral = filenameGeneral.substr(0, p);




	std::vector < std::vector<double> > tmpResultsMSD;

	// 14 columns
	tmpResultsMSD.resize(14);

	//for(int i = 0; i < 14; i++)
	//	tmpResultsMSD[i].resize(g1.size());


	for ( std::map<uint64_t, StatisticMoment>::iterator it=g1.begin(); it!=g1.end(); ++it)
	{
		tmpResultsMSD[0].push_back( it->first );
		tmpResultsMSD[1].push_back( it->second.ReturnM1() );
		tmpResultsMSD[2].push_back( it->second.ReturnM2() );
		tmpResultsMSD[3].push_back( 2.0* it->second.ReturnSigma()/sqrt(1.0*it->second.ReturnN() ));

	}

	for ( std::map<uint64_t, StatisticMoment>::iterator it=g2.begin(); it!=g2.end(); ++it)
		{
			tmpResultsMSD[4].push_back( it->second.ReturnM1() );
			tmpResultsMSD[5].push_back( it->second.ReturnM2() );
			tmpResultsMSD[6].push_back( 2.0* it->second.ReturnSigma()/sqrt(1.0*it->second.ReturnN() ));
		}

	for ( std::map<uint64_t, StatisticMoment>::iterator it=g3.begin(); it!=g3.end(); ++it)
		{
			tmpResultsMSD[7].push_back( it->second.ReturnM1() );
			tmpResultsMSD[8].push_back( it->second.ReturnM2() );
			tmpResultsMSD[9].push_back( 2.0* it->second.ReturnSigma()/sqrt(1.0*it->second.ReturnN() ));
		}

	for ( std::map<uint64_t, StatisticMoment>::iterator it=g4.begin(); it!=g4.end(); ++it)
		{
			tmpResultsMSD[10].push_back( it->second.ReturnM1() );
			tmpResultsMSD[11].push_back( it->second.ReturnM2() );
			tmpResultsMSD[12].push_back( 2.0* it->second.ReturnSigma()/sqrt(1.0*it->second.ReturnN() ));
			tmpResultsMSD[13].push_back( it->second.ReturnN() );

		}

	std::stringstream comment;
	comment <<"NrOfMonomersPerChain=" << numMonoPerObjects << "\n"
			<<"NrOfChains=" <<  numObjects <<  "\n"
			<<"NrDensity c="<<((8.0*(ingredients.getMolecules().size()))/(ingredients.getBoxX()*ingredients.getBoxY()*ingredients.getBoxZ())) << "\n"
			<<"\n"
			<<"dt  g1  (g1)^2 d(g1) g2  (g2)^2 d(g2) g3  (g3)^2 d(g3) g4  (g4)^2 d(g4) SampleSize" << "\n";


	std::stringstream ssD;	ssD << diffFrame;
	std::stringstream ssF;	ssF << startFrame;

	std::string filenameMSD = filenameGeneral +"_MSD_g1234_Time_df"+ssD.str()+"_f0_"+ssF.str()+".dat";

	std::cout  << " Write output to: " << dstdir  <<"/" << filenameMSD << std::endl;

	ResultFormattingTools::writeResultFile(dstdir+"/"+filenameMSD, this->ingredients, tmpResultsMSD, comment.str());

*/
/*
	// print results into a file
	// get the filename and path
	// find the filename without path and extensions
			std::string filenameGeneral= std::string( std::find_if( ingredients.getName().rbegin(), ingredients.getName().rend(), PathSeparator() ).base(), ingredients.getName().end() );

			std::string::size_type const p(filenameGeneral.find_last_of('.'));
			filenameGeneral = filenameGeneral.substr(0, p);


	// construct a list
	std::vector < std::vector<double> > tmpResultsHG_ZCOM;

	// we have 3 columns and row
	uint32_t columns = 2;
	uint32_t rows = HG_Statistik_TopoDistance_Walking.getNBins();



	// we have columns
	tmpResultsHG_ZCOM.resize(columns);

	// we have rows
	for(int i = 0; i < columns; i++)
		tmpResultsHG_ZCOM[i].resize(rows);

	// fill the list
				for(size_t n=0;n<HG_Statistik_TopoDistance_Walking.getNBins();n++)
				{
					tmpResultsHG_ZCOM[0][n] = HG_Statistik_TopoDistance_Walking.getLowerBoundaryOfBin(n);
					tmpResultsHG_ZCOM[1][n] = HG_Statistik_TopoDistance_Walking.getFirstMomentInBin(n);
					//tmpResultsHG_ZCOM[2][n] = HG_Statistik_TopoDistance_Walking.getFirstMomentInBin(n)/Statistic_ZCOM.ReturnM1();
					//tmpResultsHG_ZCOM[3][n] = HG_TopoDistance.GetNrInBin(n)/(1.0*nrOfSamples);

				}



	std::stringstream comment;
	comment << " File produced by analyzer Analyzer_Dijkstra_TopologicalDistance_ZCOM" << std::endl
			<< " Histogram of topological distance from center with " << ingredients.getMolecules().size() << " monomers" << std::endl
		//	<< " Mean zCOM <zCOM>=" << Statistic_ZCOM.ReturnM1() << std::endl
			<< std::endl
			<< "topological-distance <z> <z>/<zCOM>";



	//new filename
	std::string filenameHG_ZCOM = filenameGeneral + "_HG_TopologicalDistance_zMean.dat";

	ResultFormattingTools::writeResultFile(dstdir+"/"+filenameHG_ZCOM, this->ingredients, tmpResultsHG_ZCOM, comment.str());

	*/
}

#endif /*Analyzer_TopologicalPropertiesTree_H*/
