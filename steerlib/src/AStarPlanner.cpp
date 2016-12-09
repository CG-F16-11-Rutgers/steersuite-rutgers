//
// Copyright (c) 2009-2015 Glen Berseth, Mubbasir Kapadia, Shawn Singh, Petros Faloutsos, Glenn Reinman, Rahul Shome
// See license.txt for complete license.
//

#include <vector>
#include <stack>
#include <set>
#include <map>
#include <iostream>
#include <algorithm> 
#include <functional>
#include <queue>
#include <math.h>
#include "planning/AStarPlanner.h"


#define COLLISION_COST  1000
#define GRID_STEP  1
#define OBSTACLE_CLEARANCE 1
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))


namespace SteerLib
{
	AStarPlanner::AStarPlanner() {}

	AStarPlanner::~AStarPlanner() {}


	bool AStarPlanner::canBeTraversed(int id)
	{
		double traversal_cost = 0;
		int current_id = id;
		unsigned int x, z;
		gSpatialDatabase->getGridCoordinatesFromIndex(current_id, x, z);
		int x_range_min, x_range_max, z_range_min, z_range_max;

		x_range_min = MAX(x - OBSTACLE_CLEARANCE, 0);
		x_range_max = MIN(x + OBSTACLE_CLEARANCE, gSpatialDatabase->getNumCellsX());

		z_range_min = MAX(z - OBSTACLE_CLEARANCE, 0);
		z_range_max = MIN(z + OBSTACLE_CLEARANCE, gSpatialDatabase->getNumCellsZ());


		for (int i = x_range_min; i <= x_range_max; i += GRID_STEP)
		{
			for (int j = z_range_min; j <= z_range_max; j += GRID_STEP)
			{
				int index = gSpatialDatabase->getCellIndexFromGridCoords(i, j);
				traversal_cost += gSpatialDatabase->getTraversalCost(index);

			}
		}

		if (traversal_cost > COLLISION_COST)
			return false;
		return true;
	}

	//returns true if grid coord given by point can be traversed
	bool AStarPlanner::canBeTraversed1(Util::Point point) {
		return canBeTraversed(gSpatialDatabase->getCellIndexFromLocation(point));
	}

	//returns true if grid coord given by point can be traversed
	bool AStarPlanner::canBeTraversed(Util::Point point) {
		return canBeTraversed(canBeTraversed(gSpatialDatabase->getCellIndexFromLocation(point)));
	}


	Util::Point AStarPlanner::getPointFromGridIndex(int id)
	{
		Util::Point p;
		gSpatialDatabase->getLocationFromIndex(id, p);
		return p;
	}

	int getLowestFIndex(std::vector<AStarPlannerNode*> set) {
		int returnIndex;
		double bestFValue = DBL_MAX;

		for (int i = 0; i < set.size(); i++) {
			if (set[i]->f <= bestFValue) {
				bestFValue = set[i]->f;
				returnIndex = i;
			}
		}
		return returnIndex;
	}

	// checks if the Node is walkable,
	void AStarPlanner::addNeighborNode(AStarPlannerNode* parentNode, Util::Point location, std::vector<AStarPlannerNode*> &neighborNodes) {
		AStarPlannerNode* currentNode = new AStarPlannerNode(location, double(0), double(0), parentNode);
		neighborNodes.push_back(currentNode);
	}

	bool inSet(AStarPlannerNode* currentNode, std::vector<AStarPlannerNode*> set) {
		for (int i = 0; i < set.size(); i++) {
			if (currentNode->point == set[i]->point) {
				return true;
			}
		}
		return false;
	}

	//Stuff that I am going to try out here
	bool AStarPlanner::checkAddSuccessor(AStarPlannerNode* parentNode, Util::Point location, std::vector<AStarPlannerNode*> &successors, Util::Point goal) {
		if (!canBeTraversed1(location)) {
			//std::cout << "The location: " << location << " is blocked\n";
			return false;
		}

		//std::cout << "The location: " << location << " can be traversed and is added to the list of successors\n";
		AStarPlannerNode* successorNode = new AStarPlannerNode(location, parentNode->g + distanceBetween(parentNode->point, location), double(0), parentNode);
		successors.push_back(successorNode);
		return true;
	}

	std::vector<AStarPlannerNode*> AStarPlanner::getSuccessors(AStarPlannerNode* currentNode, Util::Point goal) {
		std::vector<AStarPlannerNode*> successors;

		int size = 3;
		int foo[3] = { -1,0,1 };
		for (int a = 0; a < size; a++)
		{
			for (int b = 0; b < size; b++)
			{
				if (foo[a] == 0 && foo[b] == 0)
				{
					continue;
				}
				checkAddSuccessor(currentNode,
					Util::Point(currentNode->point.x + foo[a], currentNode->point.y, currentNode->point.z + foo[b]),
					successors, goal);
			}
		}

		return successors;
	}


	bool AStarPlanner::weightedAStar(std::vector<Util::Point>& agent_path, Util::Point start,
		Util::Point goal, bool append_to_path)
	{
		double inflationFactor = 1;
		std::cout << "\nWEIGHTED A*\n";

		std::vector<AStarPlannerNode*> closedSet, openSet, path;
		openSet.push_back(new AStarPlannerNode(start, double(0), distanceBetween(start, goal), NULL));

		while (openSet.size() != 0) {
			//current node is node with lowest f-score value (computed and assigned here)
			int index = getLowestFIndex(openSet);
			AStarPlannerNode* current = openSet[index];

			//The current node is the goal. Assign closed set to agent_path
			//std::cout << "\n//////////The current position is " << current->point << "//////////\n";
			if (current->point == goal) {
				std::cout << "Target Found!!\n";
				//path.push_back(current);
				AStarPlannerNode* traceBack = current;
				while (traceBack != NULL) {
					path.push_back(traceBack);
					traceBack = traceBack->parent;
				}


				storePath(path, agent_path, goal);
				return true;
			}

			//remove current node from the openset
			openSet.erase(openSet.begin() + index);
			//add current node to closedset
			closedSet.push_back(current);


			std::vector<AStarPlannerNode*> neighborNodes = getSuccessors(current, goal);

			for (int i = 0; i < neighborNodes.size(); i++) {

				//if neighbor in closed set then continue
				if (inSet(neighborNodes[i], closedSet)) {
					continue;
				}

				//get the tentative g score
				double tentative_g = current->g + distanceBetween(current->point, neighborNodes[i]->point);

				//std::cout << "\nTentative g = " << tentative_g << " and neighborNodes[i]->g = " << neighborNodes[i]->g << "\n";
				if (tentative_g <= neighborNodes[i]->g) {
					neighborNodes[i]->g = tentative_g;
					neighborNodes[i]->f = neighborNodes[i]->g + inflationFactor*distanceBetween(neighborNodes[i]->point, goal);
					if (!inSet(neighborNodes[i], openSet)) {
						openSet.push_back(neighborNodes[i]);
					}
				}

			}

		}

		return false;
	}


	//////////////////////////////////////////////////COMPUTE PATH IS HERE////////////////////////////////////////////////////
	//Comment or uncomment the returns depending on what algorithms needs to run
	bool AStarPlanner::computePath(std::vector<Util::Point>& agent_path, Util::Point start, Util::Point goal, SteerLib::SpatialDataBaseInterface * _gSpatialDatabase, bool append_to_path)
	{
		std::cout << "\nSTART COMPUTE PATH\n";

		gSpatialDatabase = _gSpatialDatabase;

		return weightedAStar(agent_path, start, goal, append_to_path);
		//return ARAStar(agent_path, start, goal, append_to_path);
		//return ADStar(agent_path, start, goal, append_to_path);
	}
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	bool AStarPlanner::ARAStar(std::vector<Util::Point>& agent_path, Util::Point start, Util::Point goal, bool append_to_path)
	{
		double inflationFactor = 1000;
		double timeLimitInSeconds = 1;
		double decrementFactor = 1;

		time_t now;
		time_t later;


		time(&now);  // get current time; same as: now = time(NULL)  


		std::vector<SteerLib::AStarPlannerNode*> openList, closedList, inconsistentList, path;
		//std::cout << "\ncreated lists";


		SteerLib::AStarPlannerNode* goalState = new SteerLib::AStarPlannerNode(goal, DBL_MAX, DBL_MAX, NULL);
		//std::cout << "\ncreated goalState";

		//inserting start state into openList
		openList.push_back(new SteerLib::AStarPlannerNode(start, 0, 0, NULL));
		//std::cout << "\nadded start state";

		ImprovePath(openList, closedList, inconsistentList, path, goal, goalState, inflationFactor);
		if (path.empty())
		{
			std::cout << "\nNo Path exists\n";
			return false;
		}

		double condition = min(inflationFactor, goalState->g) / getMin(openList, inconsistentList, goal);
		//std::cout << "\ncondition is "<<condition<<"\n";
		time(&later);
		//std::cout << "\ndifftime(later, now) is " << difftime(later, now) << "\n";
		while (timeLimitInSeconds > difftime(later, now))
		{
			//std::cout << "\nloop \n";
			//std::cout << "\ndifftime(later, now) is " << difftime(later, now) << "\n";
			inflationFactor -= decrementFactor;
			if (inflationFactor < 1)
				inflationFactor = 1;
			//move states from incon into OPen
			while (!inconsistentList.empty())
			{
				SteerLib::AStarPlannerNode* node = inconsistentList.at(0);

				inconsistentList.pop_back();
				openList.push_back(node);
			}

			//update f value for all open
			for (int a = 0; a < openList.size(); a++)
			{
				openList.at(a)->f = openList.at(a)->g + inflationFactor*distanceBetween(openList.at(a)->point, goal);
			}
			//selectionSort(openList);

			//std::cout << "after sort";
			//printList(openList);
			path.clear();
			closedList.clear();
			fixList(openList, start);

			ImprovePath(openList, closedList, inconsistentList, path, goal, goalState, inflationFactor);
			if (inflationFactor == 1)
			{
				//std::cout << "\optimal solution was found\n";
				break;
			}
			condition = min(inflationFactor, goalState->g) / getMin(openList, inconsistentList, goal);

			time(&later);


			//if(condition <= 1 )
			//std::cout << "\ncondition is good\n";
			//else if(timeLimitInSeconds <= difftime(later, now))
			//std::cout << "\ntime is up\n";


		}
		storePath(path, agent_path, goal);

		return true;

	}
	double AStarPlanner::getMin(std::vector<SteerLib::AStarPlannerNode*> openList,
		std::vector<SteerLib::AStarPlannerNode*> inconsistentList, Util::Point goal)
	{
		double min = DBL_MAX;

		for (int a = 0; a < openList.size(); a++)
		{
			if ((openList.at(a)->g + Util::distanceBetween(openList.at(a)->point, goal)) < min)
				min = openList.at(a)->g + Util::distanceBetween(openList.at(a)->point, goal);
		}
		for (int a = 0; a < inconsistentList.size(); a++)
		{
			if ((inconsistentList.at(a)->g + Util::distanceBetween(inconsistentList.at(a)->point, goal)) < min)
				min = inconsistentList.at(a)->g + Util::distanceBetween(inconsistentList.at(a)->point, goal);
		}
		return min;

	}

	void AStarPlanner::storePath(std::vector<SteerLib::AStarPlannerNode*> &path, std::vector<Util::Point>& agent_path, Util::Point goal)
	{
		for (int i = path.size() - 1; i > 0; i--) {
			//std::cout << "list size = " << closedSet.size() << ". loop " << i << "\n";
			//std::cout << "pushing point " << path[i]->point << "\n";
			agent_path.push_back(path[i]->point);
		}
		agent_path.push_back(goal);
	}

	void AStarPlanner::ImprovePath(std::vector<SteerLib::AStarPlannerNode*> &openSet,
		std::vector<SteerLib::AStarPlannerNode*> &closedSet,
		std::vector<SteerLib::AStarPlannerNode*> &inconsistentList,
		std::vector<SteerLib::AStarPlannerNode*> &path, Util::Point goal,
		SteerLib::AStarPlannerNode* &goalState, double inflationFactor)
	{

		std::cout << "\nImprovePath*\n";

		while (openSet.size() != 0) {
			//current node is node with lowest f-score value (computed and assigned here)
			int index = getLowestFIndex(openSet);

			AStarPlannerNode* current = openSet[index];

			//The current node is the goal. Assign closed set to agent_path
			//std::cout << "\n//////////The current position is " << current->point << "//////////\n";
			if (current->point == goal) {
				std::cout << "Target Found!!\n";
				//path.push_back(current);
				AStarPlannerNode* traceBack = current;
				while (traceBack != NULL) {
					path.push_back(traceBack);
					traceBack = traceBack->parent;
				}

				goalState->g = current->g;


				return;
			}

			//remove current node from the openset
			openSet.erase(openSet.begin() + index);
			//add current node to closedset
			closedSet.push_back(current);


			std::vector<AStarPlannerNode*> neighborNodes = getSuccessors(current, goal);

			for (int i = 0; i < neighborNodes.size(); i++) {

				//if neighbor in closed set then continue
				if (inSet(neighborNodes[i], closedSet)) {
					continue;
				}

				//get the tentative g score
				double tentative_g = current->g + distanceBetween(current->point, neighborNodes[i]->point);

				//std::cout << "\nTentative g = " << tentative_g << " and neighborNodes[i]->g = " << neighborNodes[i]->g << "\n";
				if (tentative_g <= neighborNodes[i]->g) {
					neighborNodes[i]->g = tentative_g;
					neighborNodes[i]->f = neighborNodes[i]->g + inflationFactor*distanceBetween(neighborNodes[i]->point, goal);
					if (!inSet(neighborNodes[i], openSet)) {
						openSet.push_back(neighborNodes[i]);
					}
				}
			}
		}
	}
	void AStarPlanner::fixList(std::vector<SteerLib::AStarPlannerNode*> &o, Util::Point p)
	{
		o.clear();
		o.push_back(new SteerLib::AStarPlannerNode(p, 0, 0, NULL));
	}

	//AD Star stuff

	//Needs to return a 2-tuple (Done)
	std::tuple<double, double> AStarPlanner::key(AStarPlannerNode* startNode, AStarPlannerNode* node) {
		if (node->g > node->rhs) {
			return{ (node->rhs + 1/* EPS*/ * distanceBetween(startNode->point, node->point)), node->rhs };
		}
		else {
			return{ (node->g + distanceBetween(startNode->point, node->point)), node->g };
		}
	}

	//Done
	bool AStarPlanner::UpdateState(AStarPlannerNode* startNode, AStarPlannerNode* s, AStarPlannerNode* goal, std::vector<AStarPlannerNode*> openList, std::vector<AStarPlannerNode*> closedList, std::vector<AStarPlannerNode*> inconsList)
	{
		if (std::find(closedList.begin(), closedList.end(), s) == closedList.end())
		{
			s->g = DBL_MAX;
		}
		if (s->point != goal->point)
		{
			std::vector<AStarPlannerNode*> successors = getSuccessors(s, goal->point);

			float minS = DBL_MAX;
			for (int i = 0; i < successors.size(); i++) {
				float f = successors[i]->g + distanceBetween(s->point, successors[i]->point);
				if (f < minS) {
					minS = f;
				}
			}
			s->rhs = minS;
		}
		if (std::find(openList.begin(), openList.end(), s) != openList.end())
		{
			openList.erase(std::find(openList.begin(), openList.end(), s)); // remove curr from openNodes
		}
		if (s->g != s->rhs)
		{
			if (std::find(closedList.begin(), closedList.end(), s) == closedList.end())
			{
				s->key = key(startNode, s);
				openList.push_back(s);
			}
			else
			{
				inconsList.push_back(s);
			}
		}
	}

	//Done
	bool keyLessthan(AStarPlannerNode* s1, AStarPlannerNode* s2) {

		if (std::get<0>(s1->key) < std::get<0>(s2->key) || (std::get<0>(s1->key) == std::get<0>(s2->key) && (std::get<1>(s1->key) < std::get<1>(s2->key)))) {
			return true;
		}
		return false;
	}


	double epsilon;
	std::vector<AStarPlannerNode*> openList, closedList, inconsistentList;
	AStarPlannerNode* startNode;
	bool AStarPlanner::ADStar(std::vector<Util::Point>& agent_path, Util::Point start, Util::Point goal, bool append_to_path) {

		std::cout << "Starting AD*\n";

		//set values for starting node
		startNode = new AStarPlannerNode(start, DBL_MAX, double(0), NULL);
		startNode->rhs = DBL_MAX;
		epsilon = 4.0;

		//Push goal node into openlist (default value of rhs = 0)
		AStarPlannerNode* goalNode = new AStarPlannerNode(goal, DBL_MAX, double(0), NULL);
		goalNode->key = key(startNode, goalNode);
		openList.push_back(goalNode);
		computeOrImprovePath(openList, closedList, inconsistentList, startNode, goalNode);
		//PUBLISH TO AGENT_PATH

		while (true /*TODO time remains*/) {
			ADStarUpdate(agent_path, start, goal, append_to_path, goalNode);
		}

		return true;
	}

	//called if changes in edge costs detected
	bool AStarPlanner::ADStarUpdate(std::vector<Util::Point>& agent_path, Util::Point start, Util::Point goal, bool append_to_path, AStarPlannerNode* goalNode) {
		if (append_to_path) {
			for (AStarPlannerNode * node : openList, inconsistentList, closedList) {
				//UpdateState(node);
			}
			epsilon *= 2;
		}
		else {
			epsilon = max(1.0, epsilon / 2); //halve epsilon, must have minimum value of 1.0
		}
		for (int i = 0; i < inconsistentList.size(); i++) {
			AStarPlannerNode* node = inconsistentList[i];
			inconsistentList.erase(inconsistentList.begin() + i);
			openList.push_back(node);
		}
		for (AStarPlannerNode* node : openList) {
			node->key = key(startNode, node);
		}
		for (int i = 0; i < closedList.size(); i++) {
			closedList.erase(closedList.begin() + i);
		}
		computeOrImprovePath(openList, closedList, inconsistentList, startNode, goalNode);
		//PUBLISH TO AGENT_PATH
		return true;
	}

	//Done
	void AStarPlanner::computeOrImprovePath(std::vector<AStarPlannerNode*> openList, std::vector<AStarPlannerNode*> closedList, std::vector<AStarPlannerNode*> inconsList, AStarPlannerNode* startNode, AStarPlannerNode* goal) {
		int index = getLowestFIndex(openList);
		AStarPlannerNode* currentNode = openList[index];
		while (keyLessthan(currentNode, startNode) || startNode->rhs != startNode->g) {
			openList.erase(openList.begin() + index);
			if (currentNode->g > currentNode->rhs) {
				currentNode->g = currentNode->rhs;
				closedList.push_back(currentNode);

				AStarPlannerNode* traceBack = currentNode->parent;
				while (traceBack != NULL) {
					UpdateState(startNode, traceBack, goal, openList, closedList, inconsList);
					traceBack = traceBack->parent;
				}

			}
			else
			{
				currentNode->g = DBL_MAX;
				AStarPlannerNode* traceBack = currentNode;
				while (traceBack != NULL) {
					UpdateState(startNode, traceBack, goal, openList, closedList, inconsList);
					traceBack = traceBack->parent;
				}

				//for all s' contained in predecessors including s updateState(s')
			}
		}
	}

}