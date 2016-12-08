//
// Copyright (c) 2009-2015 Glen Berseth, Mubbasir Kapadia, Shawn Singh, Petros Faloutsos, Glenn Reinman
// See license.txt for complete license.
// Copyright (c) 2015 Mahyar Khayatkhoei
//

#include <algorithm>
#include <vector>
#include <util/Geometry.h>
#include <util/Curve.h>
#include <util/Color.h>
#include <util/DrawLib.h>
#include "Globals.h"

using namespace Util;

Curve::Curve(const CurvePoint& startPoint, int curveType) : type(curveType)
{
	controlPoints.push_back(startPoint);
}

Curve::Curve(const std::vector<CurvePoint>& inputPoints, int curveType) : type(curveType)
{
	controlPoints = inputPoints;
	sortControlPoints();
}

// Add one control point to the vector controlPoints
void Curve::addControlPoint(const CurvePoint& inputPoint)
{
	controlPoints.push_back(inputPoint);
	sortControlPoints();
}
// Add a vector of control points to the vector controlPoints
void Curve::addControlPoints(const std::vector<CurvePoint>& inputPoints)
{
	for (int i = 0; i < inputPoints.size(); i++)
		controlPoints.push_back(inputPoints[i]);
	sortControlPoints();
}

// Draw the curve shape on screen, usign window as step size (bigger window: less accurate shape)
void Curve::drawCurve(Color curveColor, float curveThickness, int window)
{

	#ifdef ENABLE_GUI
	
	if(checkRobust())
	{
		//Since the list is sorted last/greatest time will be the last element and the first/least time will be the first element
		float endTime=controlPoints.at(controlPoints.size()-1).time;
		float startTime=controlPoints.at(0).time;
		
	
		float timeAtCurve;
		Point currentPoint;
		Point nextPoint;
		
		//the amount of intervals the curve was split up into
		int timeIntervalSize=(int)((endTime-startTime)/window);
		
		// Move on the curve from t=0 to t=finalPoint, using window as step size, and linearly interpolate the curve points
		currentPoint=controlPoints.at(0).position;
		timeAtCurve=controlPoints.at(0).time+window;
		for(int a=0;a<timeIntervalSize-1;a+=1)
		{
			//uses timeAtCurve to compute the next point on the curve and stores it in the variable nextPoint
			if(calculatePoint(nextPoint,timeAtCurve)==false)
			{
				std::cerr << "ERROR. Calculate point failed" << std::endl;
				return;
			}
			
			//draws part of the curve (aka a straight line that connects currentPoint and nextPoint)
			DrawLib::drawLine(currentPoint,nextPoint,curveColor,curveThickness);
			
			//moves to the next part of the curve 
			timeAtCurve+=window;
			currentPoint=nextPoint;
			

		}
		//gets the last point of the curve
		nextPoint=controlPoints.at(controlPoints.size()-1).position;
		
		//draws the last part of the curve 
		DrawLib::drawLine(currentPoint,nextPoint,curveColor,curveThickness);

		
	}

	
	
	#endif
}

bool compareControlPoints(CurvePoint a, CurvePoint b)
{
	return (a.time < b.time);
}

// Sort controlPoints vector in ascending order: min-first
void Curve::sortControlPoints()
{
	std::sort(Curve::controlPoints.begin(), Curve::controlPoints.end(), compareControlPoints);
	for (int i = 0; i < Curve::controlPoints.size() - 1; i++) {
		if (Curve::controlPoints.at(i).time == Curve::controlPoints.at(i + 1).time) {
			Curve::controlPoints.erase(Curve::controlPoints.begin() + i + 1);
		}
	}
	return;
}

// Calculate the position on curve corresponding to the given time, outputPoint is the resulting position
// Note that this function should return false if the end of the curve is reached, or no next point can be found
bool Curve::calculatePoint(Point& outputPoint, float time)
{
	// Robustness: make sure there is at least two control point: start and end points
	if (!checkRobust())
		return false;

	// Define temporary parameters for calculation
	unsigned int nextPoint;

	// Find the current interval in time, supposing that controlPoints is sorted (sorting is done whenever control points are added)
	// Note that nextPoint is an integer containing the index of the next control point
	if (!findTimeInterval(nextPoint, time))
		return false;

	// Calculate position at t = time on curve given the next control point (nextPoint)
	if (type == hermiteCurve)
	{
		outputPoint = useHermiteCurve(nextPoint, time);
	}
	else if (type == catmullCurve)
	{
		outputPoint = useCatmullCurve(nextPoint, time);
	}

	// Return
	return true;
}

// Check Robustness
bool Curve::checkRobust()
{
	
		if (controlPoints.size()<2)
		{
			return false;
		}
		return true;
	
}

// Find the current time interval (i.e. index of the next control point to follow according to current time)
bool Curve::findTimeInterval(unsigned int& nextPoint, float time)
{
	//Finds the control point that has a greater time value than the given time. 
	int a;
	for (a = 0; a < controlPoints.size() && controlPoints.at(a).time <= time ; a++)
	{
		
	}
	
	if (a==controlPoints.size())
	{
		return false; 
	}	
	nextPoint = a;
	return true;
	
}

// Implement Hermite curve
Point Curve::useHermiteCurve(const unsigned int nextPoint, const float time)
{
	Point newPosition;
	float normalTime, intervalTime;

	intervalTime = controlPoints[nextPoint].time - controlPoints[nextPoint - 1].time;

	normalTime = (time - controlPoints[nextPoint - 1].time) / intervalTime;



	float f1 = 2 * pow(normalTime, 3) - 3 * pow(normalTime, 2) + 1;
	float f2 = -2 * pow(normalTime, 3) + 3 * pow(normalTime, 2);
	float f3 = pow(normalTime, 3) - 2 * pow(normalTime, 2) + normalTime;
	float f4 = pow(normalTime, 3) - pow(normalTime, 2);


	Point p1 = controlPoints[nextPoint - 1].position;
	Vector t1 = controlPoints[nextPoint - 1].tangent;
	Point p2 = controlPoints[nextPoint].position;
	Vector t2 = controlPoints[nextPoint].tangent;

	newPosition = f1*p1 + f2*p2 + f3*t1*intervalTime + f4*t2*intervalTime;


	return newPosition;
}

// Implement Catmull-Rom curve
Point Curve::useCatmullCurve(const unsigned int nextPoint, const float time)
{	
	Point newPosition;
	float normalTime, intervalTime;

	intervalTime = controlPoints[nextPoint].time - controlPoints[nextPoint - 1].time;

	normalTime = (time - controlPoints[nextPoint - 1].time) / intervalTime;



	Point p1 = controlPoints[nextPoint - 1].position;
	Point p2 = controlPoints[nextPoint].position;
	float ti1 = controlPoints[nextPoint - 1].time;
	float ti2 = controlPoints[nextPoint].time;

	Vector t1, t2;
	if (controlPoints.size() == 2)
	{
		t1 = (controlPoints[1].position - controlPoints[0].position) / (controlPoints[1].time - controlPoints[0].time);
	}
	else if (controlPoints.size() == 2)
	{
		t1 = (controlPoints[1].position - controlPoints[0].position) / (controlPoints[1].time - controlPoints[0].time);
	}
	else if (nextPoint == 1)
	{
		Point p3 = controlPoints[nextPoint + 1].position;
		float ti3 = controlPoints[nextPoint + 1].time;
		t1 = (p2 - p1) / (ti2 - ti1) * (ti3 - ti1) / (ti3 - ti2) - (p3 - p1) / (ti3 - ti1) * (ti2 - ti1) / (ti3 - ti2);
		t2 = (p2 - p1) / (ti2 - ti1) * (ti3 - ti2) / (ti3 - ti1) + (p3 - p2) / (ti3 - ti2) * (ti2 - ti1) / (ti3 - ti1);
	}
	else if (nextPoint == (controlPoints.size() - 1))
	{
		Point p0 = controlPoints[nextPoint - 2].position;
		float ti0 = controlPoints[nextPoint - 2].time;
		t1 = (p1 - p0) / (ti1 - ti0) * (ti2 - ti1) / (ti2 - ti0) + (p2 - p1) / (ti2 - ti1) * (ti1 - ti0) / (ti2 - ti0);
		t2 = (p2 - p0) / (ti2 - ti0) * (ti1 - ti0) / (ti2 - ti1) - (p1 - p2) / (ti1 - ti0) * (ti2 - ti0) / (ti2 - ti1);

	}
	else
	{
		Point p0 = controlPoints[nextPoint - 2].position;
		float ti0 = controlPoints[nextPoint - 2].time;
		Point p3 = controlPoints[nextPoint + 1].position;
		float ti3 = controlPoints[nextPoint + 1].time;
		t1 = (p1 - p0) / (ti1 - ti0) * (ti2 - ti1) / (ti2 - ti0) + (p2 - p1) / (ti2 - ti1) * (ti1 - ti0) / (ti2 - ti0);
		t2 = (p2 - p1) / (ti2 - ti1) * (ti3 - ti2) / (ti3 - ti1) + (p3 - p2) / (ti3 - ti2) * (ti2 - ti1) / (ti3 - ti1);


	}

	float f1 = 2 * pow(normalTime, 3) - 3 * pow(normalTime, 2) + 1;
	float f2 = -2 * pow(normalTime, 3) + 3 * pow(normalTime, 2);
	float f3 = pow(normalTime, 3) - 2 * pow(normalTime, 2) + normalTime;
	float f4 = pow(normalTime, 3) - pow(normalTime, 2);


	newPosition = f1*p1 + f2*p2 + f3*t1*intervalTime + f4*t2*intervalTime;

	return newPosition;
}