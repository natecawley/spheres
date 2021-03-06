#include <stdlib.h>
#include <time.h>
#include <random>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <stdexcept>

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923078164062

/* Project Overview
	This project runs simulations of spheres in a box.
	
	The spheres have some potential function associated with them, which changes with each type. 
	The code is designed so that new sphere types can be created that overload methods in the sphere 
	class, which allows very general code to be written and tested, without having to rewrite very much 
	code every time that you develop a new model sphere. The way this is implemented is using pointers, 
	which are polymorphic by nature.
	
	Author: Nate Cawley
*/


/* TODO:
	Method for determining equilibrium condition for metropolis monte carlo: energy stabilizes within some tolerance
	Method for saving sphere data for individual spheres
	Method for reading data from a string and applying it to sphere pointers
	Method for calculating energy of whole spherebox, maybe by storing energy in each individual sphere?
	Method for creating a sphere box with uniformly distributed spheres
*/

/* 	Compiler Note
	Program uses libraries from the C++ 11 Standard Libraries
	To access these, compile using -std=c++0x (c++14 should work, too)
*/
using namespace std;

// Global variables of doom
mt19937_64 gen;
double gMax = gen.max();

/* Classes defined
Basic:
    Sphere: default sphere, is a hard sphere
    SphereBox: a class that holds spheres, a bit more structure added than just an array
Currently in development:
    WellSphere: A sphere with square well
Future Plans
    PatchyOneSphere: A patchy (anisotropic square well) sphere that interacts with hard spheres (or really, treats any other sphere as hard) (needs a new perturb function)
    PatchyTwoSphere: A patchy sphere that only interacts with others of its type
*/

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/* Class definition for Sphere
	Spheres are defined with a location and radius
	Non-Default Constructor defines radius, location.
	They interact with other spheres via getDistance() and getPotential()
	A lot of the other methods store data that is then overloaded 
*/
class Sphere{
protected:
    double radius;
    vector<double> location;
	const static int sphereType = 1;
    /*
    Methods defined for sphere:
        getRadius: returns radius
        setRadius: changes radius to a set value
        setLoc: sets the position to a certain value
        getLoc: returns the location of a particle
        getDistance: gets the distance to another instance of class sphere
        getPotential: gets the distance to another instance of class sphere (uses getDistance)
        getInteractionLength: overloaded method that determines maximum length that two spheres of this type would interact
        perturb: move location by a small random bit, magnitude determined by user
        readData: read data on location, size, etc from a string
        getData: create a string with data on location, size, etc from a string
    */
public:
    Sphere ();
    Sphere (double);
    Sphere (double, vector<double>);
	Sphere (string);
    double getRadius();
    void setRadius(double);
    void setLoc(vector<double>);
    vector<double> getLoc();
    void readData(string);
    string getData();
	string getHeader();
	const static int getType();
    double getDistance(Sphere*,vector<double> bounds);
    double getPotential(Sphere*,vector<double> bounds);
    double getInteractionLength();
    void perturb(double);
};

//Default Constructor
Sphere::Sphere(){
    // No default values
}

// Constructor for set radius
Sphere::Sphere(double rad){
    radius = rad;
}

// Constructor for preset radius and position (not to be overloaded)
Sphere::Sphere(double rad, vector<double> newLoc){
    radius = rad;
    location = newLoc;
}

// Constructor from string (use read data)
Sphere::Sphere(string line){
	readData(line);
}

// Set the radius of a sphere (not to be overloaded)
void Sphere::setRadius(double rad){
    radius = rad;
}

// Return the radius of a sphere (not to be overloaded)
double Sphere::getRadius(){
    return radius;
}

// Method for storing location of a sphere (not to be overloaded)
void Sphere::setLoc(vector<double> loc){
    location = loc;
}

// Return the location of a sphere (not to be overloaded)
vector<double> Sphere::getLoc(){
    return location;
}

// Find the minimum distance between two spheres, accounting for periodic boundary.
double Sphere::getDistance(Sphere *s, vector<double> bounds){
    vector<double> otherLoc = s->getLoc();
    double normsq = 0;
    double locDiff;
    // Maybe add a throw statement in case thisLoc and otherLoc are different lengths?
    for(int i = 0; i < 3; i++)
    {
        // Calculate separation difference
        locDiff = (location[i]-otherLoc[i]);
        // If greater than bounds/2, actually closer because of periodic boundary conditions.
        if(locDiff > bounds[i]/2)
        {
            locDiff = bounds[i]-locDiff;
        }
        // Take potentially adjusted difference squared
        normsq += pow(locDiff,2);
    }
    double norm = sqrt(normsq);
    return norm;
}

// Default spheres are hard spheres, meaning they can't overlap but otherwise don't interact.
// This will get overwritten in sticky spheres, patchy spheres, etc.
double Sphere::getPotential(Sphere* s, vector<double> bounds){
    // Calculate the minimum allowable distance between the spheres
    int minDist = radius + s->getRadius();
    // If the distance is less than allowed, set potential to infinity
    if (getDistance(s,bounds) <= minDist){
        return INFINITY;
    }
    // Otherwise default is hard sphere
    else {
            return 0;
        }
}

/* Method to return interaction length
	Interaction length is determined by sphere
*/
double Sphere::getInteractionLength(){
    //Overload this if potential acts at a longer distance
    return 2*radius;
}

// This is supposed to do something to change the location of the sphere
void Sphere::perturb(double mag){
    double u,v,phi,theta,ranJump;
    vector<double> dir;
    // pick a random amount from 0 to mag
    ranJump = gen()*mag/(gMax);

    // pick a random direction (uniform distribution)
        // two random numbers u and v between 0 and 1;
        u = gen()/(gMax);
        v = gen()/(gMax);
        // some stuff i pulled from wolfram mathworld
        phi = 8*atan(1)*u;
        theta = acos(2*v-1);
        // trig stuff to get normal vector
        dir[1] = ranJump*sin(theta)*cos(phi);
        dir[2] = ranJump*sin(theta)*sin(phi);
        dir[3] = ranJump*cos(theta);


    // move a random amount from 0 to mag in length
	location[1]+=dir[1];
	location[2]+=dir[2];
	location[3]+=dir[3];
}

// Read in data and parse it into variables, then set the internal variables
void Sphere::readData(string line){
    double radius,xLoc,yLoc,zLoc;
    vector<double> locs;
    stringstream ss;
    ss << line;
    ss >> radius >> xLoc >> yLoc >> zLoc;
    setRadius(radius);
    locs.push_back(xLoc);
    locs.push_back(yLoc);
    locs.push_back(zLoc);
    setLoc(locs);
}

// Return data as a list of numbers
string Sphere::getData(){
    stringstream writeStr;
    string dataStream;
    vector<double> locs = getLoc();
    writeStr << "	" << getRadius() << "	" << locs[0] << "	" << locs[1] << "	" << locs[2];
    dataStream = writeStr.str();
    return dataStream;
}

// Return the type of the sphere
const int Sphere::getType()
{
	return sphereType;
}

// Return the header for this type of sphere
string Sphere::getHeader(){
    stringstream writeStr;
    string dataStream;
    writeStr << "SphereType: Hard Sphere";
    writeStr >> dataStream;
    return dataStream;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/* SphereBox is a box that holds spheres...
 Box is of defined size, using rectangular coordinates for now
 In rectangular box, one corner is assumed to be at the origin
*/
class SphereBox{

/*
    sphereNum: the number of spheres in the box (set by constructor)
    b1,b2,b3: the x,y,and z bounds of the box
    box: the vector that holds pointers to all of the spheres in the box
    sphereIter: an iterator to run through the elements of the vector
*/
    int sphereNum;
    vector<double> bounds;
    vector<Sphere*> box;
    int numPlaced;
/*  Methods:
    fillBox: fills the box with spheres, insures no overlap
    getNearSpheres: returns all sphere within a certain distance of a given sphere
    addTestSphere: Calculate energy from adding a test sphere
    printSpheres: Prints out a list of all spheres in box (WIP)
    readSpheres: Reads a text file with sphere data, and adds all of them to the box.
    getRandomCoords: returns a random location within the box
    getRandomSpheres: selects a random sphere
    setSphereData: Read info from a file to fill box
	askEquil: determines if the box is at equilibrium
	calcTotalEnergy: calculate total energy of the box
	getTotalEnergy: return total energy of the box
*/
public:
    // Memory retrieve
    vector<double> getBounds();
    double getMeanFreePath();
    Sphere* getRandomSphere();

    // Physics methods
    void fillBox();
    vector<Sphere *> getNearSpheres(Sphere*);
    double getEnergySpheres(Sphere*, vector<Sphere*>);
    double addTestSphere(Sphere * s);
	bool askEquil();
	void calcTotalEnergy();


    // File I/O
    void printSpheres();
    void printSpheres(string);
    void readSpheres(string);
    void setSphereData(string);

    // Empty constructor
    SphereBox();
    // Rectangular constructor
    SphereBox (vector<Sphere*>, double, double, double);

private:
    vector<double> getRandomCoords();
};

// Empty constructor for SphereBox
SphereBox::SphereBox(){
    // Empty Constructor for file I/O

}

/* Rectangular Constructor for SphereBox class
    Inputs:
        s: a vector holding pointers to spheres that do not yet have set locations, but do have set radii (all assumed to be identical)
        xMax, yMax, zMax: input rectangular bounds for the edges of the box
    Function:
        Determines if the input parameters are valid, sets private variables, determines whether or not the box can plausibly be filled with spheres
        Executes fillBox method (which defines location of spheres) if parameters are valid
*/
SphereBox::SphereBox (vector<Sphere*> s, double xMax, double yMax, double zMax){
    // Set internal values
	cout << "\nsetting internals\n";
    box = s;
	cout << "counting spheres\n";
    sphereNum = box.size();
	cout << "size = " << sphereNum;
	bounds = vector<double>(3,0);
    bounds[0] = xMax;
    bounds[1] = yMax;
    bounds[2] = zMax;
	cout << "\nconstructor ended";
}

// Uses random number generator to return a set of coordinates
vector<double> SphereBox::getRandomCoords(){
    vector<double> coords;
    coords= {((double) gen()/ gMax)*bounds[0],((double) gen()/ gMax)*bounds[1],((double) gen()/ gMax)*bounds[2]};
    return coords;
};

// return the bounds of the box
vector<double> SphereBox::getBounds(){
    return bounds;
}

// Finds all spheres within the interaction length of a given sphere
vector<Sphere*> SphereBox::getNearSpheres(Sphere* s){
    // Empty container to hold close spheres
    vector<Sphere*> nearSpheres;
    // New iterator to look through box for close spheres
    vector<double> kLoc;
    double cToCDistance;
    Sphere * holder;
    // Move sphereSlider until it hits an undefined sphere ( will be box.end() if fillBox has already completed)
    for(int k = 0; k < numPlaced; k++)
    {
        holder = box[k];
        kLoc = (holder->getLoc());
        try
        {
            cToCDistance = (*s).getDistance(holder,getBounds());
            if(cToCDistance < (*s).getInteractionLength())
            {
                nearSpheres.push_back(holder);
            }
        }
        catch(...)
        {
            cout << "The location of the checked sphere was not defined!";
        }
    }
    return nearSpheres;
}

// Calculates the energy of one sphere due to a set of other spheres.
double SphereBox::getEnergySpheres(Sphere* s, vector<Sphere*> nears){
    double ener;
    int nsize;
    Sphere* nSphere;
    ener = 0;
    nsize = nears.size();
    for(int n = 0; n < nsize; n++)
    {
        nSphere = nears[n];
        ener += s->getPotential(nSphere,getBounds());
    }
    return ener;
}

// Calculate the mean free path for perturbations
double SphereBox::getMeanFreePath(){
    //mean free path is (number density)*(scattering cross section)
    double numDens =sphereNum/(bounds[0]*bounds[1]*bounds[2]);
    double rad = (getRandomSphere())->getRadius();
    return numDens*(4*atan(1)*pow(rad,2));
}

// Determine if sphere box is at equilibrium
bool SphereBox::askEquil(){
	return 0;
}

// Return a random sphere
Sphere* SphereBox::getRandomSphere(){
    int randFactor;
    randFactor = (int) ((double) gen()/ gen.max())*sphereNum;
    return box[randFactor];
}

// Fill a box with randomly placed spheres, currently archaic, but might be worth saving
void SphereBox::fillBox(){
    bool spherePlaced;
    Sphere * place;
    Sphere * trialSphere;
    double rad;
    vector<Sphere*> overSpheres;
    // Fills the box with spheres
    for(int j = 0; j < sphereNum; j++)
    {
        place = box[j];
        rad = place->getRadius();
        trialSphere = new Sphere(rad);
        spherePlaced = false;
        bool overlap = false;
        vector<double> testCoords;

    // Attempts to assign a location to each sphere
        while(!spherePlaced)
        {
            // Generate a random set of coordinates
            testCoords = getRandomCoords();
            // Set trial sphere location
            trialSphere->setLoc(testCoords);
            double op = (trialSphere->getLoc())[0];
            // See if there is an overlap, if not, accept new location
            overSpheres = getNearSpheres(trialSphere);
            int overSize = 0;
            overSize = (int) overSpheres.size();
            if(overSize == 0)
            {
                numPlaced++;
                if(numPlaced % 100 == 0)
                {
                    cout << numPlaced << "\n";
                }
                box[j]->setLoc(testCoords);
                break;
            }

        }

    }
}

// Calculate the energy of adding a test sphere (for Widom insertion method)
double SphereBox::addTestSphere(Sphere* s){
    //initialize energy
    double energy;
    int nsize;
    Sphere* nSphere;
    // Set the test sphere to a random location
    s->setLoc(getRandomCoords());
    vector<Sphere*> nears = getNearSpheres(s);
    energy = getEnergySpheres(s,nears);
    return energy;
}

// Prints out the sphere data to the screen
void SphereBox::printSpheres(){
    cout << "\nSphere Number   radius  x   y   z \n";
    Sphere * printr;
    for(int i = 0; i < sphereNum; i++)
    {
        printr = box[i];
        cout << i+1 << "  " << printr->getRadius() << "    " << printr->getLoc()[0] << "    " << printr->getLoc()[1] << "    " << printr->getLoc()[2] <<"\n";
    }

}

// Print sphere data to a file
void SphereBox::printSpheres(string name){
    ofstream wfile;
	int numType;
	Sphere * first;
	first = box[0];
    wfile.open(name);
	// Make the header
    wfile << "Bounds" << bounds[0] << bounds[1] << bounds[2] << "\n";
    wfile << (*first).getHeader()<< "\n";
	numType = (*first).getType();
	wfile << numType << "\n";
	// Print the spheres
    Sphere * printr;
    for(int i = 0; i < sphereNum; i++)
    {
        printr = box[i];
        wfile << i+1 << (printr->getData()) <<"\n";
    }

}

//Read sphere data in from a file
void SphereBox::setSphereData(string fname){
    const int headerNum = 4;
     // Load the file name
    ifstream filename(fname);
    if(filename.is_open())
    {
        int lineCount;
        lineCount = 0;
        string line;
        Sphere * newS;
        int sphereType;
        while(getline(filename,line))
        {
            if(lineCount < headerNum)
            {
                //Switch case with sphere type
                newS = new Sphere();
                // Read in location data
                newS->readData(line);
            }
            else
            {
                // Switch-case statement with header types
                sphereType = 4;

            }
            lineCount++;
        }
    }
    else
    {
        cout << "Bad filename.";
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/* Class WellSphere 
	defines a sphere with a square well interaction
    Extension of class Sphere
*/
class WellSphere : public Sphere{
protected:
    double wellDepth;
    double wellWidth;
	const static int sphereType = 2;
public:
    // Default constructor
    WellSphere();
    // Well depth and width
    WellSphere(double, double);
    // Well depth and width, and also radius
    WellSphere(double, double, double);
    // All above, and also location
    WellSphere(double, double, double, vector<double>);
    // New potential is dependent on wellDepth
    double getPotential(Sphere* s, vector<double>);
    // New interaction length is based on wellWidth
    double getInteractionLength();
};

/* The default square well sphere is a hard sphere.
	This will be used for testing purposes.
*/
WellSphere::WellSphere(){
    wellWidth = 1;
    wellDepth = 0;
}

// Constructor that defines well depth and well width
WellSphere::WellSphere(double depth, double width){
    wellDepth = depth;
    wellWidth = width;
}

// Constructor that defines well depth and well width and also radius
WellSphere::WellSphere(double depth, double width, double rad){
    wellDepth = depth;
    wellWidth = width;
    radius = rad;
}

// Constructor that defines well depth and well width and also radius and even location!
WellSphere::WellSphere(double depth, double width, double rad, vector<double> loc){
    wellDepth = depth;
    wellWidth = width;
    radius = rad;
    location = loc;
}

// Overloaded potential function
double WellSphere::getPotential(Sphere* s, vector<double> bounds){
    double minDist;
    double dist;
    // Calculate the minimum allowable distance between the spheres
    minDist = radius + s->getRadius();
    dist = getDistance(s,bounds);
    // If the distance is less than allowed, set potential to infinity
    if (dist <= minDist){
        return INFINITY;
    }
    // Otherwise default is hard sphere
    else if (dist <= wellWidth*minDist)
    {
        return wellDepth;
    }
    else
    {
        return 0;
    }
}

// Overloaded interaction length
double WellSphere::getInteractionLength(){
    return 2*wellWidth*radius;
}

/*////////////////////////////////////////////MAIN FUNCTIONS/////////////////////////////////////////////////

	1. Create a box of spheres of some defined type
	2. Use a Metropolis Markov Chain Monte Carlo simulation to calculate the partition function
	3. Use a Widom insertion method simulation to calculate the chemical potential
	
*////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Make a list of spheres, save to file.
void makeSphereList(){
    // Input parameters
    cout << "\nInput the number of spheres:\n";
    int numSpheres;
    cin >> numSpheres;
    cout << "\nInput the sphere radius\n";
    double inputRadius;
    cin >> inputRadius;
    cout << "\nInput the side length for the box\n";
    double boxSize;
    cin >> boxSize;
	
	// Calculate grid spacing
	int sStep;
	double kStep;
	double shiftVal;
	sStep = ceil(cbrt((double) numSpheres));
	kStep = boxSize / ((double) sStep + 1);
	shiftVal = ((double) sStep+.5) * kStep;

    // Reject impossible cases
    double mpf = pow((2*(inputRadius) / boxSize),3)*numSpheres;
    if(mpf >= 1)
    {
        cout << "\nPacking fraction too high for this sort of simulation";
        return;
    }
	
	// Calculate actual packing fraction
	double pf = (4 * PI * pow(inputRadius/boxSize,3) / 3);
	cout << "Creating box with packing fraction =" << pf << "\n";

	// Make a header
	int inputType;
    cout << "Enter the type of sphere to create:\n1. Hard Sphere\n2. Square Well\n";
    cin >> inputType;
	
	//Do a different thing based on the type of sphere used
	Sphere* s;
	vector<Sphere*> sphereList;
	vector<double> locScanner (3, kStep);
	switch (inputType){
		// hard sphere
		case 1: 
		// do the loop thing
		for (int sphereAdd = 0; sphereAdd < numSpheres; sphereAdd++){
			// Make a new sphere with the input radius and location based on locScanner
			s = new Sphere(inputRadius, locScanner);
			sphereList.push_back(s);
			
			// Increment locScanner
			locScanner[0] += kStep;
			
			// If overflow, move to a new row
			if(locScanner[0] > shiftVal){
				locScanner[0] = kStep;
				locScanner[1] += kStep;
				// If overflow, move to a new slice
				if(locScanner[1] > shiftVal){
					locScanner[1] = kStep;
					locScanner[2] += kStep;
				}
			}
			cout << "iteration " << sphereAdd+1 << " completed\n";
		}
		break;
		
		// square well
		case 2:
		// prompt width and depth
		double width, depth;
		width = 1;
		depth = 0;
		cout << "Enter width and depth:";
		for (int sphereAdd = 0; sphereAdd < numSpheres; sphereAdd++){
			// Make a new sphere with the input radius and location based on locScanner
			s = new WellSphere(width, depth, inputRadius, locScanner);
			sphereList.push_back(s);
			// Increment locScanner
			locScanner[1] += kStep;
			// If overflow, move to a new row
			if(locScanner[1] > shiftVal){
				locScanner[1] = kStep;
				locScanner[2] += kStep;
				// If overflow, move to a new slice
				if(locScanner[2] > shiftVal){
					locScanner[2] = kStep;
					locScanner[3] += kStep;
				}
			}
		}
		// do the loop thing
		break;
		
		default:
		cout << "Bad input type";
		break;
	}
	
	cout << "Calling constructor";

    // Call constructor with input parameters, which automatically fills box.
    SphereBox boxOfSpheres (sphereList,boxSize,boxSize,boxSize);
	
	cout << "printing spheres";
	//boxOfSpheres.printSpheres();

    // Output the data to an input file
    string str;
    cout << "Enter a filename to write the sphere data to:\n";
    cin >> str;
    // Try to write to file
    try
    {
        boxOfSpheres.printSpheres(str);
        cout << "File written.\n";
    } // If that fails for whatever reason, print the box to screen.
    catch(int e)
    {
        boxOfSpheres.printSpheres();
        cout << "File write failed, output printed to screen.\n";
    }
	//
}

//Read in a list of spheres, and then perturb them
void metroMonteCarlo(){
    // Construct a sphere box from a pre-made file
    SphereBox box;
    Sphere * pert;
    Sphere * placeHold;
    bool equilibrium;
    string fname;
    double oldEnergy,newEnergy,energyRatio,rfactor,pertFactor,gMax;

    // Read in sphere data from file
    cout << "\nEnter the file name where the sphere data is saved.\n";
    cin >> fname;
    box.setSphereData(fname);

    // Calculate mean free path for Metropolis Monte Carlo
    pertFactor = box.getMeanFreePath();

    // Allow the box to come to equilibrium using boltzmann interaction
    while(!equilibrium){
        // Select a sphere from the box at random
        pert = box.getRandomSphere();
        // Get energy of the sphere
        oldEnergy = box.getEnergySpheres(pert, box.getNearSpheres(pert));
        // Save a copy of the sphere
        *placeHold = *pert;
        // Perturb the sphere
        pert->perturb(pertFactor);
        // Calculate the new energy
        newEnergy = box.getEnergySpheres(pert, box.getNearSpheres(pert));
        // Compare the two boltzmann factors
        energyRatio = exp(-1*newEnergy)/exp(-1*oldEnergy);
        // If new energy is unfavored, only keeps by random chance
        if(energyRatio < 1)
        {
            //generate a low resolution random number between 0 and 1
            rfactor = gen()/gMax;
            // if that number is higher than the energy ratio, reject the change
            if(rfactor > energyRatio)
            {
                // If failed, overwrite with old energy again
                *pert = *placeHold;
            }
        }
		// Ask the box if it's at equilibrium
		equilibrium = box.askEquil();
    }

    // Ask if new box should be written to file
    char input;
    cout << "Save box to file? (Y/N)";
    cin >> input;
    if (input == 'y')
    {
        cout << "Enter File to save new data to";
    }
}

// Read in a list of pre-made spheres, set the type of sphere, perform a simulation
void widomPotential(){
    double temp;
    double eTot,eAvg, eTest;
    double chemPot;
    // Request the filename where the sphere data is written
    string name;
    cout << "\nEnter the name of the file where the data is saved:\n";
    cin >> name;
    SphereBox box;
    box.setSphereData(name);
    // Inquire about simulation details (i.e. number of spheres to test with)
    int numTests;
    cout << "\nEnter the number of spheres to test with.\n";
    cin >> numTests;
    // Perform simulation
    eTot = 0;
    for(int num = 0; num < numTests; num++)
    {
        eTest = box.addTestSphere(new Sphere());
        eTot += exp(-1*eTest/temp);
    }
    eAvg = eTot/numTests;
    chemPot = -temp*log(eAvg);
    cout << "\nThe Chemical Potential of the input box is approximated to be:" << chemPot << "\n";

}

// Main
int main(){
	gen.seed(time(NULL));
    bool done = false;
    while(!done)
    {
        cout << "\nInput a number to perform a task:\n(1):Make a box of spheres.\n(2):Perform a metropolis Monte-Carlo simulation.\n(3):Perform Widom insertion method.\n(0):Quit\n";
        int input;
        cin >> input;
        switch (input)
        {
            case 1: makeSphereList();
            break;
            case 2: metroMonteCarlo();
            break;
            case 3: widomPotential();
            break;
            case 0: done = true;
            break;
            default: cout << "\nBad input.";
            break;
        }
    }
    return 0;
}
