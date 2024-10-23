/*
 * SystemBuilder.h
 *
 *  Created on: 25 авг. 2014 г.
 *      Author: yan
 */

#ifndef SystemBuilder_H_
#define SystemBuilder_H_

#include "SystemStructures.h" // Include the header containing various data structures and functors.
#include "System.h" // Include the System class.

// SystemBuilder class is responsible for building the system by adding nodes, edges, elements, and fixing nodes.
class SystemBuilder {
public:
    // Variables set by the constructor using command-line arguments.
    double dt; // Time step for simulation.
    int solve_time; // Total simulation time.

    // Default values for various physical constants (can be overridden by an XML input file).
    double defaultTau = 1.0; // Default tau (time constant) value.
    double defaultKBT = 1.0; // Default kBT (Boltzmann constant times temperature) value.
    double defaultLinear_Const = 9.0; // Default linear constant for edge springs.
    double defaultArea_Const = 10.0; // Default area constant for area springs.
    double defaultBending_Const = 4.0; // Default bending constant for bending springs.
    double defaultLJ_Eps = 0.1; // Default Lennard-Jones epsilon value.
    double defaultLJ_Rmin = 2.0; // Default Lennard-Jones minimum distance value.
    double defaultLJ_Rmax = 2.0 * 1.4; // Default Lennard-Jones maximum distance value.
    double defaultLJ_Const = 1.0; // Default Lennard-Jones constant.
    double defaultLJ_X = 0.0; // Default Lennard-Jones X value.
    double defaultLJ_Y = 0.0; // Default Lennard-Jones Y value.
    double defaultLJ_Z = -0.1; // Default Lennard-Jones Z value.

    double defaultEdgeEq = 1.0; // Default equilibrium length for edges (can be overridden by input).
    double defaultAreaEq = 0.0; // Default equilibrium area for triangles (can be overridden by input).
    double defaultAngleEq = 0.0; // Default equilibrium angle for bending triangles (can be overridden by input).

    HostSetInfoVecs hostSetInfoVecs; // Collection of host vectors for setting up the system.

public:
    // Constructor that takes the timestep and solve time as arguments.
    SystemBuilder(double timestep, int solve_time);

    // Function to add a node to the system with the specified coordinates (x, y, z).
    void addNode(double x, double y, double z);

    // Function to add non-dimensional node data (x1 - x9) for the node. This is likely related to a specific simulation case.
    void addNndata(double x1, double x2, double x3, double x4, double x5, double x6, double x7, double x8, double x9);

    // Function to add an edge between two nodes with the specified IDs (idL and idR).
    void addEdge(int idL, int idR);

    // Function to add an edge between two nodes with the specified IDs (idL and idR) and initial length (edge_initial_length).
    void addEdge(int idL, int idR, double edge_initial_length);

    // Function to add an element (triangle) to the system using three node IDs (idA, idB, and idC).
    void addElement(int idA, int idB, int idC);

    // Function to add an element (triangle) to the system using three edge IDs (idA, idB, and idC).
    void addElement2Edge(int idA, int idB, int idC);

    // Function to add an edge to element mapping by providing element ID (idA) and two edge IDs (idB and idC).
    void addEdge2Elem(int idA, int idB);

    // Function to fix a node with the specified ID (id).
    void fixNodes(int id);

    // Function to add a capsid node to the system with the specified coordinates (x, y, z).
    void addCapsidNode(double x, double y, double z);

    // Function to create the system and return a shared pointer to it.
    std::shared_ptr<System> createSystem();
};

#endif /* SystemBuilder_H_ */
//The "SystemBuilder.h" header file defines the `SystemBuilder` class, which is responsible for building a system by adding nodes, edges, elements, and fixing nodes. It also provides default values for various physical constants and functions to create the system.


//
///*
// * SystemBuilder.h
// *
// *  Created on: 25 авг. 2014 г.
// *      Author: yan
// */
