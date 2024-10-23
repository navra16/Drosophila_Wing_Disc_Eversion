
#include "System.h"
#include "SystemStructures.h"
#include "BendingTrianglesEnergy.h"

double ComputeCosTrianglesEnergy(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    BendingTriangleInfoVecs& bendingTriangleInfoVecs) {
	
    
    thrust::counting_iterator<int> elemId(0); 

	  // bendingTriangleInfoVecs.initial_angle = 1.5707963267/2.0;

    // Apply force to temporary vectors.
    bendingTriangleInfoVecs.bending_triangle_energy= 
    thrust::transform_reduce(
        thrust::make_zip_iterator(
            thrust::make_tuple(
				        elemId,
                coordInfoVecs.edges2Triangles_1.begin(),
                coordInfoVecs.edges2Triangles_2.begin(),
                coordInfoVecs.edges2Nodes_1.begin(),
                coordInfoVecs.edges2Nodes_2.begin())),
        thrust::make_zip_iterator(
            thrust::make_tuple(
				        elemId, 
                coordInfoVecs.edges2Triangles_1.begin(),
                coordInfoVecs.edges2Triangles_2.begin(),
                coordInfoVecs.edges2Nodes_1.begin(),
                coordInfoVecs.edges2Nodes_2.begin())) + coordInfoVecs.num_edges,
        CosBendingEnergyFunctor(
            bendingTriangleInfoVecs.spring_constant,
            bendingTriangleInfoVecs.initial_angle,        
            thrust::raw_pointer_cast(coordInfoVecs.nodeLocX.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeLocY.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeLocZ.data()),
            
            thrust::raw_pointer_cast(bendingTriangleInfoVecs.tempNodeIdUnreduced.data()),
            thrust::raw_pointer_cast(coordInfoVecs.triangles2Nodes_1.data()), 
            thrust::raw_pointer_cast(coordInfoVecs.triangles2Nodes_2.data()), 
            thrust::raw_pointer_cast(coordInfoVecs.triangles2Nodes_3.data())),
		    0.0, thrust::plus<double>() );
	  
    // Energy value for bending triangles returned. 
    // This could be the cosine bending potential for the Foppl-von Karman number.
    return bendingTriangleInfoVecs.bending_triangle_energy;
};

// Triangles are attached to each other using springs. 
// Since the parent cell as well as the bud are spherical or partially spherical, the triangles adjacent to each other would possess bending energy in their springs. 
