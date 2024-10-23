
#include "System.h"
#include "SystemStructures.h"
#include "BendingTriangles.h"

// Define a function named ComputeCosTriangleSprings.
// The function takes references to objects of types GeneralParams, CoordInfoVecs, and BendingTriangleInfoVecs as arguments and returns void. There is no output for this?
void ComputeCosTriangleSprings(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    BendingTriangleInfoVecs& bendingTriangleInfoVecs) {
	
    // Create a counting iterator named elemID initialized with 0.
    thrust::counting_iterator<int> elemId(0); 
    
    //Fill the tempNodeForceXReduced, tempNodeForceYReduced, tempNodeForceZReduced, tempNodeForceXUnreduced, tempNodeForceYUnreduced and tempNodeForceZUnreduced vectors with zeros. 
  	//bendingTriangleInfoVecs.initial_angle = 1.5707963267/2.0;
  	thrust::fill(bendingTriangleInfoVecs.tempNodeForceXReduced.begin(),bendingTriangleInfoVecs.tempNodeForceXReduced.end(),0.0);
  	thrust::fill(bendingTriangleInfoVecs.tempNodeForceYReduced.begin(),bendingTriangleInfoVecs.tempNodeForceYReduced.end(),0.0);
  	thrust::fill(bendingTriangleInfoVecs.tempNodeForceZReduced.begin(),bendingTriangleInfoVecs.tempNodeForceZReduced.end(),0.0);
  	thrust::fill(bendingTriangleInfoVecs.tempNodeForceXUnreduced.begin(),bendingTriangleInfoVecs.tempNodeForceXUnreduced.end(),0.0);
  	thrust::fill(bendingTriangleInfoVecs.tempNodeForceYUnreduced.begin(),bendingTriangleInfoVecs.tempNodeForceYUnreduced.end(),0.0);
  	thrust::fill(bendingTriangleInfoVecs.tempNodeForceZUnreduced.begin(),bendingTriangleInfoVecs.tempNodeForceZUnreduced.end(),0.0);

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
        //This is used to calculate the bending triangle energy. 
        CosBendingFunctor(
            generalParams.SCALE_TYPE,
            generalParams.nonuniform_wall_weakening_bend,
            generalParams.maxSpringScaler_bend,
            generalParams.scaling_pow,
            generalParams.gausssigma,
            generalParams.hilleqnconst,
            generalParams.hilleqnpow,
            thrust::raw_pointer_cast(coordInfoVecs.scaling_per_edge.data()),
            bendingTriangleInfoVecs.spring_constant,
            bendingTriangleInfoVecs.spring_constant_weak,
            thrust::raw_pointer_cast(generalParams.edges_in_upperhem.data()),
            thrust::raw_pointer_cast(generalParams.boundaries_in_upperhem.data()),
            bendingTriangleInfoVecs.initial_angle,
            bendingTriangleInfoVecs.initial_angle_bud,        
            thrust::raw_pointer_cast(coordInfoVecs.nodeLocX.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeLocY.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeLocZ.data()),
            
            thrust::raw_pointer_cast(bendingTriangleInfoVecs.tempNodeIdUnreduced.data()),
            thrust::raw_pointer_cast(bendingTriangleInfoVecs.tempNodeForceXUnreduced.data()),
            thrust::raw_pointer_cast(bendingTriangleInfoVecs.tempNodeForceYUnreduced.data()),
            thrust::raw_pointer_cast(bendingTriangleInfoVecs.tempNodeForceZUnreduced.data()),
            thrust::raw_pointer_cast(coordInfoVecs.triangles2Nodes_1.data()), 
            thrust::raw_pointer_cast(coordInfoVecs.triangles2Nodes_2.data()), 
            thrust::raw_pointer_cast(coordInfoVecs.triangles2Nodes_3.data())),
		    0.0, thrust::plus<double>() );
	
    // Calculate the bending triangle energy using the CosBendingFunctor and store the result in bendingTrianglesInfoVecs.bending_triangle_energy.
    //thrust::sort_by_key function is called to sort the temporary forced by node ID. 
    //now we have un reduced forces. Sort by id and reduce. 
    //key, then value. Each vector returns sorted		
    thrust::sort_by_key(bendingTriangleInfoVecs.tempNodeIdUnreduced.begin(), bendingTriangleInfoVecs.tempNodeIdUnreduced.begin() + (bendingTriangleInfoVecs.factor*coordInfoVecs.num_edges),
        thrust::make_zip_iterator(
            thrust::make_tuple(
                bendingTriangleInfoVecs.tempNodeForceXUnreduced.begin(),
                bendingTriangleInfoVecs.tempNodeForceYUnreduced.begin(),
                bendingTriangleInfoVecs.tempNodeForceZUnreduced.begin())), thrust::less<int>());

    // Calculate the endKey vale by reducing the forces by key (node ID). 
    int endKey = thrust::get<0>(
        thrust::reduce_by_key(
            bendingTriangleInfoVecs.tempNodeIdUnreduced.begin(), 
            bendingTriangleInfoVecs.tempNodeIdUnreduced.begin() + (bendingTriangleInfoVecs.factor*coordInfoVecs.num_edges),
        thrust::make_zip_iterator(
            thrust::make_tuple(
                bendingTriangleInfoVecs.tempNodeForceXUnreduced.begin(),
                bendingTriangleInfoVecs.tempNodeForceYUnreduced.begin(),
                bendingTriangleInfoVecs.tempNodeForceZUnreduced.begin())),
			
            bendingTriangleInfoVecs.tempNodeIdReduced.begin(),
        thrust::make_zip_iterator(
            thrust::make_tuple(
                bendingTriangleInfoVecs.tempNodeForceXReduced.begin(),
                bendingTriangleInfoVecs.tempNodeForceYReduced.begin(),
                bendingTriangleInfoVecs.tempNodeForceZReduced.begin())),
		    thrust::equal_to<int>(), CVec3Add())) - bendingTriangleInfoVecs.tempNodeIdReduced.begin();//binary_pred, binary_op 
		
     // The reduced forces are obtained and stored in bendingTriangleInfoVecs.tempNodeForceXReduced, bendingTriangleInfoVecs,tempNodeForceYReduced, bendingTriangleInfoVecs,tempNodeForceZReduced by reducing the forcces by key (node ID).
     
     // Apply the reduced forces to all nodes. 
    thrust::for_each(
        thrust::make_zip_iterator(//1st begin
            thrust::make_tuple(
                bendingTriangleInfoVecs.tempNodeIdReduced.begin(),
                bendingTriangleInfoVecs.tempNodeForceXReduced.begin(),
                bendingTriangleInfoVecs.tempNodeForceYReduced.begin(),
                bendingTriangleInfoVecs.tempNodeForceZReduced.begin())),
        thrust::make_zip_iterator(//1st end
            thrust::make_tuple(
                bendingTriangleInfoVecs.tempNodeIdReduced.begin(),
                bendingTriangleInfoVecs.tempNodeForceXReduced.begin(),
                bendingTriangleInfoVecs.tempNodeForceYReduced.begin(),
                bendingTriangleInfoVecs.tempNodeForceZReduced.begin())) + endKey,
        AddForceFunctor (
            thrust::raw_pointer_cast(coordInfoVecs.nodeForceX.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeForceY.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeForceZ.data())));


};


