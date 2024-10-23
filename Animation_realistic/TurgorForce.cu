
#include "System.h"
#include "SystemStructures.h"
#include "TurgorForce.h"

// Function to compute turgor pressure forces.
void ComputeTurgorSprings(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    AreaTriangleInfoVecs& areaTriangleInfoVecs) {

    // Create a counting iterator over element IDs (triangles).
    thrust::counting_iterator<int> elemId(0);
    
    // Initialize vectors to store the turgor pressure forces.
    thrust::fill(areaTriangleInfoVecs.tempNodeForceXReduced.begin(), areaTriangleInfoVecs.tempNodeForceXReduced.end(), 0.0);
    thrust::fill(areaTriangleInfoVecs.tempNodeForceYReduced.begin(), areaTriangleInfoVecs.tempNodeForceYReduced.end(), 0.0);
    thrust::fill(areaTriangleInfoVecs.tempNodeForceZReduced.begin(), areaTriangleInfoVecs.tempNodeForceZReduced.end(), 0.0);
    thrust::fill(areaTriangleInfoVecs.tempNodeForceXUnreduced.begin(), areaTriangleInfoVecs.tempNodeForceXUnreduced.end(), 0.0);
    thrust::fill(areaTriangleInfoVecs.tempNodeForceYUnreduced.begin(), areaTriangleInfoVecs.tempNodeForceYUnreduced.end(), 0.0);
    thrust::fill(areaTriangleInfoVecs.tempNodeForceZUnreduced.begin(), areaTriangleInfoVecs.tempNodeForceZUnreduced.end(), 0.0);

    // Calculate the un-reduced turgor pressure forces acting on each triangle.
    areaTriangleInfoVecs.dummy = thrust::transform_reduce(
        thrust::make_zip_iterator(
            thrust::make_tuple(
                elemId,
                coordInfoVecs.triangles2Nodes_1.begin(),
                coordInfoVecs.triangles2Nodes_2.begin(),
                coordInfoVecs.triangles2Nodes_3.begin(),
                coordInfoVecs.triangles2Edges_1.begin(),
                coordInfoVecs.triangles2Edges_2.begin(),
                coordInfoVecs.triangles2Edges_3.begin())),
        thrust::make_zip_iterator(
            thrust::make_tuple(
                elemId,
                coordInfoVecs.triangles2Nodes_1.begin(),
                coordInfoVecs.triangles2Nodes_2.begin(),
                coordInfoVecs.triangles2Nodes_3.begin(),
                coordInfoVecs.triangles2Edges_1.begin(),
                coordInfoVecs.triangles2Edges_2.begin(),
                coordInfoVecs.triangles2Edges_3.begin())) + coordInfoVecs.num_triangles,
        TurgorSpringFunctor(
            generalParams.volume_spring_constant,
            thrust::raw_pointer_cast(coordInfoVecs.nodeLocX.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeLocY.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeLocZ.data()),
            thrust::raw_pointer_cast(areaTriangleInfoVecs.tempNodeIdUnreduced.data()),
            thrust::raw_pointer_cast(areaTriangleInfoVecs.tempNodeForceXUnreduced.data()),
            thrust::raw_pointer_cast(areaTriangleInfoVecs.tempNodeForceYUnreduced.data()),
            thrust::raw_pointer_cast(areaTriangleInfoVecs.tempNodeForceZUnreduced.data())),
        0.0, thrust::plus<double>());

    // Sort the un-reduced forces by node ID.
    thrust::sort_by_key(areaTriangleInfoVecs.tempNodeIdUnreduced.begin(), 
                        areaTriangleInfoVecs.tempNodeIdUnreduced.begin() + (areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles),
                        thrust::make_zip_iterator(
                            thrust::make_tuple(
                                areaTriangleInfoVecs.tempNodeForceXUnreduced.begin(),
                                areaTriangleInfoVecs.tempNodeForceYUnreduced.begin(),
                                areaTriangleInfoVecs.tempNodeForceZUnreduced.begin())),
                        thrust::less<int>());

    // Reduce the un-reduced forces by node ID and store the result in reduced vectors.
    int endKey = thrust::get<0>(
        thrust::reduce_by_key(
            areaTriangleInfoVecs.tempNodeIdUnreduced.begin(), 
            areaTriangleInfoVecs.tempNodeIdUnreduced.begin() + (areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles),
            thrust::make_zip_iterator(
                thrust::make_tuple(
                    areaTriangleInfoVecs.tempNodeForceXUnreduced.begin(),
                    areaTriangleInfoVecs.tempNodeForceYUnreduced.begin(),
                    areaTriangleInfoVecs.tempNodeForceZUnreduced.begin())),
            areaTriangleInfoVecs.tempNodeIdReduced.begin(),
            thrust::make_zip_iterator(
                thrust::make_tuple(
                    areaTriangleInfoVecs.tempNodeForceXReduced.begin(),
                    areaTriangleInfoVecs.tempNodeForceYReduced.begin(),
                    areaTriangleInfoVecs.tempNodeForceZReduced.begin())),
            thrust::equal_to<int>(), CVec3Add())) - areaTriangleInfoVecs.tempNodeIdReduced.begin();

    // Apply the reduced turgor forces to all nodes.
    thrust::for_each(
        thrust::make_zip_iterator(
            thrust::make_tuple(
                areaTriangleInfoVecs.tempNodeIdReduced.begin(),
                areaTriangleInfoVecs.tempNodeForceXReduced.begin(),
                areaTriangleInfoVecs.tempNodeForceYReduced.begin(),
                areaTriangleInfoVecs.tempNodeForceZReduced.begin())),
        thrust::make_zip_iterator(
            thrust::make_tuple(
                areaTriangleInfoVecs.tempNodeIdReduced.begin(),
                areaTriangleInfoVecs.tempNodeForceXReduced.begin(),
                areaTriangleInfoVecs.tempNodeForceYReduced.begin(),
                areaTriangleInfoVecs.tempNodeForceZReduced.begin())) + endKey,
        AddForceFunctor (
            thrust::raw_pointer_cast(coordInfoVecs.nodeForceX.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeForceY.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeForceZ.data())));
            
   
};
