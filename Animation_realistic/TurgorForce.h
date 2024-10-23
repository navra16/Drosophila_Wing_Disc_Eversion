
#ifndef TURGORFORCE_H_
#define TURGORFORCE_H_

#include "SystemStructures.h"

// Function declaration to compute turgor pressure forces.
void ComputeTurgorSprings(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    AreaTriangleInfoVecs& areaTriangleInfoVecs);

// Functor to calculate turgor pressure forces.
struct TurgorSpringFunctor {
    double spring_constant; // Spring constant for the turgor pressure force.
    double* locXAddr;       // Pointer to the X-coordinates of the nodes.
    double* locYAddr;       // Pointer to the Y-coordinates of the nodes.
    double* locZAddr;       // Pointer to the Z-coordinates of the nodes.

    int* idKey;             // Pointer to store the node IDs.
    double* forceXAddr;     // Pointer to store the X-component of the turgor pressure force on nodes.
    double* forceYAddr;     // Pointer to store the Y-component of the turgor pressure force on nodes.
    double* forceZAddr;     // Pointer to store the Z-component of the turgor pressure force on nodes.

    // Constructor for the TurgorSpringFunctor.
    __host__ __device__ TurgorSpringFunctor(
        double& _spring_constant,   // Reference to the spring constant for turgor force.
        double* _locXAddr,          // Pointer to X-coordinates of nodes.
        double* _locYAddr,          // Pointer to Y-coordinates of nodes.
        double* _locZAddr,          // Pointer to Z-coordinates of nodes.
        int* _idKey,                // Pointer to store the node IDs.
        double* _forceXAddr,        // Pointer to store X-component of turgor force on nodes.
        double* _forceYAddr,        // Pointer to store Y-component of turgor force on nodes.
        double* _forceZAddr         // Pointer to store Z-component of turgor force on nodes.
    ) :
        spring_constant(_spring_constant),
        locXAddr(_locXAddr),
        locYAddr(_locYAddr),
        locZAddr(_locZAddr),
        idKey(_idKey),
        forceXAddr(_forceXAddr),
        forceYAddr(_forceYAddr),
        forceZAddr(_forceZAddr) {}

    // Functor operator that calculates turgor pressure forces on the nodes.
    __device__ double operator()(const Tuuuuuuu &u7) {
        // Test placing the ids of the nodes and then get positions.
        // double scaling_pow = 4.0;
        int counter = thrust::get<0>(u7);
        int place = 3 * counter; // Represents location in write to vector.

        int id_i = thrust::get<1>(u7); // Node ID of the first vertex of the triangle.
        int id_j = thrust::get<2>(u7); // Node ID of the second vertex of the triangle.
        int id_k = thrust::get<3>(u7); // Node ID of the third vertex of the triangle.
        int e_id_i = thrust::get<4>(u7); // Edge ID of the first edge of the triangle.
        int e_id_j = thrust::get<5>(u7); // Edge ID of the second edge of the triangle.
        int e_id_k = thrust::get<6>(u7); // Edge ID of the third edge of the triangle.

        // If the node IDs are valid (not INT_MAX) or negative, proceed with the calculations.
        if ((id_i != INT_MAX || id_i < 0) && (id_j != INT_MAX || id_j < 0) && (id_k != INT_MAX || id_k < 0)) {
            // Extract node positions.
            CVec3 ri = thrust::make_tuple<double>(locXAddr[id_i], locYAddr[id_i], locZAddr[id_i]);
            CVec3 rj = thrust::make_tuple<double>(locXAddr[id_j], locYAddr[id_j], locZAddr[id_j]);
            CVec3 rk = thrust::make_tuple<double>(locXAddr[id_k], locYAddr[id_k], locZAddr[id_k]);

            // Compute vectors between nodes.
            CVec3 rkj = CVec3_subtract(rk, rj);
            CVec3 rij = CVec3_subtract(ri, rj);

            // Calculate the area of the triangle formed by the three nodes.
            double area_current = sqrt(CVec3_dot(CVec3_cross(rkj, rij), CVec3_cross(rkj, rij))) / 2.0;

            // Compute derivatives with respect to area.
            CVec3 A = CVec3_cross(rkj, rij); // rkj must come first
            double A1 = thrust::get<0>(A);
            double A2 = thrust::get<1>(A);
            double A3 = thrust::get<2>(A);
            double magnitude = sqrt(A1 * A1 + A2 * A2 + A3 * A3);
            A1 = A1 / magnitude;
            A2 = A2 / magnitude;
            A3 = A3 / magnitude;
            
            // Calculate the turgor pressure forces acting on the nodes.
            CVec3 rj_force = thrust::make_tuple<double>((1.0 / 3.0) * area_current * A1 * spring_constant,
                                                        (1.0 / 3.0) * area_current * A2 * spring_constant,
                                                        (1.0 / 3.0) * area_current * A3 * spring_constant);
                
            CVec3 rk_force = thrust::make_tuple<double>((1.0 / 3.0) * area_current * A1 * spring_constant,
                                                        (1.0 / 3.0) * area_current * A2 * spring_constant,
                                                        (1.0 / 3.0) * area_current * A3 * spring_constant);

            CVec3 ri_force = thrust::make_tuple<double>((1.0 / 3.0) * area_current * A1 * spring_constant,
                                                        (1.0 / 3.0) * area_current * A2 * spring_constant,
                                                        (1.0 / 3.0) * area_current * A3 * spring_constant);

            // Store the computed forces and node IDs in the appropriate arrays.
            idKey[place] = id_i;
            forceXAddr[place] = thrust::get<0>(ri_force);
            forceYAddr[place] = thrust::get<1>(ri_force);
            forceZAddr[place] = thrust::get<2>(ri_force);

            idKey[place + 1] = id_j;
            forceXAddr[place + 1] = thrust::get<0>(rj_force);
            forceYAddr[place + 1] = thrust::get<1>(rj_force);
            forceZAddr[place + 1] = thrust::get<2>(rj_force);
            
            idKey[place + 2] = id_k;
            forceXAddr[place + 2] = thrust::get<0>(rk_force);
            forceYAddr[place + 2] = thrust::get<1>(rk_force);
            forceZAddr[place + 2] = thrust::get<2>(rk_force);

            // double energy = what_spring_constant / (2.0) * (area_current - area_0) * (area_current - area_0) / area_0;
            // return energy;
        }
        else {
            // double energy = 0.0;
            // return energy;
        }
    }
};
#endif // TURGORFORCE_H_
