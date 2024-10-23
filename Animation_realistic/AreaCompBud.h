#ifndef AREACOMPBUD_H_
#define AREACOMPBUD_H_ 

#include "SystemStructures.h"

// Function declaration for computing the area of a bud.
void ComputeAreaBud(
    
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs//,
    
    );

// Functor struct for computing the area of a bud.    
struct AreaBudCompFunctor {
    
    // Member variables
    // double spring_constant;  // Not used in this code snippet
    double* locXAddr; // Pointer to the X-coordinate array of node positions
    double* locYAddr; // Pointer to the Y-coordinate array of node positions
    double* locZAddr; // Pointer to the Z-coordinate array of node positions
    int* triangles_in_upperhem; // Pointer to the array indicating whether each triangle is in the upper hemisphere

    // Constructor
	  __host__ __device__ AreaBudCompFunctor(
        
        // double& _spring_constant,  // Not used in this code snippet
        double* _locXAddr,
        double* _locYAddr,
        double* _locZAddr,
        int* _triangles_in_upperhem
        
        ):

    // spring_constant(_spring_constant),  // Not used in this code snippet
    locXAddr(_locXAddr),
    locYAddr(_locYAddr),
    locZAddr(_locZAddr),
    triangles_in_upperhem(_triangles_in_upperhem){}

  	// Operator () Overload
    // This function is called for each element in the input iterator range.
    // hand in counting iterator and id of two edges and preferred length
  	__device__ double operator()(const Tuuuu &u4) {
        		
        // Extract the tuple elements
        //counter ranges from 0 to num_edges. 
        int counter = thrust::get<0>(u4); // Counter ranging from 0 to num_edges
        int r1 = thrust::get<1>(u4); // Node index 1 of the triangle
        int r2 = thrust::get<2>(u4); // Node index 2 of the triangle
        int r3 = thrust::get<3>(u4); // Node index 3 of the triangle
        
        // Check if the triangle is in the upper hemisphere and the node indices are valid
        if (triangles_in_upperhem[counter] == 1 && r1 != INT_MAX && r2 != INT_MAX && r3 != INT_MAX){
          
            // Compute the area of the triangle        
            double r1x = locXAddr[r1];
            double r1y = locYAddr[r1];
            double r1z = locZAddr[r1];
            double r2x = locXAddr[r2];
            double r2y = locYAddr[r2];
            double r2z = locZAddr[r2];
            double r3x = locXAddr[r3];
            double r3y = locYAddr[r3];
            double r3z = locZAddr[r3];
            double norm_r1r2 = sqrt((r2x-r1x)*(r2x-r1x) + (r2y-r1y)*(r2y-r1y) + (r2z-r1z)*(r2z-r1z));
            double norm_r2r3 = sqrt((r3x-r2x)*(r3x-r2x) + (r3y-r2y)*(r3y-r2y) + (r3z-r2z)*(r3z-r2z));
            double norm_r3r1 = sqrt((r3x-r1x)*(r3x-r1x) + (r3y-r1y)*(r3y-r1y) + (r3z-r1z)*(r3z-r1z));
            double s = (norm_r1r2 + norm_r2r3 + norm_r3r1)/2.0;
            double area = sqrt(s*(s-norm_r1r2)*(s-norm_r2r3)*(s-norm_r3r1));
            
            return area;
        
        }
        
        else{
            
            double area = 0.0;
            
            return area;
        }


    }
};

#endif