#include "boundary.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"

void treatBoundary(double *collideField, int* flagField, const double * const wallVelocity, int xlength){
    int x, y, z, i, index;
    int dim = xlength + 2;
    int dim_pow2 = dim * dim;
    int xlen_pow2 = (dim - 1) * dim;
    int xlen_pow3 = xlen_pow2 * dim;
    int Iindex[5]; // stores i numbers
    double cs_2 = C_S * C_S;

    /******************/
  	/** bottom, z = 0 */
  	/******************/
  	// interior boundary cells
 	for( x = 2; x < xlength; x++ )
 		for( y = 2; y < xlength; y++ )
 			for( i = 14; i < 19; i++ ){

 			     collideField[ 19 * (y * dim + x) + i] = 
 			                 collideField[19 * (LATTICEVELOCITIES[i][2] * dim_pow2 + 
 			           	                  (y + LATTICEVELOCITIES[i][1]) * dim + x + LATTICEVELOCITIES[i][0]) + 18 - i];
 		    }
    // z = 0, y = 1 & y = 0
    for( x = 2; x < xlength; x++ )
    // y = 1
		for( i = 15; i < 19; i++ ){

		     collideField[ 19 * (dim + x) + i] = 
		                 collideField[19 * (LATTICEVELOCITIES[i][2] * dim_pow2 + 
		           	                  (1 + LATTICEVELOCITIES[i][1]) * dim + x + LATTICEVELOCITIES[i][0]) + 18 - i];
	    }		    
 	// y = 0, i = 18
 	for( x = 1; x < xlength + 1; x++ ){

 		collideField[ 19 * x + 18] = collideField[19 * (LATTICEVELOCITIES[18][2] * dim_pow2 + 
 			           	                  LATTICEVELOCITIES[18][1] * dim + x + LATTICEVELOCITIES[18][0])];
 	}

 	// z = 0, y = xlength & y = xlength + 1
 	// y = xlength
 	for( x = 2; x < xlength; x++ )
			for( i = 14; i < 18; i++ ){

			     collideField[ 19 * (xlength * dim + x) + i] = 
			                 collideField[19 * (LATTICEVELOCITIES[i][2] * dim_pow2 + 
			           	                  (y + LATTICEVELOCITIES[i][1]) * dim + x + LATTICEVELOCITIES[i][0]) + 18 - i];
		    }
 	// y = xlength + 1 & i = 14
    for( x = 1; x < xlength + 1; x++ ){
    	 	collideField[ 19 * (xlength + 1 * dim + x) + 14] = 
 			                 collideField[19 * (LATTICEVELOCITIES[14][2] * dim_pow2 + 
 			           	                  (xlength + 1 + LATTICEVELOCITIES[14][1]) * dim + x + LATTICEVELOCITIES[14][0]) + 4];
    }

    // x = 1 & z = 0 & y = 2:xlength - 1
 	for( y = 2; y < xlength; y++ ){
 			// i = 14
 		    collideField[ 19 * (y * dim + 1) + 14] = 
 			                 collideField[19 * (LATTICEVELOCITIES[14][2] * dim_pow2 + 
 			           	                  (y + LATTICEVELOCITIES[14][1]) * dim + 1 + LATTICEVELOCITIES[14][0]) + 14];
 			for( i = 16; i < 19; i++ ){

 			     collideField[ 19 * (y * dim + 1) + i] = 
 			                 collideField[19 * (LATTICEVELOCITIES[i][2] * dim_pow2 + 
 			           	                  (y + LATTICEVELOCITIES[i][1]) * dim + 1 + LATTICEVELOCITIES[i][0]) + 18 - i];
 		    }
 	}
 	// x = 0 & z = 0 & y = 1:xlength & i = 17
 	for( y = 1; y < xlength + 1; y++ ){

 			     collideField[ 19 * (y * dim) + 17] = 
 			                 collideField[19 * (LATTICEVELOCITIES[17][2] * dim_pow2 + 
 			           	                  (y + LATTICEVELOCITIES[17][1]) * dim + LATTICEVELOCITIES[17][0]) + 1];
 	}

 	// x = xlength || x = xlength + 1 & z = 0
 	// x = xlength
 	for( y = 2; y < xlength; y++ ){
 		    // i =14:16
 			for( i = 14; i < 17; i++ )
 			     collideField[ 19 * (y * dim + xlength) + i] = 
 			                 collideField[19 * (LATTICEVELOCITIES[i][2] * dim_pow2 + 
 			           	                  (y + LATTICEVELOCITIES[i][1]) * dim + xlength + LATTICEVELOCITIES[i][0]) + 18 - i];
 			// i = 18
 			collideField[ 19 * (y * dim + xlength) + 18] = 
 			                 collideField[19 * (LATTICEVELOCITIES[18][2] * dim_pow2 + 
 			           	                  (y + LATTICEVELOCITIES[18][1]) * dim + xlength + LATTICEVELOCITIES[18][0])];
 	}
 	// x = xlength + 1, i = 15
 	for( y = 1; y < xlength + 1; y++ ){

 			     collideField[ 19 * (y * dim + xlength + 1) + 15] = 
 			                 collideField[19 * (LATTICEVELOCITIES[15][2] * dim_pow2 + 
 			           	                  (y + LATTICEVELOCITIES[15][1]) * dim + xlength + 1 + LATTICEVELOCITIES[15][0]) + 3];
 		    } 	
    // four corner points
 	// x = 1, y = 1, z = 0, i = 16,17,18
 	for( i = 16; i < 19; i++ ){

 			     collideField[ 19 * (dim + 1) + i] = 
 			                 collideField[19 * (LATTICEVELOCITIES[i][2] * dim_pow2 + 
 			           	                  (1 + LATTICEVELOCITIES[i][1]) * dim + 1 + LATTICEVELOCITIES[i][0]) + 18 - i];
 		    }
    // x = xlength, y = 1, z = 0, i = 15,16,18
 	for( i = 15; i < 19; i++ ){
 			if(i != 17 )
 			     collideField[ 19 * (dim + xlength) + i] = 
 			                 collideField[19 * (LATTICEVELOCITIES[i][2] * dim_pow2 + 
 			           	                  (1 + LATTICEVELOCITIES[i][1]) * dim + xlength + LATTICEVELOCITIES[i][0]) + 18 - i];
 	}
 	// x = 1, y = xlength, z = 0, i = 14,16,17
 	for( i = 14; i < 18; i++ ){
 			if(i != 15)
 			     collideField[ 19 * (xlength * dim + 1) + i] = 
 			                 collideField[19 * (LATTICEVELOCITIES[i][2] * dim_pow2 + 
 			           	                  (xlength + LATTICEVELOCITIES[i][1]) * dim + 1 + LATTICEVELOCITIES[i][0]) + 18 - i];
    }
    // x = xlength, y = xlength, z = 0, i = 14,15,16
	for( i = 14; i < 17; i++ ){

	     collideField[ 19 * (xlength * dim + xlength) + i] = 
	                 collideField[19 * (LATTICEVELOCITIES[i][2] * dim_pow2 + 
	           	                  (xlength + LATTICEVELOCITIES[i][1]) * dim + xlength + LATTICEVELOCITIES[i][0]) + 18 - i];
    }    


    /********************/
  	/** side face,x = 0 */
  	/********************/
  	// interior, i = 3, 7, 10, 13, 17 
  	Iindex[0] = 3; Iindex[1] = 7; Iindex[2] = 10; Iindex[3] = 13; Iindex[4] = 17;
    for( y = 2; y < xlength; y++ )
 		for( z = 2; z < xlength; z++ ){
 			for( i = 0; i < 5; i++ )
 				collideField[ 19 * (z * dim_pow2 + y * dim) + Iindex[i] ] = collideField[ 19 * ( (z + LATTICEVELOCITIES[Iindex[i]][2]) * dim_pow2 + 
 					                                          (y + LATTICEVELOCITIES[Iindex[i]][1]) * dim + LATTICEVELOCITIES[Iindex[i]][0] ) + 18 - Iindex[i]];
 	    }

 	// y = 1, z = 2: xlength - 1, i = 3,10,13,17
	for( z = 2; z < xlength; z++ ){
		for( i = 0; i < 5; i++ )
			if(i != 1)
			collideField[ 19 * (z * dim_pow2 + dim) + Iindex[i] ] = collideField[ 19 * ( (z + LATTICEVELOCITIES[Iindex[i]][2]) * dim_pow2 + 
				                                          (1 + LATTICEVELOCITIES[Iindex[i]][1]) * dim + LATTICEVELOCITIES[Iindex[i]][0] ) + 18 - Iindex[i]];
    }    

    // y = xlength, z = 2: xlength - 1, i = 3,7,10,17
	for( z = 2; z < xlength; z++ ){
		for( i = 0; i < 5; i++ )
			if(1 != 3)
			collideField[ 19 * (z * dim_pow2 + xlength * dim) + Iindex[i] ] = collideField[ 19 * ( (z + LATTICEVELOCITIES[Iindex[i]][2]) * dim_pow2 + 
				                                          (xlength + LATTICEVELOCITIES[Iindex[i]][1]) * dim + LATTICEVELOCITIES[Iindex[i]][0] ) + 18 - Iindex[i]];
    }

    // z = 1, y = 2:xlength-1, i = 7,10,13,17
    for( y = 2; y < xlength; y++ ){
 			for( i = 1; i < 5; i++ )
 				collideField[ 19 * (dim_pow2 + y * dim) + Iindex[i] ] = collideField[ 19 * ( (1 + LATTICEVELOCITIES[Iindex[i]][2]) * dim_pow2 + 
 					                                          (y + LATTICEVELOCITIES[Iindex[i]][1]) * dim + LATTICEVELOCITIES[Iindex[i]][0] ) + 18 - Iindex[i]];
 	    }
 	// z = xlength, y = 2:xlength-1, i = 3, 7, 10, 13
    for( y = 2; y < xlength; y++ ){
 			for( i = 0; i < 4; i++ )
 				collideField[ 19 * (xlength * dim_pow2 + y * dim) + Iindex[i] ] = collideField[ 19 * ( (xlength + LATTICEVELOCITIES[Iindex[i]][2]) * dim_pow2 + 
 					                                          (y + LATTICEVELOCITIES[Iindex[i]][1]) * dim + LATTICEVELOCITIES[Iindex[i]][0] ) + 18 - Iindex[i]];
 	    }
 	// four corner cells at x = 0
 	// y = 1, z = 1, x = 0, i = 10,13,17{
	for( i = 2; i < 5; i++ ){
		collideField[ 19 * (dim_pow2 + dim) + Iindex[i] ] = collideField[ 19 * ( (1 + LATTICEVELOCITIES[Iindex[i]][2]) * dim_pow2 + 
			                                          (1 + LATTICEVELOCITIES[Iindex[i]][1]) * dim + LATTICEVELOCITIES[Iindex[i]][0] ) + 18 - Iindex[i]];
    }
    // y = xlength, z = 1, x = 0, i = 7,10,17
	for( i = 0; i < 5; i++ )
		if( i == 1 || i == 2 || i == 4){
		    collideField[ 19 * (dim_pow2 + xlength * dim) + Iindex[i] ] = collideField[ 19 * ( (1 + LATTICEVELOCITIES[Iindex[i]][2]) * dim_pow2 + 
			                                          (xlength + LATTICEVELOCITIES[Iindex[i]][1]) * dim + LATTICEVELOCITIES[Iindex[i]][0] ) + 18 - Iindex[i]];
	    }
	// z= xlength, y = 1, x = 0; i = 3,10,13
	for( i = 0; i < 4; i++ )
		if(i != 1){
		collideField[ 19 * (xlength * dim_pow2 + dim) + Iindex[i] ] = collideField[ 19 * ( (xlength + LATTICEVELOCITIES[Iindex[i]][2]) * dim_pow2 + 
			                                          (1 + LATTICEVELOCITIES[Iindex[i]][1]) * dim + LATTICEVELOCITIES[Iindex[i]][0] ) + 18 - Iindex[i]];
	}
	// z = xlength, y = xlength; i = 3,7,10,
	for( i = 0; i < 3; i++ ){
		collideField[ 19 * (xlength * dim_pow2 + xlength * dim) + Iindex[i] ] = collideField[ 19 * ( (xlength + LATTICEVELOCITIES[Iindex[i]][2]) * dim_pow2 + 
			                                          (xlength + LATTICEVELOCITIES[Iindex[i]][1]) * dim + LATTICEVELOCITIES[Iindex[i]][0] ) + 18 - Iindex[i]];
	}

 	/********************/
  	/** side face,y = 0 */
  	/********************/
  	// interior, i = 4,11,12,13,18
	Iindex[0] = 4; Iindex[1] = 11; Iindex[2] = 12; Iindex[3] = 13; Iindex[4] = 18;
 	for( x = 2; x < xlength; x++ )
 		for( z = 2; z < xlength; z++ )
 			for( i = 0; i < 5; i++ ){
 				collideField[ 19 * (z * dim_pow2 + x) + Iindex[i] ] = collideField[ 19 * ((z + LATTICEVELOCITIES[Iindex[i]][2]) * dim_pow2 + 
 																	LATTICEVELOCITIES[Iindex[i]][1] * dim + x + LATTICEVELOCITIES[Iindex[i]][0]) + 18 - Iindex[i]];
 			}
 	// z = 1, x = 2:xlength-1, i = 11,12,13,18
 	for( x = 2; x < xlength; x++ )
		for( i = 1; i < 5; i++ ){
			collideField[ 19 * (dim_pow2 + x) + Iindex[i] ] = collideField[ 19 * ((1 + LATTICEVELOCITIES[Iindex[i]][2]) * dim_pow2 + 
																LATTICEVELOCITIES[Iindex[i]][1] * dim + x + LATTICEVELOCITIES[Iindex[i]][0]) + 18 - Iindex[i]];
		}
	// z = xlength,x = 2:xlength-1, i = 4,11,12,13
 	for( x = 2; x < xlength; x++ )
 			for( i = 0; i < 4; i++ ){
 				collideField[ 19 * (xlength * dim_pow2 + x) + Iindex[i] ] = collideField[ 19 * ((xlength + LATTICEVELOCITIES[Iindex[i]][2]) * dim_pow2 + 
 																	LATTICEVELOCITIES[Iindex[i]][1] * dim + x + LATTICEVELOCITIES[Iindex[i]][0]) + 18 - Iindex[i]];
 			}
 	// x = 1, z = 2:xlength-1, i = 4,12,13,18
	for( z = 2; z < xlength; z++ )
		for( i = 0; i < 5; i++ )
			if(i != 1){
			collideField[ 19 * (z * dim_pow2 + 1) + Iindex[i] ] = collideField[ 19 * ((z + LATTICEVELOCITIES[Iindex[i]][2]) * dim_pow2 + 
																LATTICEVELOCITIES[Iindex[i]][1] * dim + 1 + LATTICEVELOCITIES[Iindex[i]][0]) + 18 - Iindex[i]];
			}
	// x = xlength,z = 2:xlength-1, i = 4,11,12,18
	for( z = 2; z < xlength; z++ )
		for( i = 0; i < 5; i++ )
			if(i != 3){
			   collideField[ 19 * (z * dim_pow2 + xlength) + Iindex[i] ] = collideField[ 19 * ((z + LATTICEVELOCITIES[Iindex[i]][2]) * dim_pow2 + 
																LATTICEVELOCITIES[Iindex[i]][1] * dim + xlength + LATTICEVELOCITIES[Iindex[i]][0]) + 18 - Iindex[i]];
			}
	// four corner cells
	// x=1,z=1,i = 11,12,18
	for( i = 1; i < 5; i++ )
		if(i != 3){
		collideField[ 19 * (dim_pow2 + 1) + Iindex[i] ] = collideField[ 19 * ((1 + LATTICEVELOCITIES[Iindex[i]][2]) * dim_pow2 + 
															LATTICEVELOCITIES[Iindex[i]][1] * dim + 1 + LATTICEVELOCITIES[Iindex[i]][0]) + 18 - Iindex[i]];
	}
	// x=1, z= xlength, i = 4,12,13
	for( i = 0; i < 4; i++ )
		if(i != 1){
		    collideField[ 19 * (xlength * dim_pow2 + 1) + Iindex[i] ] = collideField[ 19 * ((xlength + LATTICEVELOCITIES[Iindex[i]][2]) * dim_pow2 + 
															LATTICEVELOCITIES[Iindex[i]][1] * dim + 1 + LATTICEVELOCITIES[Iindex[i]][0]) + 18 - Iindex[i]];
		}
	// x = xlength, z = 1, i = 11,12,18
	for( i = 1; i < 5; i++ )
		if(i != 3){
		collideField[ 19 * (dim_pow2 + xlength) + Iindex[i] ] = collideField[ 19 * ((1 + LATTICEVELOCITIES[Iindex[i]][2]) * dim_pow2 + 
															LATTICEVELOCITIES[Iindex[i]][1] * dim + xlength + LATTICEVELOCITIES[Iindex[i]][0]) + 18 - Iindex[i]];
	    }
	// x=xlength, z=xlength, i = 4,11,12
	for( i = 0; i < 3; i++ ){
		collideField[ 19 * (xlength * dim_pow2 + xlength) + Iindex[i] ] = collideField[ 19 * ((xlength + LATTICEVELOCITIES[Iindex[i]][2]) * dim_pow2 + 
															LATTICEVELOCITIES[Iindex[i]][1] * dim + xlength + LATTICEVELOCITIES[Iindex[i]][0]) + 18 - Iindex[i]];
	}

 	/******************************/
  	/** side face,x = xlength + 1 */
  	/******************************/
  	//interior, i = 1,5,8,11,15
    Iindex[0] = 1; Iindex[1] = 5; Iindex[2] = 8; Iindex[3] = 11; Iindex[4] = 15;
 	for( y = 2; y < xlength; y++ )
 		for( z = 2; z < xlength; z++ )
 			for( i = 0; i < 5; i++ ){
				collideField[ 19 * (z * dim_pow2 + y * dim + dim - 1) + Iindex[i] ] = collideField[ 19 * ( (z + LATTICEVELOCITIES[Iindex[i]][2]) * dim_pow2 +
																(y + LATTICEVELOCITIES[Iindex[i]][1]) * dim + dim - 1 + LATTICEVELOCITIES[Iindex[i]][0]) + 18 - Iindex[i]]; 				
 			}

 	// z=1,y=2:xlength-1, i = 5,8,11,15
  	for( y = 2; y < xlength; y++ )
 			for( i = 1; i < 5; i++ ){
				collideField[ 19 * (dim_pow2 + y * dim + dim - 1) + Iindex[i] ] = collideField[ 19 * ( (1 + LATTICEVELOCITIES[Iindex[i]][2]) * dim_pow2 +
																(y + LATTICEVELOCITIES[Iindex[i]][1]) * dim + dim - 1 + LATTICEVELOCITIES[Iindex[i]][0]) + 18 - Iindex[i]]; 				
 			}	
 	// z=xlength,y=2:xlength-1, i =1,5,8,11
 	for( y = 2; y < xlength; y++ )
		for( i = 0; i < 4; i++ ){
		    collideField[ 19 * (xlength * dim_pow2 + y * dim + dim - 1) + Iindex[i] ] = collideField[ 19 * ( (xlength + LATTICEVELOCITIES[Iindex[i]][2]) * dim_pow2 +
														(y + LATTICEVELOCITIES[Iindex[i]][1]) * dim + dim - 1 + LATTICEVELOCITIES[Iindex[i]][0]) + 18 - Iindex[i]]; 				
		}
	// y=1, y=2:xlength-1, i = 1,8,11,15
	for( z = 2; z < xlength; z++ )
		for( i = 0; i < 5; i++ )
			if(i != 1){
	    	   collideField[ 19 * (z * dim_pow2 + dim + dim - 1) + Iindex[i] ] = collideField[ 19 * ( (z + LATTICEVELOCITIES[Iindex[i]][2]) * dim_pow2 +
														(1 + LATTICEVELOCITIES[Iindex[i]][1]) * dim + dim - 1 + LATTICEVELOCITIES[Iindex[i]][0]) + 18 - Iindex[i]]; 				
		    }
    // y=xlength, z=2:xlength-1, i =1,5,8,15
	for( z = 2; z < xlength; z++ )
		for( i = 0; i < 5; i++ )
			if(i != 3){
		    collideField[ 19 * (z * dim_pow2 + xlength * dim + dim - 1) + Iindex[i] ] = collideField[ 19 * ( (z + LATTICEVELOCITIES[Iindex[i]][2]) * dim_pow2 +
														(xlength + LATTICEVELOCITIES[Iindex[i]][1]) * dim + dim - 1 + LATTICEVELOCITIES[Iindex[i]][0]) + 18 - Iindex[i]]; 				
		    }
	// four corner cells
    // y=1,z=1, i = 8,11,15
	for( i = 2; i < 5; i++ ){
	    collideField[ 19 * (dim_pow2 + dim + dim - 1) + Iindex[i] ] = collideField[ 19 * ( (1 + LATTICEVELOCITIES[Iindex[i]][2]) * dim_pow2 +
													(1 + LATTICEVELOCITIES[Iindex[i]][1]) * dim + dim - 1 + LATTICEVELOCITIES[Iindex[i]][0]) + 18 - Iindex[i]]; 				
	}   
	// y=1, z=xlength, i= 1,8,11
	for( i = 0; i < 4; i++ )
		if(i != 1){
	    collideField[ 19 * (xlength * dim_pow2 + dim + dim - 1) + Iindex[i] ] = collideField[ 19 * ( (xlength + LATTICEVELOCITIES[Iindex[i]][2]) * dim_pow2 +
													(1 + LATTICEVELOCITIES[Iindex[i]][1]) * dim + dim - 1 + LATTICEVELOCITIES[Iindex[i]][0]) + 18 - Iindex[i]]; 				
	}
	// y=xlength,z=1, i= 5,8,15
	for( i = 1; i < 5; i++ )
		if(i != 2){
	      collideField[ 19 * (dim_pow2 + xlength * dim + dim - 1) + Iindex[i] ] = collideField[ 19 * ( (1 + LATTICEVELOCITIES[Iindex[i]][2]) * dim_pow2 +
													(xlength + LATTICEVELOCITIES[Iindex[i]][1]) * dim + dim - 1 + LATTICEVELOCITIES[Iindex[i]][0]) + 18 - Iindex[i]]; 				
	    }

	// y=xlength,z=xlength, i=1,5,8
	for( i = 0; i < 3; i++ ){
	    collideField[ 19 * (xlength * dim_pow2 + xlength * dim + dim - 1) + Iindex[i] ] = collideField[ 19 * ( (xlength + LATTICEVELOCITIES[Iindex[i]][2]) * dim_pow2 +
													(xlength + LATTICEVELOCITIES[Iindex[i]][1]) * dim + dim - 1 + LATTICEVELOCITIES[Iindex[i]][0]) + 18 - Iindex[i]]; 				
	}

 	/**************************/
  	/** side face,y = dim - 1 */
  	/**************************/
  	//interior, i = 0,5,6,7,14
    Iindex[0] = 0; Iindex[1] = 5; Iindex[2] = 6; Iindex[3] = 7; Iindex[4] = 14;

 	for( x = 2; x < xlength; x++ )
 		for( z = 2; z < xlength; z++ )
 			for( i = 0; i < 5; i++ ){
 				collideField[ 19 * (z * dim_pow2 + xlen_pow2 + x) + Iindex[i] ] = collideField[ 19 * ( (z + LATTICEVELOCITIES[Iindex[i]][2]) * dim_pow2 +
 																	(dim - 1 + LATTICEVELOCITIES[Iindex[i]][1]) * dim + x + LATTICEVELOCITIES[Iindex[i]][0]) + 18 - Iindex[i] ];

 			}

 	// x=1, z=2:xlength-1, i=0,6,7,14
	for( z = 2; z < xlength; z++ )
		for( i = 0; i < 5; i++ )
			if(i != 1){
			collideField[ 19 * (z * dim_pow2 + xlen_pow2 + 1) + Iindex[i] ] = collideField[ 19 * ( (z + LATTICEVELOCITIES[Iindex[i]][2]) * dim_pow2 +
																(dim - 1 + LATTICEVELOCITIES[Iindex[i]][1]) * dim + 1 + LATTICEVELOCITIES[Iindex[i]][0]) + 18 - Iindex[i] ];

		    }
    // x=xlength,z=2:xlength-1, i=0,5,6,14
	for( z = 2; z < xlength; z++ )
		for( i = 0; i < 5; i++ )
			if(i != 3){
			  collideField[ 19 * (z * dim_pow2 + xlen_pow2 + xlength) + Iindex[i] ] = collideField[ 19 * ( (z + LATTICEVELOCITIES[Iindex[i]][2]) * dim_pow2 +
																(dim - 1 + LATTICEVELOCITIES[Iindex[i]][1]) * dim + xlength + LATTICEVELOCITIES[Iindex[i]][0]) + 18 - Iindex[i] ];

		    }
	// z=1, x=2:xlength-1, i=5,6,7,14
 	for( x = 2; x < xlength; x++ )
 			for( i = 1; i < 5; i++ ){
 				collideField[ 19 * (dim_pow2 + xlen_pow2 + x) + Iindex[i] ] = collideField[ 19 * ( (1 + LATTICEVELOCITIES[Iindex[i]][2]) * dim_pow2 +
 																	(dim - 1 + LATTICEVELOCITIES[Iindex[i]][1]) * dim + x + LATTICEVELOCITIES[Iindex[i]][0]) + 18 - Iindex[i] ];

 			}
 	// z=xlength,x=2:xlength-1, i=0,5,6,7
  	for( x = 2; x < xlength; x++ )
 			for( i = 0; i < 4; i++ ){
 				collideField[ 19 * (xlength * dim_pow2 + xlen_pow2 + x) + Iindex[i] ] = collideField[ 19 * ( (xlength + LATTICEVELOCITIES[Iindex[i]][2]) * dim_pow2 +
 																	(dim - 1 + LATTICEVELOCITIES[Iindex[i]][1]) * dim + x + LATTICEVELOCITIES[Iindex[i]][0]) + 18 - Iindex[i] ];

 			}
 	//four corner cells:
 	// x=1, z=1, i= 6,7,14
	for( i = 2; i < 5; i++ ){
		collideField[ 19 * (dim_pow2 + xlen_pow2 + 1) + Iindex[i] ] = collideField[ 19 * ( (1 + LATTICEVELOCITIES[Iindex[i]][2]) * dim_pow2 +
															(dim - 1 + LATTICEVELOCITIES[Iindex[i]][1]) * dim + 1 + LATTICEVELOCITIES[Iindex[i]][0]) + 18 - Iindex[i] ];

	}
	// x=xlength,z=1,i=5,6,14
	for( i = 1; i < 5; i++ )
		if(i != 3){
		collideField[ 19 * (dim_pow2 + xlen_pow2 + xlength) + Iindex[i] ] = collideField[ 19 * ( (1 + LATTICEVELOCITIES[Iindex[i]][2]) * dim_pow2 +
															(dim - 1 + LATTICEVELOCITIES[Iindex[i]][1]) * dim + xlength + LATTICEVELOCITIES[Iindex[i]][0]) + 18 - Iindex[i] ];

		}
	// x=xlength,z=xlength, i=0,5,6,
	for( i = 0; i < 4; i++ ){
		collideField[ 19 * (xlength * dim_pow2 + xlen_pow2 + xlength) + Iindex[i] ] = collideField[ 19 * ( (xlength + LATTICEVELOCITIES[Iindex[i]][2]) * dim_pow2 +
															(dim - 1 + LATTICEVELOCITIES[Iindex[i]][1]) * dim + xlength + LATTICEVELOCITIES[Iindex[i]][0]) + 18 - Iindex[i] ];

	}
	// x=1, z=xlength, i=0,6,7,
	for( i = 0; i < 4; i++ )
		if(i != 1){
		collideField[ 19 * (xlength * dim_pow2 + xlen_pow2 + 1) + Iindex[i] ] = collideField[ 19 * ( (xlength + LATTICEVELOCITIES[Iindex[i]][2]) * dim_pow2 +
															(dim - 1 + LATTICEVELOCITIES[Iindex[i]][1]) * dim + 1 + LATTICEVELOCITIES[Iindex[i]][0]) + 18 - Iindex[i] ];

	}

 	/**************************/
  	/** four side boudary bars*/
  	/**************************/
 	// x = y = 0, z = 1 : xlength, i = 13
	for( z = 1; z < xlength + 1; z++ ){
		    collideField[ 19 * (z * dim_pow2) + 13 ] = collideField[ 19 * ( (z + LATTICEVELOCITIES[13][2]) * dim_pow2 +
															LATTICEVELOCITIES[13][1] * dim + LATTICEVELOCITIES[13][0]) + 5 ]; 				
		}

	// x = xlength + 1, y = 0, z = 1 : xlength, i = 11
	for( z = 1; z < dim - 1; z++ ){

			collideField[ 19 * (z * dim_pow2 + xlength + 1) + 11 ] = collideField[ 19 * ((z + LATTICEVELOCITIES[11][2]) * dim_pow2 + 
																LATTICEVELOCITIES[11][1] * dim + xlength + 1 + LATTICEVELOCITIES[11][0]) + 7];
		}

 	// x = 0, y = xlength + 1, z = 1 : xlength, i = 7
	for( z = 1; z < dim - 1; z++ ){
			collideField[ 19 * (z * dim_pow2 + xlen_pow2) + 7 ] = collideField[ 19 * ( (z + LATTICEVELOCITIES[7][2]) * dim_pow2 +
																(dim - 1 + LATTICEVELOCITIES[7][1]) * dim + LATTICEVELOCITIES[7][0]) + 11];

		}
	// x = xlength + 1, y = xlength + 1, z = 1 : xlength, i = 5
	for( z = 1; z < dim - 1; z++ ){
			collideField[ 19 * (z * dim_pow2 + xlen_pow2 + xlength + 1) + 5 ] = collideField[ 19 * ( (z + LATTICEVELOCITIES[5][2]) * dim_pow2 +
															(dim - 1 + LATTICEVELOCITIES[5][1]) * dim + xlength + 1 + LATTICEVELOCITIES[5][0]) + 13 ];

		}	

    
 	/************************************/
 	/** top z = dim - 1, moving wall!!! */
 	/************************************/
	// compute c_i * u_wall / cs^2
	double val[19], density;
	for(i = 0; i < 19; i++){
		val[i] = wallVelocity[0] * LATTICEVELOCITIES[i][0] + wallVelocity[1] * LATTICEVELOCITIES[i][1] + wallVelocity[2] * LATTICEVELOCITIES[i][2];
		val[i] = val[i] / cs_2;
	}		

	// interior boundary cells, i = 0,1,2,3,4
 	for( x = 2; x < xlength; x++ )
 		for( y = 2; y < xlength; y++ )
 			for( i = 0; i < 5; i++ ){

 				index = (dim - 1 + LATTICEVELOCITIES[i][2]) * dim_pow2 + (y + LATTICEVELOCITIES[i][1]) * dim + x + LATTICEVELOCITIES[i][0];

 				computeDensity(&(collideField[19 * index]), &density);

 				collideField[ 19 * (xlen_pow3 + y * dim + x) + i] = collideField[19 * index + 18 - i] + 2 * LATTICEWEIGHTS[i] * density * val[i];
 			}

 	// z = dim - 1, y = 1, i = 1,2,3,4
 	for( x = 2; x < xlength; x++ )
 			for( i = 1; i < 5; i++ ){

 				index = (dim - 1 + LATTICEVELOCITIES[i][2]) * dim_pow2 + (1 + LATTICEVELOCITIES[i][1]) * dim + x + LATTICEVELOCITIES[i][0];

 				computeDensity(&(collideField[19*index]), &density);

 				collideField[ 19 * (xlen_pow3 + dim + x) + i] = collideField[19 * index + 18 - i] + 2 * LATTICEWEIGHTS[i] * density * val[i];
 			}

 	// y = 0, z = dim - 1, i = 4
 	for( x = 1; x < xlength + 1; x++ ){

 				index = (dim - 1 + LATTICEVELOCITIES[4][2]) * dim_pow2 + LATTICEVELOCITIES[4][1] * dim + x + LATTICEVELOCITIES[4][0];

 				computeDensity(&(collideField[19*index]), &density);

 				collideField[ 19 * (xlen_pow3 + x) + 4] = collideField[19 * index + 14] + 2 * LATTICEWEIGHTS[4] * density * val[4];
 			}

 	// y = xlength, z = dim - 1, i = 0,1,2,3
 	for( x = 2; x < xlength; x++ )
 			for( i = 0; i < 4; i++ ){

 				index = (dim - 1 + LATTICEVELOCITIES[i][2]) * dim_pow2 + (xlength + LATTICEVELOCITIES[i][1]) * dim + x + LATTICEVELOCITIES[i][0];

 				computeDensity(&(collideField[19*index]), &density);

 				collideField[ 19 * (xlen_pow3 + xlength * dim + x) + i] = collideField[19 * index + 18 - i] + 2 * LATTICEWEIGHTS[i] * density * val[i];
 			}

 	// y = xlength+1, z = dim - 1, i = 0
 	for( x = 1; x < xlength + 1; x++ ){

 				index = (dim - 1 + LATTICEVELOCITIES[0][2]) * dim_pow2 + (xlength + 1 + LATTICEVELOCITIES[0][1]) * dim + x + LATTICEVELOCITIES[0][0];

 				computeDensity(&(collideField[19*index]), &density);

 				collideField[ 19 * (xlen_pow3 + (xlength + 1) * dim + x)] = collideField[19 * index + 18] + 2 * LATTICEWEIGHTS[0] * density * val[0];
 			}

 	// x = 1, z = dim - 1, i = 0,2,3,4
	for( y = 2; y < xlength; y++ )
		for( i = 0; i < 5; i++ )
			if(i != 1){

			index = (dim - 1 + LATTICEVELOCITIES[i][2]) * dim_pow2 + (y + LATTICEVELOCITIES[i][1]) * dim + 1 + LATTICEVELOCITIES[i][0];

			computeDensity(&(collideField[19*index]), &density);

			collideField[ 19 * (xlen_pow3 + y * dim + 1) + i] = collideField[19 * index + 18 - i] + 2 * LATTICEWEIGHTS[i] * density * val[i];
		    }

    // x = 0, z = dim - 1, i = 3
	for( y = 1; y < xlength + 1; y++ ){

			index = (dim - 1 + LATTICEVELOCITIES[3][2]) * dim_pow2 + (y + LATTICEVELOCITIES[3][1]) * dim + LATTICEVELOCITIES[3][0];

			computeDensity(&(collideField[19*index]), &density);

			collideField[ 19 * (xlen_pow3 + y * dim) + 3] = collideField[19 * index + 15] + 2 * LATTICEWEIGHTS[3] * density * val[3];
			}

    // x = xlength, z = dim - 1, i = 0,1,2,4
	for( y = 2; y < xlength; y++ )
		for( i = 0; i < 5; i++ )
			if( i != 3){

			index = (dim - 1 + LATTICEVELOCITIES[i][2]) * dim_pow2 + (y + LATTICEVELOCITIES[i][1]) * dim + xlength + LATTICEVELOCITIES[i][0];

			computeDensity(&(collideField[19*index]), &density);

			collideField[ 19 * (xlen_pow3 + y * dim + xlength) + i] = collideField[19 * index + 18 - i] + 2 * LATTICEWEIGHTS[i] * density * val[i];
			}

	// x = xlength + 1, z = dim - 1, i = 1
	for( y = 1; y < xlength + 1; y++ ){

			index = (dim - 1 + LATTICEVELOCITIES[1][2]) * dim_pow2 + (y + LATTICEVELOCITIES[1][1]) * dim + xlength + 1 + LATTICEVELOCITIES[1][0];

			computeDensity(&(collideField[19*index]), &density);

			collideField[ 19 * (xlen_pow3 + y * dim + xlength + 1) + 1] = collideField[19 * index + 17] + 2 * LATTICEWEIGHTS[1] * density * val[1];
		}

	// four corner cells, z = dim - 1
	// y = 1, x = 1, z = dim - 1, i = 2,3,4
	for( i = 2; i < 5; i++ ){

		index = (dim - 1 + LATTICEVELOCITIES[i][2]) * dim_pow2 + (1 + LATTICEVELOCITIES[i][1]) * dim + 1 + LATTICEVELOCITIES[i][0];

		computeDensity(&(collideField[19*index]), &density);

		collideField[ 19 * (xlen_pow3 + dim + 1) + i] = collideField[19 * index + 18 - i] + 2 * LATTICEWEIGHTS[i] * density * val[i];
	}	
	// y = xlength, x = 1, z = dim - 1, i = 0,2,3
	for( i = 0; i < 4; i++ )
		if(i != 1){

		index = (dim - 1 + LATTICEVELOCITIES[i][2]) * dim_pow2 + (xlength + LATTICEVELOCITIES[i][1]) * dim + 1 + LATTICEVELOCITIES[i][0];

		computeDensity(&(collideField[19*index]), &density);

		collideField[ 19 * (xlen_pow3 + xlength * dim + 1) + i] = collideField[19 * index + 18 - i] + 2 * LATTICEWEIGHTS[i] * density * val[i];
		}
	// y = 1, x = xlength, z = dim - 1, i = 1,2,4
	for( i = 1; i < 5; i++ )
		if(i != 3){

		index = (dim - 1 + LATTICEVELOCITIES[i][2]) * dim_pow2 + (1 + LATTICEVELOCITIES[i][1]) * dim + xlength + LATTICEVELOCITIES[i][0];

		computeDensity(&(collideField[19*index]), &density);

		collideField[ 19 * (xlen_pow3 + dim + xlength) + i] = collideField[19 * index + 18 - i] + 2 * LATTICEWEIGHTS[i] * density * val[i];
		}
	// y = xlength, x = xlength, z = dim - 1, i = 0,1, 2
	for( i = 0; i < 3; i++ ){

		index = (dim - 1 + LATTICEVELOCITIES[i][2]) * dim_pow2 + (xlength + LATTICEVELOCITIES[i][1]) * dim + xlength + LATTICEVELOCITIES[i][0];

		computeDensity(&(collideField[19*index]), &density);

		collideField[ 19 * (xlen_pow3 + xlength * dim + xlength) + i] = collideField[19 * index + 18 - i] + 2 * LATTICEWEIGHTS[i] * density * val[i];
	}



}

