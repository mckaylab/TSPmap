#include <stdio.h>
#include <stdlib.h>
#include <math.h>   // for round function.

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>	// Need this so the compiler recognizes the AS_NUMERIC function.

#include <Rcpp.h>
using namespace Rcpp;

/*
 * rf.c
 *
 *  Created on: Jan 12, 2015
 *      Author: Zachariah A. Allen
 */

/*
FOR USE WITH R:

Create a shared object file from this file using:
R CMD SHLIB rf.c

Then load the object file into R using:
dyn.load("/Users/zallen/Documents/msProject/Rpackage/TSPmap/src/rf.so")
*/


// Find duplicates in the raw data.
// Exact duplicates are markers which are the same in every cell, including missing cells (indicated by 0 values).
// Similar markers are measured by ignoring any spot where either marker is missing data.
// This function takes a similarity threshold as a parameter, which is used to decide how similar two markers have to be in order to be considered similar enough that one will be removed.  The default value of threshold is -1.  The default value means that the user wishes to only identify exact duplicate markers, meaning that both markers must be missing data in all the same spots.  Otherwise, only the spots where neither marker is missing data is considered.
// NOTE: a marker may be identified as a duplicate more than once, so the return list will have duplicate values.
// Eliminate these in R using the unique() function within the finddups function.

// [[Rcpp::export]]
SEXP findDups(SEXP input, SEXP threshold)
{
	printf("\n\nFinding duplicate markers.\n");

	int i, j, k;
	int nprot = 0;

    float similarity = 0;

    // When comparing two columns of the input matrix:
    // same keeps track of the number of cells that are the same.
    // missing keeps track of the number of cells where at least one of the columns is zero (which corresponds to missing data in the original data matrix).
    // missingA and missingB keep track of how many cells are missing in the two columns that we are comparing.  If they are duplicates, we removed the one with the higher number of missing cells.
	int same, missing, missingA, missingB;

    // Convert threshold parameter to an integer.
    float thresh;

    thresh = *REAL(threshold) * 10;
    //printf("-----> The threshold is %f\n", thresh);

	// Get dimensions of input matrix so we can create the output matrix.
	// dims[0] is the numbers of A/B markers, dims[1] is the total number of markers.
	int *dims;
	dims = INTEGER(GET_DIM(input));

	// Shortcut for dims[0] since we will use it a lot.
	int c = dims[0];
	//printf("c = %i, dims[1] = %i\n", c, dims[1]);

	// create return matrix.
	SEXP ans;
	PROTECT(ans = allocVector(REALSXP, dims[1]));
	nprot++;

	// Set all elements of ans to -1 so we know which elements to ignore when we return to R.
	for(i = 0; i < dims[1]; i++)
		REAL(ans)[i] = -1;

	// Step through all markers.
	for(i = 0; i < dims[1]-1; i++)
	//for(i = 1516; i < 1517; i++)
	{
        //printf("--------> HERE WE ARE #1, i = %i <----------\n", i);

		if(i%100 == 0)
		{
			printf(".");
			if(i % 500 == 0)
				printf(" Processed %i of %i markers\n", i, dims[1]);
		}

		// Step through all markers from i to the end.
		for(j = i+1; j < dims[1]; j++)
		//for(j = 1532; j < 1533; j++)
		{
			same = 0;
            missing = 0;
            missingA = 0;
            missingB = 0;

			//printf("HERE WE ARE #2, i, j = %i, %i\n", i, j);

			// Step through each individual A/B values and compare them.
			// Start at k = 3 to ignore the header lines.
			for(k = 0; k < dims[0]; k++)
			{
				//printf("k = %i, marker1 = %i, marker2 = %i\n", k, INTEGER(input)[i*c + k], INTEGER(input)[j*c + k]);
				// If only one value is missing we don't count that as being the same
				// because two markers could be missing values in different spots and we don't
				// want to assume that they will have the same recom fractions when compared
				// with all other markers.
				// Therefore both elements must be the same to count toward the total.
				if(INTEGER(input)[i*c + k] == INTEGER(input)[j*c + k])
                {
                    same++;


                    //printf("same = %i\n", same);

                    // If both values are zero, increase the missing conter.
                    if(INTEGER(input)[i*c + k] == 0)
                        missing++;
                }
                // Otherwise, if only one of the markers is zero, increase the missing counter for that marker.
                else if(INTEGER(input)[i*c + k] == 0)
                {
                    missingA++;
                }
                else if(INTEGER(input)[j*c + k] == 0)
                {
                    missingB ++;
                }
			}

            // If thresh == -1 (default value), we don't care about similarity, only exact duplicate markers.
            if(thresh < 0)
            {
                //printf("i = %i, j = %i, same = %i, c = %i\n", i, j, same, c);
                // If all values are the same, we have found a duplicate.
                if(same == c)
                {
                    //if(i <= 10)
                        //printf("found a duplicate at i = %i, j = %i\n", i, j);

                    // Store the value i in position j, this way R will know the indices of the original markers (i) and can extract the indices of the duplicates that need to be removed simply by finding all indices that are not -1.
                    // Only insert the value of i if another i has not already been put in that slot.
                    // Actually we use the value i+1 since R starts numbering array elements at 1.
                    if(REAL(ans)[j] == -1)
                        REAL(ans)[j] = i+1;
                }
                //printf("ans = %f\n", REAL(ans)[i*dims[1] + j]);
                //printf("ic + j = %i\n", i*dims[1] + j);
            }
            else
            {
                // If threshold is not the default value, we want to check if the number of cells that are the same (only counting cells where neither marker is missing data) is above the threshold value.

                similarity = round((float)same*1000 / (float)c);

                //printf("missing =  %i, same = %i, c = %i, thresh = %f, test = %f\n", missing, same, c, thresh, (float)same*100.0/(float)(c));
                //printf("similarity =  %f, thresh = %f, test = %f\n", similarity, thresh, (float)same*100.0/(float)(c));

                if(similarity >= thresh)
                {
                    // Store the value i in position j, this way R will know the indices of the original markers (i) and can extract the indices of the duplicates that need to be removed simply by finding all indices that are not -1.
                    // Only insert the value of i if another i has not already been put in that slot.
                    // Actually we use the value i+1 since R starts numbering array elements at 1.

                    // We want to keep the marker that has more data.  If they are missing the same number of cells, we don't care which one we keep.
                    if(missingA <= missingB)
                    {
                        if(REAL(ans)[j] == -1)
                            REAL(ans)[j] = i+1;
                    }
                    else
                    {
                        if(REAL(ans)[i] == -1)
                            REAL(ans)[i] = j+1;
                    }
                }
            }
		}
	}

	UNPROTECT(nprot);
	return ans;
}

// end findDups() =============================================================================
// ============================================================================================


// Here we want to find markers which have 90% similarity in the recombination frequency matrix.
SEXP findSimilars(SEXP input)
{
	printf("\n\nFinding similar markers.\n");

	int i, j, k;
	int nprot = 0;
	int same;

	// Get dimensions of input matrix so we can create the output matrix.
	// dims[0] is the numbers of A/B markers, dims[1] is the total number of markers.
	int *dims;
	dims = INTEGER(GET_DIM(input));

	// Shortcut for dims[0] since we will use it a lot.
	int c = dims[0];
	//printf("c = %i, dims[1] = %i\n", c, dims[1]);

	// create return matrix.
	SEXP ans;
	PROTECT(ans = allocVector(REALSXP, dims[1]));
	nprot++;

	// Set all elements of ans to -1 so we know which elements to ignore when we return to R.
	for(i = 0; i < dims[1]; i++)
		REAL(ans)[i] = -1;

	// Step through all markers.
	// Start at i = 1 to ignore the marker name column.
	for(i = 0; i < dims[1]-1; i++)
	//for(i = 4; i < 5; i++)		// TESTING
	{
		if(i%100 == 0)
		{
			printf(".");
			if(i % 500 == 0)
				printf(" Processed %i markers\n", i);
		}

		// Step through all markers from i to the end.
		for(j = i+1; j < dims[1]; j++)
		//for(j = i+1; j < i+2; j++)		// TESTING
		{
			same = 0;

			//printf("HERE WE ARE #2, i, j = %i, %i\n", i, j);

			// Step through each individual A/B values and compare them.
			// Start at k = 3 to ignore the header lines.
			for(k = 0; k < dims[0]; k++)
			{
				//if(k < 100)
					//printf("k = %i, marker1 = %f, marker2 = %f\n", k, REAL(input)[i*c + k], REAL(input)[j*c + k]);
				// If both values are equal (or very close to equal), increase the same marker by 1.

				//printf("k = %i, first = %f, second = %f\n", k, REAL(input)[i*c + k], REAL(input)[j*c + k]);
				//printf("difference %f\n", fabs(REAL(input)[i*c + k] - REAL(input)[j*c + k]));
				if(fabs(REAL(input)[i*c + k] - REAL(input)[j*c + k]) <= 0.0001)
				{
					same++;
					//if(k < 100)
						//printf("found a match at k = %i, difference is %f\n", k, fabs(REAL(input)[i*c + k] - REAL(input)[j*c + k]));
				}
			}

			//printf("i = %i, j = %i, same = %i, c = %i\n", i, j, same, c);
			// If 90% of the values are the same, we have found a similar marker.
			// Remember that the length of the marker list is 2 less that the dimension of the input matrix.
			if(same >= 0.9*(c-2))
			{
				// We want to remove the marker with the lesser number of valid individuals.
				// This is because the marker with more info is more valuable, so we discard the marker with less info.
				// The number of valid individuals is stored in the diagonal elements of the input matrix.

				if(REAL(input)[i*c + i] < REAL(input)[j*c + j])
				{
					REAL(ans)[i] = i;
				}
				else
				{
					REAL(ans)[j] = j;
				}
				//printf("found a similar at i = %i, j = %i, same = %i\n", i, j, same);

			}
			//printf("ans = %f\n", REAL(ans)[i*dims[1] + j]);
			//printf("ic + j = %i\n", i*dims[1] + j);
		}
	}

	UNPROTECT(nprot);
	return ans;
}

// end findSimilars() =============================================================================
// ============================================================================================



// Compute the entire recombination frequency matrix.
// Values of zero represent missing data.
// NOTE: this returns the FULL matrix, not just the lower triangle.
// correctionFlag - this flag controls whether or not the data correction term is used when computing the recombination frequency.
SEXP computeRF(SEXP input, SEXP correctionFlag)
{
	printf("\n\nComputing recombination frequency matrix.\n");

	int i, j, k;
	int nprot = 0;

	// THIS CAUSES AN ERROR.
	//PROTECT(input = AS_NUMERIC(input));
	//nprot++;

	// Get dimensions of input matrix so we can create the output matrix.
	// dims[0] is the numbers of A/B markers, dims[1] is the total number of markers.
	int *dims;
	dims = INTEGER(GET_DIM(input));

    // Convert correctionFlag parameter to a boolean.
    float corrFlag;

    corrFlag = *REAL(correctionFlag);

	// create return matrix as a square matrix with both dimensions equal to the number of markers.
	SEXP ans;
	PROTECT(ans = allocMatrix(REALSXP, dims[1], dims[1]));
	nprot++;

	int diff, missing;

	// Shortcut for dims[0] since we will use it a lot.
	int c = dims[0];
	printf("c = %i, dims[1] = %i\n", c, dims[1]);

	// Step through all markers.
	// Start at i = 1 to ignore the marker name column.
	for(i = 0; i < dims[1]-1; i++)
	//for(i = 1; i < 2; i++)
	{
		if(i%100 == 0)
		{
			printf(".");
			if(i % 500 == 0)
				printf(" Processed %i of %i markers\n", i, dims[1]);
		}

		// Step through all markers from i to the end.
		for(j = i+1; j < dims[1]; j++)
		//for(j = 29; j < 30; j++)
		{
			diff = 0;
			missing = 0;

			//printf("HERE WE ARE #2, i, j = %i, %i\n", i, j);

			// Step through each individual A/B values and compare them.
			// Start at k = 3 to ignore the header lines.
			for(k = 0; k < dims[0]; k++)
			{
                //printf("k = %i, marker1 = %i, marker2 = %i\n", k, INTEGER(input)[i*c + k], INTEGER(input)[j*c + k]);
                if(INTEGER(input)[i*c + k] == 0 || INTEGER(input)[j*c + k] == 0)
                {
                    missing++;
                }
                else if(INTEGER(input)[i*c + k] != INTEGER(input)[j*c + k])
                {
                    diff++;
                    //printf("found a difference, diff = %i\n", diff);
                }
			}

			// Put the computed value into the rf matrix.

            // If corrFlag is true, use the data correction term.  This accounts for missing data by splitting the difference between the rf value when all missing calls are the same and the rf value when all missing calls are different.
            if(corrFlag)
                REAL(ans)[i*dims[1] + j] = (diff + (float)missing/2)/(float)c;
            else
                REAL(ans)[i*dims[1] + j] = diff/(float)(c-missing);

			// We also want to put the value in cell [j,i] so that we have the full matrix.
			// We will need the full matrix when finding similar markers (see findSimilars function above).
			REAL(ans)[j*dims[1] + i] = REAL(ans)[i*dims[1] + j];
            //printf("missing = %i, diff = %i, c = %i\n", missing, diff, c);
            //printf("ans = %f\n", REAL(ans)[i*dims[1] + j]);
			//printf("ic + j = %i\n", i*dims[1] + j);
            //return(ans);
		}

        //printf("ans = %f\n", REAL(ans)[i*dims[1] + j]);
        //printf("missing = %i, diff = %i\n", missing, diff);

	}



	UNPROTECT(nprot);
	return ans;
}

// end computeRF() =============================================================================
// ============================================================================================

