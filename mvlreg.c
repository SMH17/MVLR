/****************************************************************************
 *
 * Multivariate Linear Regression
 *
 * S. Marano (2013)
 *
 ****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <xmmintrin.h>




float frand() {
    float r = (float) rand();
    return r/RAND_MAX;
}


/*
 *
 * 	Matrix Methods needed to calculate the Signed Beta
 * 	============
 *
 *  Signed Beta is used to determine the coefficients that reduce the extimation 
 *  square error err(B) and is defined in matrix form as ((X^T * X)^-1)X^T*y.
 *
 */


float **allocate_mem_matrix(int m, int n) {
    printf("\nMEMORY ALLOCATION...\n");
    int i;
    float **memoryArray = malloc(m*sizeof(*memoryArray));
    for(i=0; i<m; i++) {
        memoryArray[i]=malloc(n*sizeof(**memoryArray));
    }
    return memoryArray;
}

void free_mem_matrix(int m, float **memoryArray) {
   printf("\nMEMORY DEALLOCATION...\n");
   int i;
    for (i = 0; i < m; i++) {
        free(memoryArray[i]);
    }
    free(memoryArray);
}

float **transposeMatrix(int m, int n, float **matrix) {
    printf("\nTRANSPOSE MATRIX...\n");
    int i,j;
    float **result = allocate_mem_matrix(n,m);
    for(i = 0; i < m; i++) {
        for(j = 0; j < n; j++) {
            result[j][i] = matrix[i][j];
        }
    }
    return result;
}

float **multiplyMatrix(int m1, int n1, float **matrix1,int m2, int n2, float **matrix2) {
    printf("\nMULTIPLY MATRICES...\n");
    int i,j,k;
    float **result = allocate_mem_matrix(m1,n2);
    for (i = 0; i < m1; i++) {
        for (j = 0; j < n2; j++) {
            for (k = 0; k < n1; k++) {
                result[i][j] = result[i][j] + matrix1[i][k] * matrix2[k][j];
            }
        }
    }
    return result;
}

float **identityMatrix(int n) {
    printf("\nIDENTITY MATRIX CREATION...\n");
    int i,j;
    float **b= allocate_mem_matrix(n,n);
    for(i=0; i<n; i++) {
        for(j=0; j<n; j++) {
            if(i==j) {
                b[i][j]=1;
            }
            else {
                b[i][j]=0;
            }
        }
    }
    return b;
}

float **invertMatrix(int n, float **matrix) {
    printf("\nINVERT MATRIX...\n");
    float tmp0=0,tmp1=0,tmp2=0,tmp3=0,tmp4=0,tmp5=0;
    int m=0,i=0,j=0,p=0,q=0;

    //get an identity matrix of the same size
    float **a= identityMatrix(n);

    for(i=0; i<n; i++) {
        tmp1=matrix[i][i];
        if(tmp1<0)
            tmp1=tmp1*(-1);
        p=i;
        for(j=i+1; j<n; j++) {
            if(matrix[j][i]<0)
                tmp0=matrix[j][i]*(-1);
            else
                tmp0=matrix[j][i];
            if(tmp1<0)
                tmp1=tmp1*(-1);
            if(tmp0>tmp1) {
                p=j;
                tmp1=matrix[j][i];
            }
        }
        //exchange rows in both matrices
        for(j=0; j<n; j++) {
            tmp2=matrix[i][j];
            matrix[i][j]=matrix[p][j];
            matrix[p][j]=tmp2;
            tmp3=a[i][j];
            a[i][j]=a[p][j];
            a[p][j]=tmp3;
        }
        //divide row by matrix[i][i]]
        tmp4=matrix[i][i];
        if(tmp4==0) {
            printf("\nMATRIX CANNOT BE INVERTED!!!\n");
            return NULL;
        }
        for(j=0; j<n; j++) {
            matrix[i][j]=(float)matrix[i][j]/tmp4;
            a[i][j]=(float)a[i][j]/tmp4;
        }
        //turn matrix[][] in an identity matrix in order to get the inverse of a[][]
        for(q=0; q<n; q++) {
            if(q==i)
                continue;
            tmp5=matrix[q][i];
            for(j=0; j<n; j++) {
                matrix[q][j]=matrix[q][j]-(tmp5*matrix[i][j]);
                a[q][j]=a[q][j]-(tmp5*a[i][j]);
            }
        }
    }
    return a;
}

float **getLastColumn(int m, int n, float **matrix) {
    int i,j;
    float **y= allocate_mem_matrix(m,1);
    for(i=0; i<m; i++) {
        y[i][0]=matrix[i][n-1];
    }
    return y;
}

/*
 *
 * 	random_input
 * 	============
 *
 *	Generate a random matrix m x (n+1) to use as input for linreg.
 *
 */

float **random_input(int m, int n) {
    int i,j;
    float x, y;
    float **a= allocate_mem_matrix(m,n);

    for(i=0; i<m; i++) {
        y = 0;
        for(j=0; j<n; j++) {
            if(j < n-1) {
                a[i][j]= frand();
            }
            else {
                a[i][j]=1.0;
                y+=frand();
            }
        }
        a[i][j]=y*(1+frand());
    }
    return a;
}

/*
 *
 * 	load_input
 * 	===========
 *
 *	Load matrix m x (n+1) from file.
 * 	The matrix is used as input for linreg.
 *  (Note: The matrix is assumed as row-major order)
 *
 */
float **load_input(char* filename, int *m, int *n) {
    FILE* fp;
    int rows, cols, status;

    fp = fopen(filename, "rb");
    status = fread(&cols, sizeof(int), 1, fp);
    status = fread(&rows, sizeof(int), 1, fp);
    float **Xy = allocate_mem_matrix(rows,cols);
    status = fread(Xy, sizeof(float), rows*cols, fp);
    fclose(fp);

    *m = rows;
    *n = cols-1;
    return Xy;
}


/*
 *
 * 	load_input
 * 	===========
 *
 *	Write matrix m x (n+1).
 *  (Note: The matrix is assumed as row-major order)
 *
 */
void save_input(char* filename, float **Xy, int m, int n) {
    FILE* fp;
    int status;

    fp = fopen(filename, "wb");
    n++;
    status = fwrite(&n, sizeof(int), 1, fp);
    status = fwrite(&m, sizeof(int), 1, fp);
    status = fwrite(Xy, sizeof(float), m*n, fp);
    fclose(fp);
}


/*
 *	linreg
 * 	======
 *
 * Xy is a matrix of dimension m x (n + 1) in which the first n
 * columns represent the independent variables x1, ..., xn
 * and the last column represents the dependent variable y;
 * each row contains an observation (in the first n + 1 columns)
 * and the associated value of the dependent variable (in the last column).
 *
 * The output is an array of n columns containing the coefficients
 * of the regression hyperplane
 *
 */


float **linreg(float **Xy, int m, int n) {
    printf("\nSIGNED BETA CALCULATION...\n");
    //create columnsX to select n-1 column of X
    int columnsX=n-1;
    //extract Y as latest column of Xy
    float **Y = getLastColumn(m, n, Xy);
    float **transposedmatr = transposeMatrix(m, columnsX, Xy);
    float **prod = multiplyMatrix(columnsX, m, transposedmatr, m, columnsX, Xy); 
    float **temp = invertMatrix(columnsX, prod);
    float **temp2 = multiplyMatrix(columnsX,columnsX,temp, columnsX,m,transposedmatr);
    float **beta=multiplyMatrix(columnsX,m,temp2, m,1,Y);
    //free no more used memory
    free_mem_matrix(columnsX, transposedmatr);
    free_mem_matrix(columnsX, temp);
    free_mem_matrix(columnsX, temp2);
    free_mem_matrix(m, Y);
    return beta;
}

/*
 *
 *	Calculate regression error.
 *
 */
float error(float **Xy, float **beta, int m, int n) {
    int i, j;
    float err = 0, de, yp;

    for (i = 0; i < m; i++) {
        yp = 0;
        for (j = 0; j < n-1; j++)
            yp += Xy[i][j] * beta[j][0];
        de = fabs(Xy[i][j]-yp);
        err += (de*de);
    }
    return err/(m*n);
}


/*
 *
 *Test performances.
 *
 */
void main(int argc, char** argv) {
    int m = 1000;
    int n = 1000;
    float **Xy;
    float **beta;

    char* filename = "";
    int silent = 0, display = 0, fsave = 0;
    int i;

    srandom(time(NULL));

    int par = 1;
    while (par < argc) {
        if (strcmp(argv[par],"-l") == 0) {
            par++;
            if (par < argc) {
                filename = argv[par];
                par++;
            }
        } else if (strcmp(argv[par],"-r") == 0) {
            par++;
            if (par < argc) {
                m = atoi(argv[par]);
                par++;
                if (par < argc) {
                    n = atoi(argv[par]);
                    par++;
                }
            }
        } else if (strcmp(argv[par],"-f") == 0) {
            fsave = 1;
            par++;
        } else if (strcmp(argv[par],"-s") == 0) {
            silent = 1;
            par++;
        } else if (strcmp(argv[par],"-d") == 0) {
            display = 1;
            par++;
        } else
            par++;
    }

    if (!silent) {
        printf("Usage: %s [-l <file_name>][-r <observations(m)> <variables(n)>]\n", argv[0]);
        printf("\nParameters:\n");
        printf("\t-l <file_name> : reads from disk the augmented m x (n+1) matrix\n");
        printf("\t-r <observations(m)> <variables(n)>: randomly generates a m x (n+1) matrix\n");
        printf("\t-f : writes on disk the augmented m x (n+1) matrix (creates the file last.mat)\n");
        printf("\t-d : displays input and output\n");
        printf("\t-s : silent\n");

        printf("\nSolving a regression problem with %d observations and %d variables...\n", m, n);
    }

    if (strlen(filename) == 0)
        Xy = random_input(m,n);
    else
        Xy = load_input(filename, &m, &n);

    if (!silent && display) {
        printf("\nInput augmented matrix:\n");
        for (i = 0; i < m*(n+1); i++) {
            if (i % (n+1) == 0)
                printf("\n");
            printf("%f ", Xy[i]);
        }
        printf("\n");
    }

    clock_t t = clock();
    beta = linreg(Xy,m,n);
    t = clock() - t;

    if (!silent)
        printf("\nExecution time = %.3f seconds\n", ((float)t)/CLOCKS_PER_SEC);
    else
        printf("%.3f\n", ((float)t)/CLOCKS_PER_SEC);

    if (!silent && display) {
        printf("\nOutput coefficient vector:\n");
        for (i = 0; i < n; i++)
            printf("%.3f ", beta[i]);
        printf("\n");
    }

    float err = error(Xy,beta,m,n);
    if (!silent)
        printf("\nThe error is %f.\n", err);
    else
        printf("%f\n", err);

    if (strlen(filename) == 0 && fsave)
        save_input("last.mat",Xy,m,n);
}