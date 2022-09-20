#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void generateId(double** A, int N);
double sumRow(double* row, int N, int squared);
void matToThePow(double** M, int N, double exponent, int diag);
double** matSum(double** A, double** B, int N, int sign);
double** getMatrix(int N, int d);
double* getVector(int N);
void freeMatrix(double** m, int k);
void freeVector(double* m);
int searchIndex(double* arr, int N, double val);
void copyRows(double* dst, double* src, int N);
void printM(double** M, int N);
double euclidSum(double *x1, double*x2, int d);
double** wamF(double **vect, int N, int d);
double** ddgF(double **W, int N);
double** lnormF(double** W, double **D, int N);
double** getJacobiMatrix(double **V, double**A, int N);
double** getATag(double **A,int N, double c, double s, int i, int j);
int checkConversion(double** ATag, double offSquareA, int N );
double offSquare(double **A, int N);
double findT(double theta);
double** JacobiF(double** A, int N);
int eigenGap(double* eigenValues, int N);
void sortMatrixColumns(double** V, int N, double* A);
double** obtainLargestK(double **V, int N, int K);
void formTfromU(double** U, int N, int K);
int main(int argc, char **argv);
void navigator(char* goal, double** mat, int N, int d, double** ret);
void copyCol(double** dst, double** src, int N, int dstIndex, int srcIndex);
int compare( const void* a, const void* b);

/*Utility functions*/
/*Get matrix A and turnd it to the identity matrix*/
void generateId(double** A, int N)
{
    int i, j;
    for(i = 0; i < N; i++)
    {
        for(j = 0; j < N; j++)
        {
            if(i == j)
                A[i][j] = 1;
            else
                A[i][j] = 0;
        }
    }
}

/*Sum an array. If squared = 1 then sum it's squeres*/
double sumRow(double* row, int N, int squared)
{
    int i;
    double sum;
    sum = 0;
    if(squared == 0)
    {
        for(i = 0; i < N; i++)
        {
            sum += row[i];
        }
    }
    else
    {
        for(i = 0; i < N; i++)
        {
            sum += pow(row[i],2);
        }
    }
    return sum;
}
/*Take matrix M to the pow exponent. If M is a diagonal matrix then refers only to the diagonal*/
void matToThePow(double** M, int N, double exponent, int diag)
{
    int i, j;
    if(diag == 0)
    {
        for(i = 0; i < N; i++)
        {
            for(j = 0; j < N; j++)
            {
                M[i][j] = pow(M[i][j], exponent);
            }
        }
    }
    else
    {
        for(j = 0; j < N; j++)
        {
            M[j][j] = pow(M[j][j], exponent);
        }
    }
}

/*Sums matrix A and B. if Sign = 1 regualr sum, else subtract.*/
double** matSum(double** A, double** B, int N, int sign)
{
    int i,j;
    double **C;
    C = getMatrix(N,N);
    for(i = 0; i < N; i++)
    {
        for(j = 0; j < N; j++)
        {
            if(sign == 1)
                C[i][j] = A[i][j] + B[i][j];
            else
                C[i][j] = A[i][j] - B[i][j];
        }
    }
    return C;
}
/*Generates an NXd matrix of doubles*/
double ** getMatrix(int N, int d)
{
    double **data;
    int i, j;
    data = (double**) calloc(N, sizeof(double*));
    if(data == NULL)
        return data;
    for(i = 0; i < N; i++)
    {
        data[i] = (double*)calloc(d,sizeof(double));
        if(data[i] == NULL)
        {
            for(j = 0 ;j <= i-1; j++)
            {
                free(data[j]);
            }
            free(data);
            return NULL;
        }
    }
    return data;
}
/*Free a matrix with k rows*/
void freeMatrix(double** m, int k) {
    int i;
    for (i=0; i<k; i++) {
        free(m[i]);
    }
    free(m);
}
/*Generated an N-dimentional vector*/
double * getVector(int N)
{
    double *data;
    data = (double*) calloc(N, sizeof(double*));
    if(data == NULL)
        return data;
    return data;
}
/*Free a vector*/
void freeVector(double* m) {
    free(m);
}
/*Search the index of val in arr*/
int searchIndex(double* arr, int N, double val)
{
    int i;
    for(i = 0; i<N; i++)
    {
        if(arr[i] == val)
        {
            return i;
        }
    }
    return -1;
}
/*Copy row src to row dst*/
void copyRows(double* dst, double* src, int N)
{
    int i;
    for(i = 0; i < N ; i++)
    {
        dst[i] = src[i];
    }
}

/*Copying cols from matrix to matrix*/
void copyCol(double** dst, double** src, int N, int dstIndex, int srcIndex)
{
    int i;
    for(i = 0; i < N ; i++)
    {
        dst[i][dstIndex] = src[i][srcIndex];
    }
}
/*Prints a matrix in the desired format*/
void printM(double** M, int N)
{
    int i, j;
    for(i = 0; i < N; i++)
    {
        for(j = 0; j < N-1; j++)
        {
            printf("%.4f,", M[i][j]);
        }
        printf("%.4f\n", M[i][N-1]);
    }
}
/*Euclid sum of 2  vectors in d-dimention*/
double euclidSum(double *x1, double*x2, int d)
{
    int i;
    double sum;
    sum = 0.0;
    for( i = 0; i < d; i++)
    {
        sum += pow((x1[i] - x2[i]), 2);
    }
    return sqrt(sum);
}
/*End utility functions*/

/*Wam*/
double** wamF(double **vect, int N, int d)
{
    int i, j;
    double** W = getMatrix(N, N);
    if(W == NULL)
    {
        printf("An error has occurred");
        exit(1);
    }

    for(i = 0; i < N; i++)
    {
        for(j = 0; j < N; j++)
        {
            if(j == i)
            {
                W[i][j] = 0;
            }
            else
            {
                W[i][j] = exp((-(euclidSum(vect[i], vect[j], d)/2)));
            }
        }
    }
    return W;

}
/*Ddg*/
double** ddgF(double **W, int N)
{
    int i;

    double** D = getMatrix(N, N);
    if(D == NULL)
    {
        printf("An error has occurred");
    }
    for(i = 0; i < N; i++)
    {
        D[i][i] = sumRow(W[i], N, 0);
    }
    return D;
}
/*Lnorm*/
double** lnormF(double** W, double **D, int N)
{
    double **I, **L; 
    int i, j;
    I = getMatrix(N,N);
    if(I == NULL)
    {
        printf("An error has occurred");
        exit(1);
    }
    generateId(I, N);
    matToThePow(D, N, -0.5, 1);
   /*From the symmetry of D, we can switch the matMul D^-0.5 * W * D^-0.5 operations to something more efficient*/
    for(i = 0; i < N; i++)
    {
        for(j = 0; j < N; j++)
        {
            W[i][j] *= D[i][i];
            W[i][j] *= D[j][j];
        }
    }

    L = matSum(I, W, N, -1);
    freeMatrix(I, N);
    freeMatrix(W, N);
    freeMatrix(D, N);
    return L;
}

/*Jacobi functions.
Returns a (N+1)X(N) matrix, where the first row is the eigenValues and
the rest is a (N)X(N) matrix of eigenVectors.
*/
/*Union the output of jacobiF to (N+1)X(N) matrix*/
double** getJacobiMatrix(double **V, double**A, int N)
{
    int i;
    double** U;

    U = getMatrix(N+1, N);
    for(i = 0; i < N; i++)
    {
        U[0][i] = A[i][i];
    }
    for(i = 1; i < N+1; i++)
    {
        copyRows(U[i], V[i-1], N);
    }
    freeMatrix(V, N);
    freeMatrix(A, N);
    return U;
}

/*Fast way to get A'*/
double** getATag(double **A,int N, double c, double s, int i, int j) {
    int r;
    double **ATag;
    ATag = getMatrix(N, N);
    for (r = 0; r < N; r++)
    {
        copyRows(ATag[r], A[r], N);
    }
    for(r = 0; r < N; r++)
    {
        if(r != j && r != i)
        {
            ATag[r][i] = (c*A[r][i]) - (s*A[r][j]);
            ATag[i][r] = ATag[r][i];
            ATag[r][j] = (c*A[r][j]) + (s*A[r][i]);
            ATag[j][r] = ATag[r][j];
        }
    }
    ATag[i][i] = (c*c)*(A[i][i]) + (s*s)*(A[j][j]) - 2*(s*c*A[i][j]);
    ATag[j][j] = (s*s)*(A[i][i]) + (c*c)*(A[j][j]) + 2*(s*c*A[i][j]);
    ATag[i][j] =  0.0;
    ATag[j][i] = 0.0;

    for (r = 0; r < N; r++)
    {
        copyRows(A[r], ATag[r], N);
    }
    freeMatrix(ATag, N);
    return A;
}
/*Check conversion*/
int checkConversion(double** ATag, double offSquareA, int N )
{
    double offSquareATag, eps;

    eps = pow(10,-5);
    offSquareATag = offSquare(ATag, N);
    return (offSquareA - offSquareATag) <= eps;
}

/*Sum of off diag elements*/
double offSquare(double **A, int N)
{
    int i, j;
    double sum;

    sum = 0.0;
    for(i = 0; i < N; i++)
    {
        for(j = 0; j < N; j++)
        {
            if(i != j)
            {
                sum += pow(A[i][j], 2);
            }
        }
    }
    return sum;
}
/*Find the t parameter from theta to construct c and s*/
double findT(double theta)
{
    int sign;
    double base, t;
    if(theta < 0)
    {
        sign = -1;
    }
    else
    {
        sign = 1;
    }
    base = fabs(theta) + sqrt(((theta*theta) + 1));
    t = sign / base;
    return t;
}
/*Jacobi algorithm*/
double** JacobiF(double** A, int N)
{
    int i,j, iMax, jMax, countIter;
    double theta, t, c, s, maxVal, offSquareA;
    double ** P, **V, temp1, temp2;

    V = getMatrix(N, N);
    P = getMatrix(N,N);
    generateId(V, N);
    countIter = 0;

    while(countIter < 100)
    {
        maxVal = 0.0;
        iMax = 0;
        jMax = 0;
        countIter++;
        /*Find pivot*/
        for(i = 0; i < N; i++)
        {
            for(j = i+1; j < N ; j++)
            {
                if(maxVal < fabs(A[i][j]) && i!=j)
                {
                    iMax = i;
                    jMax = j;
                    maxVal = fabs(A[i][j]);
                }
            }
        }
        /*Couldn't find a non-zero pivot implies matrix is diagonal*/
        if(maxVal == 0.0)
        {
            break;
        }
        /*Find c and s*/
        theta = (A[jMax][jMax] - A[iMax][iMax]) / (2*A[iMax][jMax]);
        t = findT(theta);
        c = 1/(sqrt(t*t + 1));
        s = t*c;

        /*Generate P*/
        generateId(P, N);
        P[iMax][iMax] = c;
        P[jMax][jMax] = c;
        P[iMax][jMax] = s;
        P[jMax][iMax] = -s;

        /*Update V. From the special structure of P we can switch the matMul operation with a closed formula*/
        for(i = 0; i < N; i++)
        {
            temp1 = V[i][iMax];
            temp2 = V[i][jMax];
            V[i][iMax] = c*temp1 - s*temp2;
            V[i][jMax] = c*temp2 + s*temp1;
        }
        /*Check conversion*/
        offSquareA = offSquare(A, N);
        A = getATag(A, N, c, s, iMax, jMax);
        if(checkConversion(A, offSquareA, N))
        {
            break;
        }
    }
    freeMatrix(P, N);
    /*Constract (N+1)X(N) matrix from the diag of A and V*/
    return getJacobiMatrix(V,A,N);
}
/*End of Jacobi*/
 
/*Full spk related functions*/
/*EigenGap huristic*/
int eigenGap(double* eigenValues, int N)
{
   double deltaI, temp;
    int i, iter, maxIndex;

    iter = N/2;
    if(N == 1)
    {
        printf("An error has occured");
        exit(1);
    }
    maxIndex = 1;
    deltaI = fabs(eigenValues[0] - eigenValues[1]);
    for(i = 1; i < iter; i++)
    {
        temp = fabs(eigenValues[i] - eigenValues[i+1]);
        if(deltaI < temp)
        {
            deltaI = temp;
            maxIndex = i+1;
        }
    }
    return maxIndex;
}
/*Comparator function for qsort*/
int compare( const void* a, const void* b)
{
    double double_a, double_b;
    double_a = *((double*)a);
    double_b = *((double*)b);
    if ( double_a == double_b ) return 0;
    else if ( double_a < double_b ) return 1;
    else return -1;
}
/*Sort matrix V's cols by corresponding eigenvalues in A*/
void sortMatrixColumns(double** V, int N, double* A)
{
    int i, col;
    double* copy, **newV;

    copy = getVector(N);
    newV = getMatrix(N, N);

    /*Generate a sorted copy of the eigenValues vector*/
    copyRows(copy, A, N);
    qsort(copy,N,sizeof(double),compare);

    /*Sort V cols using the sorted copy vector*/
    for(i = 0; i < N; i++)
    {
        /*find ith largest eigenValue col index in the original A vector*/
        col = searchIndex(A, N, copy[i]);
        if(col == -1)
        {
            printf("An error has occurred");
            exit(1);
        }
        /*copy the col of the ith largest value to the ith location*/
        copyCol(newV, V, N, i, col);
    }

    /*The desired sorted V is in newV matrix so copy value by value*/
    for(i = 0; i < N; i++)
    {
        copyRows(V[i], newV[i], N);
    }
    freeMatrix(newV, N);
    freeVector(copy);
}

/*Generate U, assume V is sorted in decreasing order by eigenValues*/
double** obtainLargestK(double **V, int N, int K)
{
     int i;
    double** U;

    U = getMatrix(N, K);
    for(i = 0; i < N; i++)
    {
        copyRows(U[i], V[i], K);
    }
    freeMatrix(V, N);
    return U;
}
/*Forming T from matrix U*/
void formTfromU(double** U, int N, int K)
{
    int i, j;
    double sum;

    for(i = 0; i < N; i++)
    {
        sum = sqrt(sumRow(U[i], K, 1));
        /*sum == 0.0 implies row of 0's so we can do nothing*/
        if(sum != 0.0)
        {
            for(j = 0; j < K; j++)
            {
                U[i][j] = U[i][j] / sum;
            }
        }
    }
}
/*main function*/
int main(int argc, char **argv)
{
    int i, j, N, d;
    double **X, num, **ret;
    FILE *ifp_datapoints;
    char c;
    const char *jacobi;

    jacobi = "jacobi";
    d = 1;
    N = 0;
    ifp_datapoints = NULL;
    i = 0;
    j = 0;

    if(argc != 3)
    {
        printf("Invalid Input!");
        exit(1);
    }
    ifp_datapoints = fopen(argv[2], "r");
    /*Read file. get d (dimention) and N (num of datapoints)*/
    if(ifp_datapoints == NULL)
    {
        printf("An error has occurred");
        exit(1);
    }
   for (c = getc(ifp_datapoints); c != '\n'; c = getc(ifp_datapoints))
   {
    if(c == ',')
        d = d+1;
   }
    fseek(ifp_datapoints, 0, 0);
   for (c = getc(ifp_datapoints); c != EOF; c = getc(ifp_datapoints))
   {
        if (c == '\n')
        {
            N = N + 1;
        }
   }
    fseek(ifp_datapoints,0 , 0);
    X = getMatrix(N, d);

    i = 0;
    j = 0;
    while(fscanf(ifp_datapoints, "%lf", &num) == 1)
    {
        X[i][j] = num;
        c = getc(ifp_datapoints);
        if(c != ',')
        {
            i++;
            j = 0;
        }
        else
            j++;
    }  
    if(fclose(ifp_datapoints) != 0)
    {
        printf("An Error Has Occurred");
        exit(1);
    }
    if(!strcmp(argv[1], jacobi))
    {
        ret = getMatrix(N+1, N);
    }
    else
    {
        ret = getMatrix(N, N);
    }

    /*call navigator function*/
    navigator(argv[1], X, N, d, ret);

    if(!strcmp(argv[1], jacobi))
    {
        for(i = 0; i < N+1; i++)
        {
            for(j = 0; j < N-1; j++)
            {
                printf("%.4f,", ret[i][j]);
            }
            printf("%.4f\n", ret[i][N-1]);
        }
        freeMatrix(ret, N+1);
    }
    else
    {
        printM(ret, N);
        freeMatrix(ret, N);
        freeMatrix(X, N);
    }
    return 0;
}

/*Navigate through goals*/
void navigator(char* goal, double** mat, int N, int d, double** ret)
{
    int i;
    double **W, **D, **Lnorm, **temp;
    const char *wam, *ddg, *lnorm, *jacobi;
    wam = "wam";
    ddg = "ddg";
    lnorm = "lnorm";
    jacobi = "jacobi";

    if(!strcmp(goal,wam) ) {
        W = wamF(mat, N, d);
        for(i = 0; i < N; i++)
        {
            copyRows(ret[i], W[i], N);
        }
        freeMatrix(W, N);
    }
    else if(!strcmp(goal,ddg))
    {
        W = wamF(mat, N, d);
        D = ddgF(W, N);
         for(i = 0; i < N; i++)
        {
            copyRows(ret[i], D[i], N);
        }
        freeMatrix(W, N);
        freeMatrix(D, N);
    }

    else if(!strcmp(goal,lnorm))
    {
        W = wamF(mat, N, d);
        D = ddgF(W, N);
        Lnorm = lnormF(W, D, N);
         for(i = 0; i < N; i++)
        {
            copyRows(ret[i], Lnorm[i], N);
        }
        freeMatrix(Lnorm, N);
    }
    else if(!strcmp(goal,jacobi))
    {
        temp = JacobiF(mat, N);
         for(i = 0; i < N+1; i++)
        {
            copyRows(ret[i], temp[i], N);
        }
        freeMatrix(temp, N+1);
    }
    else
    {
        printf("Invalid Input!");
        exit(1);
    }
}