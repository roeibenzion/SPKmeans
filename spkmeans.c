#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int main(int argc, char **argv);
void navigator(char* goal, double** mat, int N, int d, int K);
double ** wamF(double **vect, int N, int d);
double** ddgF(double **W, int N);
double** lnormF(double** W, double **D, int N);
double** JacobiF(double** A, int N);
int eigenGap(double* eigenValues, int N);
double** sortMatrixColumns(double** V, int N, double* A);
double** obtainLargestK(double **V, int N, int K);
double** formTfromU(double** U, int N, int K);
double findT(double theta);
double** getATag(double **A,int N, double c, double s, int i, int j);
int checkConverstion(double** ATag, double offSquareA, int N );
double offSquare(double **A, int N);
double sumRow(double* row, int N, int squared);
void matToThePow(double** M, int N, double exponent, int diag);
double** matMul(double **A, double **B, int N);
double** matSum(double** A, double** B, int N, int sign);
double ** getMatrix(int N, int d);
void freeMatrix(double** m, int k);
double *getVector(int N);
void freeVector(double** m, int k);
int compare( const void* a, const void* b);
int searchIndex(double* arr, int N, int val);
void copyRows(double* dst, double* src, int N);
int checkDiag(double** A, int N);
void printM(double** M, int N);
void printV(double* V, int N);
double** getJacobiMatrix(double **V, double**A, int N);


int main(int argc, char **argv)
{
    int i, j, N, d, offset;
    double **X, num;
    FILE *ifp_datapoints;
    char c;

    offset = 0;
    d = 0;
    N = 0;
    ifp_datapoints = NULL;

    if(argc != 3)
    {
        printf("Invalid Input!");
        exit(1);
    }

    ifp_datapoints = fopen("C:\\Users\\roeib\\CLionProjects\\spkmeans\\test_jacobi_5_input.txt", "r");
    if(ifp_datapoints == NULL)
    {
        printf("An error has occured");
        exit(1);
    }
    //get d
    while((c = getc(ifp_datapoints)) != '\n')
    {
        if(c == ',')
            d++;
        offset++;
    }
    d += 1;
    fseek(ifp_datapoints, 0, 0);
    offset = 0;

    //get N
    while(fscanf(ifp_datapoints, "%lf", &num) == 1)
    {
        c = getc(ifp_datapoints);
        if(c == '\n')
            N++;
        offset++;
    }
    fseek(ifp_datapoints,0 , 0);
    i = 0;
    j = 0;
    X = getMatrix(N, d);
    while(fscanf(ifp_datapoints, "%lf", &num) == 1)
    {
        X[i][j] = num;
        c = getc(ifp_datapoints);
        if(c == '\n')
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
    navigator(argv[1], X, N, d, -1);
}
void navigator(char* goal, double** mat, int N, int d, int K)
{
    int i, j;
    double **W, **D, **Lnorm, **temp, **V, *eigenValues, **U, **T;
    char* spk, *wam, *ddg, *lnorm, *jacobi;
    FILE* ifp;

    ifp = fopen("C:\\Users\\roeib\\CLionProjects\\spkmeans\\out.txt", "w");
    spk = "spk";
    wam = "wam";
    ddg = "ddg";
    lnorm = "lnorm";
    jacobi = "jacobi";

    if(!strcmp(goal,wam) ) {
        W = wamF(mat, N, d);
        printM(W, N);
    }
    else if(!strcmp(goal,ddg))
    {
        W = wamF(mat, N, d);
        D = ddgF(W, N);
        printM(D, N);
    }

    else if(!strcmp(goal,lnorm))
    {
        W = wamF(mat, N, d);
        D = ddgF(W, N);
        Lnorm = lnormF(W, D, N);
        printM(Lnorm, N);
    }

    else if(!strcmp(goal,jacobi))
    {
        U = JacobiF(mat, N);
        for(i = 0; i< N; i++)
        {
            printV((U[i]), N);
            printf(("\n"));
        }
        printV((U[N]), N);
    }
    else if(!strcmp(goal, spk))
    {
        W = wamF(mat, N, d);
        D = ddgF(W, N);
        Lnorm = lnormF(W, D, N);
        temp = JacobiF(Lnorm, N);
        copyRows(eigenValues, temp[0], N);
        for(i = 0; i < N; i++)
        {
            copyRows(V[i], temp[i+1], N);
        }
        V = sortMatrixColumns(V, N, eigenValues);
        if(K < 1)
            {
                K = eigenGap(eigenValues, N);
            }
        U = obtainLargestK(V, N, K);
        T = formTfromU(U, N, K);
        //return to python and connect to KMEANS from HW2.
    }
    else
    {
        printf("Invalid Input!");
        exit(1);
    }
}

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
double** wamF(double **vect, int N, int d)
{
    int i, j, k;

    double** W = getMatrix(N, N);
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

double** ddgF(double **W, int N)
{
    int i;
    double sum;

    double** D = getMatrix(N, N);
    for(i = 0; i < N; i++)
    {
        D[i][i] = sumRow(W[i], N, 0);
    }
    return D;
}

double** lnormF(double** W, double **D, int N)
{
    double **I, **L, **C;
    int i, j;
    I = getMatrix(N,N);
    C = getMatrix(N,N);
    for(i = 0; i < N; i++)
    {
        I[i][i] = 1;
    }
    matToThePow(D, N, -0.5, 1);
    C = matMul(D, W, N);
    C = matMul(C, D, N);
    L = matSum(I, C, N, -1);
    return L;
}

//Returns a (N+1)X(N) matrix, where the first row is the eigenValues and
// the rest is a (N)X(N) matrix of eigenVectors.
double** JacobiF(double** A, int N)
{
    int i,j, iMax, jMax, countIter, firstIter;
    double theta, t, c, s, maxVal, offSquareA;
    double ** P, **V;

   V = getMatrix(N, N);
   P = getMatrix(N,N);
   for(i = 0; i < N; i++)
   {
    P[i][i] = 1;
    V[i][i] = 1;
   }

   countIter = 0;
   firstIter = 1;
   maxVal = 0.0;
   iMax = 0;
   jMax = 0;

    printM(A, N);
   while(countIter < 100)
   {
       countIter++;
       if(checkDiag(A, N))
       {
           break;
       }
       for(i = 0; i < N; i++)
        {
        for(j = 0; j < N ; j++)
        {
            if(i == j)
            {
                continue;
            }
            else if(maxVal < A[i][j])
            {
                iMax = i;
                jMax = j;
                maxVal = A[i][j];
            }
        }
    }
    theta = (A[jMax][jMax] - A[iMax][iMax]) / (2*A[iMax][jMax]);
    t = findT(theta);
    c = 1/(sqrt(t*t + 1));
    s = t*c;

    P[iMax][iMax] = c;
    P[jMax][jMax] = c;
    P[iMax][jMax] = s;
    P[jMax][iMax] = -s;

    offSquareA = offSquare(A, N);
    A = getATag(A, N, c, s, iMax, jMax);
    if(checkConverstion(A, offSquareA, N))
    {
        break;
    }
    V = matMul(V, P, N);
   }
    printf("%d\n", countIter);
   return getJacobiMatrix(V, A,N);
}

double** getJacobiMatrix(double **V, double**A, int N)
{
    int i,j;
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
    return U;
}
int eigenGap(double* eigenValues, int N)
{
    double deltaI, temp;
    int i, iter, maxIndex;

    iter = N/2;
    if(iter == 1)
    {
        return 1;
    }
    i = 1;
    deltaI = fabs(eigenValues[i] - eigenValues[i+1]);
    maxIndex = 1;
    for(i = 2; i < iter; i++)
    {
        temp = fabs(eigenValues[i] - eigenValues[i+1]);
        if(deltaI < temp)
        {
            maxIndex = i;
            deltaI = temp;
        }
    }
    return maxIndex;
}

double** sortMatrixColumns(double** V, int N, double* A)
{
    int i, j, row;
    double* copy, **newV;

    copy = getVector(N);
    newV = getMatrix(N, N);
    for(i = 0; i < N; i++)
    {
        copy[i] = A[i];
    }
    qsort(copy,N,sizeof(double),compare);

    for(i = 0; i < N; i++)
    {
        row = searchIndex(copy, N, A[i]);
        copyRows(newV[i], V[row], N);
    }
    return newV;
}

double** obtainLargestK(double **V, int N, int K)
{
    int i, j;
    double** U;

    U = getMatrix(N, K);
    for(i = 0; i < K; i++)
    {
       for(j = 0; j < N; j++)
       {
        U[i][j] = V[i][j];
       }
    }
    return U;
}
double** formTfromU(double** U, int N, int K)
{
    int i, j;
    double sum, **T;

    T = getMatrix(N, K);
    for(i = 0; i < N; i++)
    {
        sum = sqrt(sumRow(U[i], K, 1));
        for(j = 0; j < K; j++)
        {
            T[i][j] = U[i][j] / sum;
        }
    }
    return T;
}
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
    base = fabs(theta) + sqrt((pow(theta,2)) + 1);
    t = sign / base;
    return t;
}

double** getATag(double **A,int N, double c, double s, int i, int j)
{
    int r;
    double **ATag;

    ATag = getMatrix(N,N);
    for(r = 0; r < N; r++)
    {
        if(r != j && r != i)
        {
            ATag[r][i] = c*A[r][i] - s*A[r][j];
            ATag[i][r] = ATag[r][i];
            ATag[r][j] = c*A[r][j] + s*A[r][i];
            ATag[j][r] = ATag[r][j];
        }
    }
    ATag[i][i] = pow(c,2)*(A[i][i]) + pow(s,2)*(A[j][j]) - 2*s*c*A[i][j];
    ATag[j][j] = pow(c,2)*(A[i][i]) + pow(c,2)*(A[j][j]) + 2*s*c*A[i][j];
    ATag[i][j] =  0;
    ATag[j][i] = 0;
    return ATag;
}

int checkConverstion(double** ATag, double offSquareA, int N )
{
    int r;
    double offSquareATag, eps;

    eps = pow(10,-5);
    offSquareATag = offSquare(ATag, N);
    return (offSquareA - offSquareATag) <= eps;
}

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

double** matMul(double **A, double **B, int N)
{
    int i,j,k;
    double **C, sum;
    C = getMatrix(N,N);
    for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {
        sum = 0.0;
        for (k = 0; k < N; k++) {
          sum = sum + A[i][k]*B[k][j];
        }
        C[i][j] = sum;
      }
    }
    return C;
}

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

void freeMatrix(double** m, int k) {
    int i;
    for (i=0; i<k; i++) {
        free(m[i]);
    }
    free(m);
}

double * getVector(int N)
{
    double *data;
    int i;
    data = (double*) calloc(N, sizeof(double*));
    if(data == NULL)
        return data;
    return data;
}

void freeVector(double** m, int k) {
    free(m);
}

int compare( const void* a, const void* b)
{
     double double_a, double_b;
     double_a = *((double*)a);
     double_b = *((double*)b);
     if ( double_a == double_b ) return 0;
     else if ( double_a < double_b ) return -1;
     else return 1;
}

int searchIndex(double* arr, int N, int val)
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

void copyRows(double* dst, double* src, int N)
{
    int i;
    for(i = 0; i < N ; i++)
    {
        dst[i] = src[i];
    }
}

int checkDiag(double** A, int N)
{
    int i, j;

    for(i = 0; i < N; i++)
    {
        for(j = 0; j < N; j++)
        {
            if(j != i && A[i][j] != 0)
            {
                return 0;
            }
        }
    }
    return 1;
}

void printM(double** M, int N)
{
    int i, j;
    for(i = 0; i < N ; i++)
    {
        for(j = 0; j < N-1; j++)
        {
            printf("%0.4f,", M[i][j]);
        }
        printf("%0.4f\n", M[i][N-1]);
    }
}

void printV(double* V, int N)
{
    int i;
    for(i = 0; i < N-1 ; i++)
    {
        printf("%0.4f,", V[i]);
    }
    printf("%0.4f", V[N-1]);
}
