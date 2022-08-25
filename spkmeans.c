#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int main(int argc, char **argv)
{
    int i, j, N, d, offset;
    double **X, num;
    FILE *ifp_datapoints;
    char*c;

    offset = 0;
    d = 0;
    N = 0; 
    ifp_datapoints = NULL;

    if(argc != 2)
    {
        printf("Invalid Input!");
        exit(1);
    }

    ifp_datapoints = fopen(argv[2], "r");
    //get d
    while((c = getc(ifp_datapoints)) != '\n')
    {
        if(c == ',')
            d++;
        offset++;
    }
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

    X = getMatrix(N, d);
    for(i = 0; i < N; i++)
    {
        for(j = 0; j < d; j++)
        {
            fscanf(ifp_datapoints, "%lf", &num);
            X[i][j] = num;
        }
    }
    navigator(argv[1], X, N, d, -1, 0);
}
void navigator(char* goal, double** mat, int N, int d, int K, int pyOrC)
{
    int i;
    double **W, **D, **Lnorm, **temp, **V, *eigenValues, **U, **T;
    char* spk, *wam, *ddg, *lnorm, *jacobi;

    spk = "spk";
    wam = "wam";
    ddg = "ddg";
    lnorm = "lnorm";
    jacobi = "jacobi";

    if(strcmp(goal,wam))
    {
        W = wamF(mat, N, D);
        if(pyOrC)
        {
            //return to python
        }
        else
        {
            printM(W, N);
        }
    }

    else if(strcmp(goal,ddg))
    {
        W = wamF(mat, N, D);
        D = ddgF(W, N);
        if(pyOrC)
        {
            //return to python
        }
        else
        {
            printM(D, N);
        }
    }

    else if(strcmp(goal,lnorm))
    {
        W = wamF(mat, N, D);
        D = ddgF(W, N);
        Lnorm = lnormF(W, D, N);
        if(pyOrC)
        {
            //return to python
        }
        else
        {
            printM(Lnorm, N);
        }
    }

    else if(strcmp(goal,jacobi))
    {
        temp = JacobiF(mat, N);
        if(pyOrC)
        {
            //return to python
        }
        else
        {
            printV(temp[0], N);
            U = getMatrix(N, N);
            for(i = 0; i < N; i++)
            {
                copyRows(U[i], temp[i+1], N);
            }
            printM(U, N);
        }
    }
    else if(strcmp(goal, spk))
    {
        W = wamF(mat, N, D);
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

double ** wamF(double **vect, int N, int d)
{
    int i, j;
    double c;
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
                c = sumVectorsDiff(vect[i], vect[j], d);
                W[i][j] = exp(-(c/2));
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
    double ** P, **V, **eigenValuesVectors;

    maxVal = 0.0;
    countIter = 0;
    firstIter = 1;
    P = getMatrix(N, N);
    for(i = 0; i < N; i++)
        {
            P[i][i] = 1;
        }

    while(countIter < 100)
    {
        countIter++;
        if(checkDiag(A, N) == 1)
            {
                break;
            }
        for(i = 0; i < N; i++)
        {
            for(j = 0; j < N; j++)
            {
                maxVal = getMax(maxVal, abs(A[i][j])); 
                iMax = i;
                jMax = j;
            }
        }
        theta = (A[jMax][jMax] - A[iMax][iMax]) / (2*(A[iMax][jMax]));
        t = findT(theta);
        c = 1/(sqrt(pow(t,2) + 1));
        s = t*c;
        P[iMax][jMax] = s;
        P[jMax][iMax] = -s;
        P[iMax][iMax] = c;
        P[jMax][jMax] = c;
        if(firstIter == 1)
        {
            firstIter = 0;
            V = P;
        }
        else
        {
            V = matMul(V, P, N);
        }
        offSquareA = offSquare(A, N, iMax, jMax);
        A = getATag(A, N, c, s, iMax, jMax);
        if(checkConverstion(A,offSquareA, N, iMax, jMax))
        {
            break;
        }
    }
    eigenValuesVectors = getMatrix(N+1,N);
    for(i = 0; i < N; i++)
    {
        eigenValuesVectors[0][i] = A[i][i];
    }
    for(i = 1; i < N+1; i++)
    {
        copyRows(eigenValuesVectors[i], V[i-1], N);
    }
    return eigenValuesVectors;
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
    deltaI = abs(eigenValues[i] - eigenValues[i+1]);
    maxIndex = 1;
    for(i = 2; i < iter; i++)
    {
        temp = abs(eigenValues[i] - eigenValues[i+1]);
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
    base = abs(theta) + sqrt((pow(theta,2)) + 1);
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
    ATag[i][j] = 0;
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
        for(j = i+1; j < N; j++)
        {
            sum += 2*pow(A[i][j],2);
        }
    }
    return sum;
}
double sumVectorsDiff(double* x1, double* x2, int d)
{
    int i;
    double c, sum;
    
    for(i = 0; i < d; i++)
    {
        c = x1[i] - x2[i];
        c = pow(c,2);
        sum += c;
    }
    sum = sqrt(sum);
    return sum;
}

double sumVectors(double* x1, double* x2, int d, int sign)
{
    int i;
    double sum;
    
    if(sign != 0)
    {
        for(i = 0; i < d; i++)
        {
            sum += x1[i] + x2[i];
        }
    }
    else
    {
        for(i = 0; i < d; i++)
        {
            sum += x1[i] -  x2[i];
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

double ** getVector(int N)
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

int getMax(int x, int y)
{
    if(x>y)
        return x;
    return y;
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

int searchIndex(int* arr, int N, int val)
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

int getD(FILE* input)
{
    int d, offset;
    char c;

    while((c = getc(input)) != '\n')
    {
        if(c == ',')
            d++;
        offset++;
    }
}