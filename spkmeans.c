#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


void navigator(char* goal, double *vect, int N, int d)
{
    double **W, **D, **Lnorm;
    char* spk, *wam, *ddg, *lnorm, *jacobi;

    spk = "spk";
    wam = "wam";
    ddg = "ddg";
    lnorm = "lnorm";
    jacobi = "jacobi";

    if(strcmp(goal,wam))
    {
        W = wamF(vect, N, D);
    }

    else if(strcmp(goal,ddg))
    {
        W = wamF(vect, N, D);
        D = ddgF(W, N);
    }

    else if(strcmp(goal,lnorm))
    {
        W = wamF(vect, N, D);
        D = ddgF(W, N);
        Lnorm = lnormF(W, D, N);
    }

    else if(strcmp(goal,jacobi))
    {
        W = wamF(vect, N, D);
        D = ddgF(W, N);
        Lnorm = lnormF(W, D, N);
        //whats the best way to do that?
    }
    else if(strcmp(goal, spk))
    {
        W = wamF(vect, N, D);
        D = ddgF(W, N);
        Lnorm = lnormF(W, D, N);
        //whats the best way to do that?
    }
    else
    {
        printf("Invalid Input!");
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
                c = sumVectors(vect[i], vect[j], d);
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

double** JacobiF(double** A, int N)
{
    int i,j, iMax, jMax, countIter, firstIter;
    double theta, t, c, s, maxVal, offSquareA;
    maxVal = 0.0;
    double ** P, **V;

    countIter = 0;
    firstIter = 0;
    P = getMatrix(N, N);
    for(i = 0; i < N; i++)
        {
            P[i][i] = 1;
        }

    while(countIter < 100)
    {
        countIter++;
        for(i = 0; i < N; i++)
        {
            for(j = 0; j < N; j++)
            {
                maxVal = getMax(maxVal, A[i][j]); 
                iMax = i;
                jMax = j;
            }
        }
        theta = (A[jMax][jMax] - A[iMax][iMax]) / 2*(A[iMax][jMax]);
        t = findT(theta);
        c = 1/sqrt(pow(t,2) + 1);
        s = t*c;
        P[iMax][jMax] = s;
        P[jMax][iMax] = -s;
        P[i][i] = c;
        P[j][j] = c;
        if(firstIter == 0)
        {
            firstIter = 1;
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
    return V;
}

int eigenGap(double** V, int N)
{
    double deltaI, temp;
    int i, iter, maxIndex;

    iter = N/2;
    if(iter == 1)
    {
        return 1;
    }
    i = 1;
    deltaI = abs(V[i][i] - V[i+1][i+1]);
    maxIndex = 1;
    for(i = 2; i < iter; i++)
    {
        temp = abs(V[i][i] - V[i+1][i+1]);
        if(deltaI < temp)
        {
            maxIndex = i;
            deltaI = temp;
        }
    }
    return maxIndex;
}

double** sortMatrixColumns(double** V, int N, double** A)
{
    int i, j;
    double* copy, **newV;

    copy = getVector(N);
    newV = getMatrix(N, N);
    for(i = 0; i < N; i++)
    {
        copy[i] = A[i][i];
    }
    qsort(copy,N,sizeof(double),compare);

    for(i = 0; i < N; i++)
    {
        /////here
    }
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
            ATag[r][j] = c*A[r][j] + s*A[r][i];
        }
    }
    ATag[i][i] = pow(c,2)*(A[i][i]) + pow(s,2)*(A[j][j]) - 2*s*c*A[i][j];
    ATag[j][j] = pow(c,2)*(A[i][i]) + pow(c,2)*(A[j][j]) + 2*s*c*A[i][j];
    ATag[i][j] = 0;
}

int checkConverstion(double** ATag, double offSquareA, int N, int i, int j)
{
    int r;
    double offSquareATag, eps;

    eps = pow(10,-5);
    offSquare(ATag, N, i, j);
    return (offSquareA - offSquareATag) <= eps;
}

double offSquare(double **A, int N, int i, int j)
{
    int r;
    double sum;

    for(r = 0; r < N; r++)
    {
        sum += pow(A[r][i],2) + pow(A[r][j],2);
    }
    return sum;
}
double sumVectors(double* x1, double* x2, int d)
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
    int i,j;
    double **C;
    C = getMatrix(N,N);
    for(i = 0; i < N; i++)
    {
        for(j = 0; j < N; j++)
        {
            C[i][j] = A[i][j] * B[i][j];
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