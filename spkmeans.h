#ifndef SPKMEANS_H_
#define SPKMEANS_H_

void generateId(double** A, int N);
double sumRow(double* row, int N, int squared);
void matToThePow(double** M, int N, double exponent, int diag);
double** matMul(double **A, double **B, int N);
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
void navigator(char* goal, double** mat, int N, int d, double** ret, int K);
void copyCol(double** dst, double** src, int N, int dstIndex, int srcIndex);
int compare( const void* a, const void* b);

#endif 