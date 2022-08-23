#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

static PyObject* fit(PyObject* self, PyObject* args);
PyObject* opt(int, double**, double**, int, int, int, double);
void init_clusters(double ** clusters,int K, int d, int* sizes);
double* update_centeroid(double *cluster, int cluster_size, int d);
int find_min(double* dataPoint, double** centeroids, int K, int d);
double calculate_dist(double* point1, double* point2, int d);
void sum_vectors(double* vect1, double* vect2, int d);
int decimal(char* c);
double ** getMatrix(int N, int d);
void freeMatrix(double** m, int k);

static PyObject* fit(PyObject* self, PyObject* args)
{
    Py_ssize_t i; 
    Py_ssize_t j;
    Py_ssize_t p;
    int N; 
    int d;
    int K;
    int max_iter; 
    double epsilon;
    PyObject *centeroids_indexesP;
    PyObject *centeroidsP; 
    PyObject *datapointsP;
    double **centeroids; 
    double **datapoints; 
    
    if(!PyArg_ParseTuple(args, "iiiidOOO:fit", &N, &d, &K, &max_iter,&epsilon, &centeroids_indexesP, &centeroidsP, &datapointsP))
    {
        return NULL;
    }

    datapoints = getMatrix(N,d);
    if (datapoints == NULL){
        printf("An Error Has Occurred");
        return Py_BuildValue("");
    }
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < d; j++)
        {
            datapoints[i][j] = PyFloat_AS_DOUBLE(PyList_GetItem(PyList_GetItem(datapointsP, i),j));
            if (PyErr_Occurred())
            {
            puts("An Error Has Occurred");
            } 
        }
    }
    centeroids = getMatrix(K, d);
    if (centeroids == NULL){
        printf("An Error Has Occurred");
        freeMatrix(centeroids ,K);
        return Py_BuildValue("");
    }
    for (i = 0; i < K; i++)
    {
        p = PyLong_AS_LONG(PyList_GetItem(centeroids_indexesP, i));
        for (j = 0; j < d; j++)
        {
            centeroids[i][j] = PyFloat_AS_DOUBLE(PyList_GetItem(PyList_GetItem(datapointsP, p),j));
            if (PyErr_Occurred())
            {
            puts("An Error Has Occurred");
            } 
        }
    }
    return opt(K, datapoints, centeroids, max_iter, N, d, epsilon);
}

static PyMethodDef capiMethods[] = {
    {"fit",
    (PyCFunction) fit,
    METH_VARARGS,
    PyDoc_STR("mykmeanssp")},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef mykmeanssp = { 
    PyModuleDef_HEAD_INIT, 
   "mykmeanssp",
   NULL,
   -1,
   capiMethods 
};

PyMODINIT_FUNC
PyInit_mykmeanssp(void)
{
    PyObject* kmeans;
    kmeans = PyModule_Create(&mykmeanssp);
    if(!kmeans)
    {
        return NULL;
    }
    return kmeans;
}

PyObject* opt(int K, double** dataPoints, double** centeroids, int max_iter, int N, int d, double epsilon)
{
    Py_ssize_t q, p;
    PyObject * returnedCenteroids;
    PyObject* temp;
    int reps, done, i, min_cluster, j;
    double* new_centeroid;
    double **clusters;
    int *cluster_sizes;
    new_centeroid = (double*) calloc(d, sizeof(double));

    clusters = getMatrix(K, d);
    cluster_sizes = (int*) calloc(K, sizeof(int));
    reps = 0;
    done = 0;

    while(done == 0 && reps < max_iter)
    {
        init_clusters(clusters, K, d, cluster_sizes);
        done = 1;
        for(i = 0; i < N; i++)
        {
            min_cluster = find_min(dataPoints[i], centeroids, K, d);
            sum_vectors(clusters[min_cluster], dataPoints[i], d);
            cluster_sizes[min_cluster]++;
        }
        for (i = 0; i < K; i++)
        {
            if(cluster_sizes[i] > 0)
                new_centeroid = update_centeroid(clusters[i], cluster_sizes[i], d);

            if(calculate_dist(new_centeroid, centeroids[i], d) >= epsilon)
                done = 0;
            for(j = 0; j < d; j++)
            {
                centeroids[i][j] = new_centeroid[j];
            }
        }
        reps++;
    }

    returnedCenteroids = PyList_New(K);
    for(q = 0; q < K; q++)
    {
        temp = PyList_New(d);
        for(p = 0; p < d; p++)
        {
            PyList_SetItem(temp, p, PyFloat_FromDouble(centeroids[q][p]));
        }
        PyList_SetItem(returnedCenteroids, q, temp);
    }
    freeMatrix(centeroids, K);
    free(new_centeroid);
    freeMatrix(dataPoints,K);
    free(cluster_sizes);
    return returnedCenteroids;
}
double* update_centeroid(double *cluster, int cluster_size, int d)
{
    int i;
    for (i = 0; i < d; i++)
    {
        cluster[i] = (double)(cluster[i] / cluster_size);
    }
    return cluster;
}

int find_min(double* dataPoint, double** centeroids, int K, int d)
{
    int index_min = 0, i;
    double min, dist;
    min = calculate_dist(centeroids[0], dataPoint, d);
    for(i = 1; i < K; i++)
    {
        dist = calculate_dist(centeroids[i], dataPoint, d);
        if(dist < min)
        {
            index_min = i;
            min = dist;
        }
    }
    return index_min;
}
double calculate_dist(double* point1, double* point2, int d)
{
    int i;
    double sum;
    sum = 0.0;
    for(i = 0; i < d; i++)
    {
        sum += pow((point1[i] - point2[i]) , 2);
    }
    sum = pow(sum, 0.5);
    return sum;
}

void sum_vectors(double* vect1, double* vect2, int d)
{
    int i;
    for(i = 0; i < d; i++)
    {
        vect1[i] += vect2[i];
    }
}

void init_clusters(double ** clusters,int K, int d, int* sizes)
{
    int i,j;
    for(i = 0; i < K; i++)
    {
        sizes[i] = 0;
        for(j = 0; j < d; j++)
        {
            clusters[i][j] = 0.0;
        }
    }
}

int decimal(char* c)
{
    int i = 0;
    if(c[0] == '-')
        i++;
    for(; c[i] != 0; i++)
    {
        if(c[i] > '9' || c[i] < '0')
            return 0;
    }
    return 1;
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

