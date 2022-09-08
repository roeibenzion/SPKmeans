#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <math.h>
#include<string.h>
#include "spkmeans.h"

/*Kmeans*/
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
/*End Kmeans*/
/*Full SPK*/
PyObject * full_spk(double** mat, int N, int d, int K)
{
    Py_ssize_t p;
    Py_ssize_t q;
    int i;
    double **W, **D, **Lnorm, **temp, *eigenValues, **V, **U;
    PyObject* ret, *tempP;
        W = wamF(mat, N, d);
        D = ddgF(W, N);
        Lnorm = lnormF(W, D, N);
        temp = JacobiF(Lnorm, N);
        eigenValues = getVector(N);
        copyRows(eigenValues, temp[0], N);
        V = getMatrix(N, N);
        for(i = 0; i < N; i++)
        {
            copyRows(V[i], temp[i+1], N);
        }
        sortMatrixColumns(V, N, eigenValues);
        qsort(eigenValues,N,sizeof(double),compare);
        if(K < 1)
        {
            K = eigenGap(eigenValues, N);
        }
        U = obtainLargestK(V, N, K);
        formTfromU(U, N, K);
        ret = PyList_New(N);
        for(q = 0; q < N; q++)
        {
            tempP = PyList_New(K);
            for(p = 0; p < K; p++)
            {
                PyList_SetItem(tempP, p, PyFloat_FromDouble(U[q][p]));
            }
            PyList_SetItem(ret, q, tempP);
        }
        freeMatrix(U, N);
        freeVector(eigenValues);
        return ret;
}
static PyObject* fit(PyObject* self, PyObject* args)
{
    double ** datapoints, **centeroids, **ret;
    PyObject *goalP,* datapointsP, *centeroidsP, *centeroids_indexesP, *returnedMat, *temp;
    const char *kmeans, *jacobi;
    char* goal; 
    int N, d, K, t;
    Py_ssize_t i; 
    Py_ssize_t j;
    Py_ssize_t p;
    Py_ssize_t q;
    
    jacobi = "jacobi";
    kmeans = "kmeans";

    if(!PyArg_ParseTuple(args, "OOOOiii:fit", &datapointsP, &centeroidsP,&centeroids_indexesP, &goalP, &N, &d, &K))
    {
        return NULL;
    }

    datapoints = getMatrix(N,d);
    if (datapoints == NULL)
    {
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
    t = N+1;
    goal = strtok(PyBytes_AS_STRING(PyUnicode_AsEncodedString(PyObject_Repr(goalP), "utf-8", "~E~")), "'");
    if(!strcmp(goal, kmeans))
    {
        centeroids = getMatrix(K, d);
        if (centeroids == NULL)
        {
            printf("An error has occurred");
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
                puts("An error has occurred");
                } 
            }
        }
        return opt(K, datapoints, centeroids, 300, N, d, 0);
    }
    else if(!strcmp(goal, jacobi))
    {
        ret = getMatrix(t, N);
        navigator(goal,datapoints,N, d, ret);
        returnedMat = PyList_New(t);
        for(q = 0; q < N+1; q++)
        {
            temp = PyList_New(N);
            for(p = 0; p < N; p++)
            {
                PyList_SetItem(temp, p, PyFloat_FromDouble(ret[q][p]));
            }
            PyList_SetItem(returnedMat, q, temp);
        }
        freeMatrix(ret, N+1);
    }
    else if(!strcmp(goal, "spk"))
    {
       return full_spk(datapoints, N, d, K);
    }
    else
    {
         ret = getMatrix(N, N);
       navigator(goal,datapoints,N, d, ret);
       returnedMat = PyList_New(N);
        for(q = 0; q < N; q++)
        {
            temp = PyList_New(N);
            for(p = 0; p < N; p++)
            {
                PyList_SetItem(temp, p, PyFloat_FromDouble(ret[q][p]));
            }
            PyList_SetItem(returnedMat, q, temp);
        }
        freeMatrix(datapoints, N);
        freeMatrix(ret,N);
    }
    return returnedMat;
}
/*Cython API*/
static PyMethodDef capiMethods[] = {
    {"fit",
    (PyCFunction) fit,
    METH_VARARGS,
    PyDoc_STR("myspkmeans")},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef myspkmeans = { 
    PyModuleDef_HEAD_INIT, 
   "myspkmeans",
   NULL,
   -1,
   capiMethods 
};

PyMODINIT_FUNC
PyInit_myspkmeans(void)
{
    PyObject* m;
    m = PyModule_Create(&myspkmeans);
    if(!m)
    {
        return NULL;
    }
    return m;
}
