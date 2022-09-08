import math
import sys
import numpy as np
import pandas as pd
import myspkmeans as spk

###OLD KMEANS HW2###
def getMinDl(xl: np.ndarray, centeroids: np.ndarray, i:int) -> int:
    min_d = np.sum(np.power((xl-centeroids[0]),2)) 
    for y in centeroids:
        min_d = min((min_d, np.sum(np.power((xl-y),2))))
        
    return min_d

def getProb(D: list) -> np.ndarray:
    P = np.array([])
    s = sum(D)
    for l in range(0,len(D)):
        P = np.append(P, (D[l]/s))
    return P

def kMeanspp(df:np.ndarray, indexes:np.ndarray, N: int, k: int, d: int, max_iter: int): 
    centeroids_indexes = []
    np.random.seed(0)

    centeroids_indexes.append(np.random.choice(a = np.arange(N), size = 1))
    centeroids = np.array(df[centeroids_indexes[0]])
    for i in range(1,k):
        D = []
        for l in range (0,N):
            D.append(getMinDl(df[l], centeroids, i))
        P = getProb(D)
        centeroids_indexes.append(np.random.choice(a = np.arange(N) ,size = 1, p = P))
        centeroids = np.append(centeroids, df[centeroids_indexes[i]], axis=0)
    c = [[0 for x in range(d)] for y in range(k)]
    centeroids = np.array(centeroids)
    for row in range(len(centeroids)):
        for col in range(len(centeroids[row])):
            c[row][col] = centeroids[row][col] 

    datapoints = [[0 for x in range(d)] for y in range(N)]
    for row in range(len(df)):
        for col in range(len(df[row])):
            datapoints[row][col] = df[row][col]
    
    for x in range(len(centeroids_indexes)):
        if(x < len(centeroids_indexes)-1):
                print(indexes[centeroids_indexes[x][0]] ,end=',')
        else:
                print(indexes[centeroids_indexes[x][0]])
        centeroids_indexes[x] = centeroids_indexes[x][0]
    centeroids = spk.fit(datapoints,c, centeroids_indexes,"kmeans",N,d,k)
    arr = np.array(centeroids)
    np.reshape(arr, (k,d))
    for x in range(len(arr)):
        for y in range(len(arr[x])):
            if(y < len(arr[x])-1):
                print('%.4f' % arr[x][y] ,end=',')
            else:
                print('%.4f' % arr[x][y])

###PRINT###
def print_matrix(arr):
      for x in range(len(arr)):
        for y in range(len(arr[x])):
            if(y < len(arr[x])-1):
                print('%.4f' % arr[x][y] ,end=',')
            else:
                print('%.4f' % arr[x][y])

###NEW###
k, goal, file = (int(0), '', '')
input = sys.argv

k = int(input[1])
goal = str(input[2])
file = str(input[3])

if(len(sys.argv) != 4):
   print("invalid input!")
   exit(1)


try:
    df = pd.read_csv(file, header=None)
    N = len(df)
    d = len(df.columns)
    df = df.to_numpy()
    df = df.astype(float)
    if(k >= N or k < 0 or goal not in ["spk", "wam", "ddg", "lnorm", "jacobi"]):
        print("invalid input!")
        exit(1)
    indexes = list(range(0, N))
    c = [[0 for x in range(d)] for y in range(N)]
    for row in range(N):
        for col in range(d):
            c[row][col] = df[row][col]
    if(goal == "spk"):
        T_matrix = spk.fit(c, None, None, goal,N, d, k)
        T_matrix = np.array(T_matrix)
        k = np.atleast_2d(T_matrix).shape[1]
        kMeanspp(T_matrix, indexes, N, k, k, 300)
    else:
        matrix = spk.fit(c, None, None, goal,N, d, k)
        print_matrix(matrix)
    
except:
    print("An error has occuredPy")
    exit(1)


    