import math
import sys
import numpy as np
import pandas as pd
import mykmeanssp as km

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
    
    centeroids = km.fit(N, d, k, max_iter,eps, centeroids_indexes, c, datapoints)
    arr = np.array(centeroids)
    np.reshape(arr, (k,d))
    for x in range(len(arr)):
        for y in range(len(arr[x])):
            if(y < len(arr[x])-1):
                print('%.4f' % arr[x][y] ,end=',')
            else:
                print('%.4f' % arr[x][y])

k, max_iter, eps, file1, file2 = (int(0), int(300), float(0), '', '')
input = sys.argv

if(len(sys.argv) == 6):
    k = int(input[1])
    max_iter = int(input[2])
    eps = float(input[3])
    file1 = input[4]
    file2 = input[5]
elif(len(sys.argv) == 5):
    k = int(input[1])
    eps = float(input[2])
    file1 = input[3]
    file2 = input[4]
else:
   print("invalid input!")
   exit(1)

if(int(k) < 1 or int(max_iter) < 1):
    print("invalid input!")
    exit(1)
try:
    df1 = pd.read_csv(file1, header=None)
    df2 = pd.read_csv(file2, header=None)
    df = pd.merge(df1, df2, on=0, how='inner')
    df = df.sort_values(0)
    indexes = np.copy(df[0])
    indexes = indexes.astype(int)
    df = df.drop(columns=0)
    df = df.to_numpy()
    N = len(df)
    d = len(df[1])
    if(k > N):
        print("invalid input!")
        exit(1)
    kMeanspp(df, indexes, N, k, d, max_iter)
except:
    print("An error has occured")
    exit(1)



    