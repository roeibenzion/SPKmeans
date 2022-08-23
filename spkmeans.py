import imp
import math
import sys
import numpy as np
import pandas as pd

def wan(X: np.array, N: int, d: int) -> np.array:
    return
    #send to C

def ddg(X: np.array, N: int, d: int) -> np.array:
    return 
    #send to C
def Lnorm(X: np.array, N: int, d: int) -> np.array:
    return 
    #send to C
def Jacobi(M: np.array, N: int) -> np.array:
    #gets a matrix
    return
    #send to C
def SPK(X: np.array, N: int, d: int):
    #send to C
    #send the returned value to Kmeans
    return

k, goal, file = (int(0), '', '')
input = sys.argv

if(len(sys.argv) != 3):
   print("invalid input!")
   exit(1)

df = pd.read_csv(file, header=None)
df = df.to_numpy()
N = len(df)
d = len(df[1])

if(k >= N):
    print("invalid input!")
    exit(1)
# if(k == 0) - notify C
#notify C what's the goal 
k = input[1]
goal = input[2]
file = input[3]

df = pd.read_csv(file, header=None)
df = df.to_numpy()
N = len(df)
d = len(df[1])
if(k > N):
    print("invalid input!")
    exit(1)

if(goal == "wan"):
    W = wan(df, N, d)
elif(goal == "ddg"):
    D = ddg(df, N, d)
elif(goal == "Lnorm"):
    L = ddg(df, N, d)
elif(goal == "Jacobi"):
    V = Jacobi(df, N)
elif(goal == "SPK"):
    SPK(df, N, d)


