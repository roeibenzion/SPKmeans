//
// Created by saar barak on 9/8/2022.
//


def kMeanspp(df: np.ndarray, indexes: np.ndarray, N: int, k: int, d: int, max_iter: int):
        np.random.seed(0)
        rand = np.random.choice(N)
        newMu = df[rand]
        mu = [newMu]
        choices = [rand]
        D = [0 for i in range(N)]
        P = [0 for i in range(N)]
        i = 1
        while i < k:
                for l in range(N):
                        minDist = float("inf")  # find min in array algorithm
                        for j in range(i):
                                curDist = dist(mtx[l], mu[j])
                                if curDist < minDist:
                                        minDist = curDist
                        D[l] = minDist
                sumD = sum(D)
                for l in range(N):
                        P[l] = D[l] / sumD
                randI = np.random.choice(a=range(N), p=P)
                newMu = df[randI]
                choices.append(randI)
                mu.append(newMu)
                i += 1
        print(str(choices)[1:-1].replace(" ", ""))  # print the chosen indexes

        # mu = [list(float(mu[i][j]) for j in range(k)) for i in range(k)]
        c = [[0 for x in range(d)] for y in range(k)]
        centeroids = np.array(choices)
        for row in range(len(centeroids)):
                for col in range(len(centeroids[row])):
                        c[row][col] = centeroids[row][col]  ###

        datapoints = [[0 for x in range(d)] for y in range(N)]
        for row in range(len(df)):
                for col in range(len(df[row])):
                        datapoints[row][col] = df[row][col]

        centeroids = spk.fit(datapoints,c, mu,"kmeans",N,d,k)
        arr = np.array(centeroids)
        print_matrix(np.reshape(arr, (k,d)))



def dist(x: list, u: list) -> float:
        return (sum((x[i] - u[i]) ** 2 for i in range(len(x))))