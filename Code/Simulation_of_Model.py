from numpy import *
from math import *
import csv

xPrecision = 0.01
dPrecision = 0.1
lPrecision = 0.1
pPrecision = 0.05
wPrecision = 0.1
wdRatio = 0.9

def getY(d, l, w, x):
        return (l * d * w + l * x * sqrt(x * x + d * d - w * w)) / (x * x + d * d) - x

def G(d, l, p, w):
    realP = sqrt(w * w + l * l - d * d)
    if p > realP:
        p = realP
    x = 0
    result = -1
    while x < p + xPrecision:
        y = getY(d, l, w, x)
        if y > result:
            result = y
        x += xPrecision
    return result



# getDataAndResults

def getDataAndResults():

    dataFile = open("Data & Results.csv", 'w')
    datas = []
    results = []

    d = 0.8
    while d < 1.6:
        l = d + 0.1
        while l < 2.1:
            w = 0.1
            while w < d:
                p = 0.25
                pMax = sqrt(w * w + l * l - d * d)
                while p < pMax:


                    result = G(d, l, p, w)
                    data = [d, l, p, w, 1]
                    datas.append(data)
                    results.append(result)

                    p += pPrecision
                w += wPrecision
            l += lPrecision
        d += dPrecision

    for i in range(len(datas)):
        for j in range(5):
            dataFile.write(str(datas[i][j]))
            dataFile.write(", ")
        dataFile.write(str(results[i]))
        dataFile.write("\n")
    dataFile.close()



# leastSquareSolve

def leastSquareSolve(dataName, residualName):
    dataFile = open(dataName, encoding = "UTF-8")
    dataReader = csv.reader(dataFile)
    datas = array(list(dataReader), dtype = float)
    matrixShape = datas.shape
    AArray, bArray = hsplit(datas, [matrixShape[1] - 1])
    A = mat(AArray)
    b = mat(bArray)
    ATA = matmul(A.T, A)
    ATb = matmul(A.T, b)
    x = linalg.solve(ATA, ATb)
    for i in range(matrixShape[1] - 1):
        x[i, 0] = int(x[i, 0] * 1000) / 1000
    print(x)
    print("\n")

    residualCount = [0, 0, 0, 0] # r>=10, 5<=r<10, 1<=r<5, 0<=r<1
    residualFile = open(residualName, "w", encoding='UTF-8')
    residualVector = b - matmul(A, x)
    for i in range(matrixShape[0]):
        for j in range(matrixShape[1]):
            residualFile.write(str(datas[i, j]) + ", ")
        residualFile.write(str(residualVector[i, 0]) + "\n")
        if abs(residualVector[i, 0]) >= 0.1:
            residualCount[0] += 1
        elif abs(residualVector[i, 0]) >= 0.05:
            residualCount[1] += 1
        elif abs(residualVector[i, 0]) >= 0.01:
            residualCount[2] += 1
        else:
            residualCount[3] += 1

    print(str(matrixShape[0]))
    for i in range(4):
        print(str(residualCount[i]) + ", " + str(residualCount[i] * 100 / matrixShape[0]), "%")
    residualFile.close()



# getQuadraticData

def getQuadraticData():
    dataFile = open("Data & Results.csv", encoding = "UTF-8")
    dataReader = csv.reader(dataFile)
    datas = array(list(dataReader), dtype = float)
    dataFile.close()

    quadraticDatas = []
    for data in datas:
        d = data[0]
        l = data[1]
        p = data[2]
        w = data[3]
        result = data[5]
        quadraticData = [d * d, l * l, p * p, w * w, d * l, d * p, d * w, l * p, l * w, p * w, d, l, p, w, 1, result]
        quadraticDatas.append(quadraticData)

    quadraticDataFile = open("Quadratic Simulation Data.csv", 'w')
    for i in range(len(quadraticDatas)):
        for j in range(15):
            quadraticDataFile.write(str(quadraticDatas[i][j]))
            quadraticDataFile.write(", ")
        quadraticDataFile.write(str(quadraticDatas[i][15]))
        quadraticDataFile.write("\n")
    quadraticDataFile.close()



# compareSmallWResults

def compareSmallWResults():
    
    w = 0.1
    wresults = []

    d = 0.8
    while d < 1.6:
        l = d + 0.1
        while l < 2.1:
            
            p = 0.25
            pMax = sqrt(w * w + l * l - d * d)
            while p < pMax:

                wresults.append(G(d, l, p, 0) - G(d, l, p, w))

                p += pPrecision
            l += lPrecision
        d += dPrecision

    wCompareFile = open("W Zero Approximation Error.csv", 'w')
    for i in range(len(wresults)):
        wCompareFile.write(str(wresults[i]) + "\n")
    wCompareFile.close()



# getSepaData

def getSepaData():

    dataFile = open("Normal Data.csv", 'w')
    specialFile = open("Special Data.csv", 'w')
    datas = []
    results = []

    d = 0.8
    while d < 1.6:
        l = d + 0.1
        while l < 2.1:
            w = 0.1
            while w < l and w < d:
                p = 0.25
                pMax = sqrt(w * w + l * l - d * d)
                while p < pMax:
                    result = G(d, l, p, w)
                    if w / d >= 0.95:
                        specialFile.write(str(d * d) + ", " + str(l * l) + ", " + str(p * p) + ", " + str(w * w) + ", " + str(d * l) + ", " + str(d * p) + ", " + str(d * w) + ", " + str(l * p) + ", " + str(l * w) + ", " + str(p * w) + ", " + str(d) + ", " + str(l) + ", " + str(p) + ", " + str(w) + ", " + "1" + ", " + str(result) + "\n")
                    else:
                        data = [d * d, l * l, p * p, w * w, d * l, d * p, d * w, l * p, l * w, p * w, d, l, p, w, 1]
                        datas.append(data)
                        results.append(result)

                    p += pPrecision
                w += wPrecision
            l += lPrecision
        d += dPrecision

    for i in range(len(datas)):
        for j in range(15):
            dataFile.write(str(datas[i][j]))
            dataFile.write(", ")
        dataFile.write(str(results[i]))
        dataFile.write("\n")
    dataFile.close()



# specialLeastSquareSolve

def specialLeastSquareSolve(normalDataName, specialDataName, residualName):
    dataFile = open(normalDataName, encoding = "UTF-8")
    dataReader = csv.reader(dataFile)
    datas = array(list(dataReader), dtype = float)
    print(datas)
    matrixShape = datas.shape
    AArray, bArray = hsplit(datas, [matrixShape[1] - 1])
    A = mat(AArray)
    b = mat(bArray)
    ATA = matmul(A.T, A)
    print(ATA)
    ATb = matmul(A.T, b)
    x = linalg.solve(ATA, ATb)
    for i in range(matrixShape[1] - 1):
        x[i, 0] = int(x[i, 0] * 1000) / 1000
    x[matrixShape[1] - 2, 0] = 0
    print(x)
    print("\n")

    specialFile = open(specialDataName, encoding = 'UTF-8')
    specialReader = csv.reader(specialFile)
    specialDatas = array(list(specialReader), dtype = float)

    residualCount = [0, 0, 0, 0]
    residualFile = open(residualName, "w", encoding='UTF-8')
    residualVector = b - matmul(A, x)
    for i in range(matrixShape[0]):
        for j in range(matrixShape[1]):
            residualFile.write(str(datas[i, j]) + ", ")
        residualFile.write(str(residualVector[i, 0]) + "\n")
        if abs(residualVector[i, 0]) >= 0.1:
            residualCount[0] += 1
        elif abs(residualVector[i, 0]) >= 0.05:
            residualCount[1] += 1
        elif abs(residualVector[i, 0]) >= 0.01:
            residualCount[2] += 1
        else:
            residualCount[3] += 1
    for i in range(specialDatas.shape[0]):
        for j in range(matrixShape[1]):
            residualFile.write(str(specialDatas[i, j]) + ", ")
        residualFile.write(str(specialDatas[i, 15] - specialDatas[i, 11]) + "\n")
        if abs(specialDatas[i, 15] - specialDatas[i, 11]) >= 0.1:
            residualCount[0] += 1
        elif abs(specialDatas[i, 15] - specialDatas[i, 11]) >= 0.05:
            residualCount[1] += 1
        elif abs(specialDatas[i, 15] - specialDatas[i, 11]) >= 0.01:
            residualCount[2] += 1
        else:
            residualCount[3] += 1
    residualFile.close()

    dataNum = residualCount[0] + residualCount[1] + residualCount[2] + residualCount[3]
    for i in range(4):
        print(str(residualCount[i]) + ", " + str(residualCount[i] * 100 / dataNum), "%")
    residualFile.close()



instruction = int(input("WHAT: "))
if instruction == 1:
    getDataAndResults()
elif instruction == 2:
    leastSquareSolve("Data & Results.csv", "Linear Simulation Residual.csv")
elif instruction == 3:
    getQuadraticData()
elif instruction == 4:
    leastSquareSolve("Quadratic Simulation Data.csv", "Quadratic Simulation Residual.csv")
elif instruction == 5:
    compareSmallWResults()
elif instruction == 6:
    getSepaData()
elif instruction == 7:
    specialLeastSquareSolve("Normal Data.csv", "Special Data.csv", "Case Simulation Residual.csv")