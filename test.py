import os
import sys
import subprocess
from decimal import *
getcontext().prec = 8


def performance():
    os.system("echo workers,width,height,precision,iterations,time")
    repeat = 1
    for workers in [44]:
        for scale in [4096, 4096*2, 4096*4]:
            for r in range(repeat):
                width = height = scale
                precision = 0.001
                command = "./main {} {} {} {} 0".format(
                    workers, width, height, precision)
                subprocess.run(command, shell=True)


def getGrid(s):
    data = []
    # Build 2D array of all values from string
    for rowString in s.split("\n"):
        row = rowString.split(" ")
        if(len(row) == 1):  # Skip just spaces or empty index's
            continue
        # print(row)
        data.append([])
        for x in row:
            if(x != ""):
                data[-1].append(Decimal(x))
    return data


def checkBounds(out):
    stdOut = out.stdout.decode("utf-8").split(",")
    width = int(stdOut[1])
    height = int(stdOut[2])

    iterations = int(stdOut[4])

    endGridString = stdOut[6]
    startGridString = stdOut[7]

    endGrid = getGrid(endGridString)
    startGrid = getGrid(startGridString)

    for x in range(width):
        for y in range(height):

            endCell = endGrid[x][y]  # Current cell
            startCell = startGrid[x][y]
            # Program is set to have bounderies at 0 and 1 that are constant
            if(x == 0 or y == 0 or x == width-1 or y == height-1):
                if(endCell != startCell):
                    return {"success": False, "reason": "Boundery value changed: {}->{}".format(startCell, endCell), "data": endGrid}

    return {"success": True, "reason": "", "data": endGrid, "iterations": iterations}


def compareOut(grid0, grid1, precision):
    totalDifference = 0

    for i in range(len(grid0)):
        for j in range(len(grid0[i])):
            cell0 = grid0[i][j]
            cell1 = grid1[i][j]
            totalDifference += abs(cell0 - cell1)

    if(totalDifference > precision):
        print(totalDifference)
        return False
    return True


def validate():
    width = height = 256
    precision = 0.000001
    sequentialData = []
    sequentialIterations = 0

    for workers in [1, 2, 4, 8, 16, 32, 44]:
        command = "./main {} {} {} {} 1".format(
            workers, width, height, precision)
        checkNum  = 0
        while(checkNum < 50):
            processOutput = subprocess.run(
                command, shell=True, capture_output=True)
            #Check the bounderies 
            check = checkBounds(processOutput)
            if(check["success"] != False):
                #Set the results of the sequential test
                if(workers == 1):
                    sequentialData = check["data"]
                    sequentialIterations = check["iterations"]
                    break
                else:
                    #Check against sequential version
                    if(not compareOut(sequentialData, check["data"], 0.00001) and abs(check["iterations"] - sequentialIterations) != 0):
                        check["success"] = False
                        check["reason"] = "Sequential and parallel programs output differ by too much"
                        break
            else:
                break
            checkNum += 1



        print("Threads: {}, Dim: {}x{}, Success: {}".format(workers,
                                                            width, height, check["success"]), end="")
        if(check["success"] == False):
            print("Reason: {}, Data: ".format(check["reason"]))
            print(check["data"][:1][:8])
        else:
            print("")


def main(mode):
    compileExitCode = os.system("gcc main.c -o main -lpthread -O3")
    if(compileExitCode != 0):
        print("Failed to compile")
        return
    if(mode == "performance"):
        performance()
    elif(mode == "validate"):
        validate()


if(__name__ == "__main__"):
    main(sys.argv[1])
