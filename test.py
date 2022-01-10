import os
import sys
import subprocess
from decimal import *
getcontext().prec = 8



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
                data[-1].append(float(x))
                
    return data


def checkBounds(out):
    stdOut = out.stdout.decode("utf-8").split(",")
    width = int(stdOut[1])
    height = int(stdOut[2])
    iterations = int(stdOut[4])
    gridString = stdOut[6]
    
    grid = getGrid(gridString)
    
    for x in range(width):
        for y in range(height):
            #print( grid[y][x], x, y, end=" ")
            if(x == 0 or (y == 0 and x != width -1)):
                if(abs(grid[y][x] - 1.0) > 0.01):
                    return {"success": False, "reason": "Boundery value changed: 1 -> {} at x={}, y={}".format(grid[y][x], x, y), "data": grid}
            if(x == width - 1 or (y == height - 1 and x != 0)):
                if(abs(grid[y][x] - (-1.0)) > 0.01):
                    return {"success": False, "reason": "Boundery value changed: -1 -> {} at x={}, y={}".format(grid[y][x], x, y), "data": grid}
        
        #print("")


    return {"success": True, "reason": "", "data": grid, "iterations": iterations}


def compareOut(grid0, grid1, precision):
    totalDifference = 0
    for i in range(len(grid0)):
        for j in range(len(grid0[i])):
            cell0 = grid0[i][j]
            cell1 = grid1[i][j]
            totalDifference += abs(cell0 - cell1)

    if(totalDifference > precision):
        print("Total diff", totalDifference)
        return False
    return True


def validate():
    width = height = 256
    precision = 0.001

    #compile sequential version
    compileExitCode = os.system("gcc ./seq_main.c -o ./bin/seq_main")
    if(compileExitCode != 0):
        print("Could not compiler sequential version")
        return
    
    # Get sequential data
    seqDataString = subprocess.run("./bin/seq_main {} {} {}".format(width, height, precision), shell=True, capture_output=True)
    allSeqData = checkBounds(seqDataString)
    sequentialData = allSeqData["data"]
    sequentialIterations = allSeqData["iterations"]


    for workers in [2, 4, 8, 16, 32, 64, 128]:
        command = "mpirun -np {} --oversubscribe ./bin/main {} {} {}".format(
            workers, width, height, precision)
        checkNum  = 0
        while(checkNum < 1):
            processOutput = subprocess.run(
                command, shell=True, capture_output=True)
            #Check the bounderies 
            check = checkBounds(processOutput)
            if(check["success"] == True):
                #Check against sequential version
                if(not compareOut(sequentialData, check["data"], 0.00001) or abs(check["iterations"] - sequentialIterations) != 0):
                    check["success"] = False
                    check["reason"] = "Sequential and parallel programs output differ by too much"
                    break
            else:
                break
            checkNum += 1



        print("Workers: {}, Dim: {}x{}, Success: {}".format(workers,
                                                            width, height, check["success"]), end="")
        if(check["success"] == False):
            print(", Reason: {}".format(check["reason"]))
        else:
            print("")


def main(mode):
    compileExitCode = os.system("./scripts/compile.sh")
    if(compileExitCode != 0):
        print("Failed to compile")
        return
    if(mode == "performance"):
        #performance()
        pass
    elif(mode == "validate"):
        validate()


if(__name__ == "__main__"):
    main(sys.argv[1])
