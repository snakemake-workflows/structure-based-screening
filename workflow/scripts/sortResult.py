"""sorting and outputting the N>1 best results from a docking run"""
import subprocess
import gzip
import os
from math import ceil
inPath = snakemake.input[0]
outFile = snakemake.output[0]
num = float(snakemake.config["RESULT_NUMBER"])

count=0
lineN=0
ids=0
swapFlg=False
scoreList=[]
strucList=[]
tempList=[]

getLigNum = "zgrep -c 'REMARK RECEPTOR' " + inPath
ligand_num = int(subprocess.check_output(getLigNum, shell=True))

if num > 1:
    listSize = num
else:
    total=float(ligand_num)
    listSize = ceil(num*total)

with gzip.open(inPath, 'rt', encoding='utf-8') as inFile:
    for line in inFile:
        if "REMARK RECEPTOR" in line:
            lineN=0
            if count>0:
                if count<=listSize:
                    strucList.append(tempList)
                elif swapFlg:
                    strucList[ids]=tempList
                    swapFlg=False
            count=count+1
            tempList=[]
        if lineN==3:
            strs=line.split()
            curValue=float(strs[3])
            if count<=listSize:
                scoreList.append(curValue)
            else:
                maxValue=max(scoreList)
                if curValue<maxValue:
                    swapFlg=True
                    ids=scoreList.index(maxValue)
                    scoreList[ids]=curValue
        tempList.append(line)
        lineN=lineN+1
    scoreDict={}

    for index, item in enumerate(scoreList):
        scoreDict[index]=item

    sortList=sorted(scoreDict, key=scoreDict.get)
    with open(outFile, 'w', encoding='utf-8') as outF:
        for index in sortList:
            for line in strucList[index]:
                outF.write(line)
os.chmod(outFile, 0o400) ## TODO: remove after included in Snakemake
