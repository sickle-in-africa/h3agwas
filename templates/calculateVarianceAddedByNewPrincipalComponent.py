#!/usr/bin/env python
outputStream = open("pvalueAndVarianceAndVarianceAddedTable.csv", "w")
print >> outputStream, "#PC,pvalue,r-squared,r-squared_added"
with open("${pvalueAndVarianceTable}") as file:
    rsquared = 0.0
    for line in file:
        if line.startswith('#'):
            continue
        linearray = line.rstrip().split(",")
        linearray.append(float(linearray[2])-rsquared)
        print >> outputStream, ",".join(str(e) for e in linearray)
        rsquared = float(linearray[2])
outputStream.close()
