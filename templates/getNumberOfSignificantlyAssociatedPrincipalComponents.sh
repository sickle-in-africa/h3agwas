#!/usr/bin/env bash

# This script needs to be cleaned. It works but is very hard to read.
# It basically reads the file ${pvalueTable} which stores principe 
# components and their p values. The first column is the ide of each
# principal component, and the second colmn is the p value of each
# principal component. The script takes in a p value threshold:
#
#   params.ancestry.maxPrincipalComponentAssociationPvalue
#
# and return the least significant principal component that has a 
# p value below this threshold.

awk -F ',' '\$2 <= ${params.ancestry.maxPrincipalComponentAssociationPvalue} {print \$1}' ${pvalueTable} > temp.dat
numberOfPrincipalComponents=`awk -F ',' 'BEGIN{a=   0}{if (\$1>0+a) a=\$1} END{print a}' temp.dat`