__author__ = 'yzhu14'
#March 2017
from albinoBasics import *

outdir = 'output'
nNew=16
start ="0"
whichFileN = sys.argv[1]
whichFileP = sys.argv[2]
nCommentLines = 2
note=""

########################################################################
print(nNew,"energy bins")
ntestpointsN = sum(1 for line in open(whichFileN))
ntestpointsP = sum(1 for line in open(whichFileP))
assert ntestpointsN == ntestpointsP
print("Number of points in trajectory:",ntestpointsN)
# I just happen t know that there are seven comment lines.
os.makedirs(outdir, exist_ok=True)
Morebins(outdir, nNew,whichFileP,whichFileN,note,ntestpointsN)
input(outdir,nNew,note,nCommentLines)
    
