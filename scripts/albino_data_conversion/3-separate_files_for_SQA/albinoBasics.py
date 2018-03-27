__author__ = 'yzhu14'
import numpy as np
from unit_const import *
import os
import sys
from scipy.interpolate import interp1d
#Input from Albino Perego's nu_potential_V2, v_density/v_potential
#out put neutrino self interaciton potentials of 8/n bins for neutrino ocsillation calculation
#use tracers from Dirk's paper, get neutrino density files from Albino Perego's nu_potential_V2
#and use this script prepare input for neutrino oscillation code
# energy Bin interpolation from 8 bins from Albino's simulation(nu_potential_v2.0)
# For neutrino capture rate calculation

#######################################################################
######################################### getnewbins #################################################################
#   Yonglin Zhu 2017
############################################################################
def getNewE(nNew):
    if nNew ==1:
        Enew =np.array([np.exp(np.log(2 * M) + (np.log(37.48 * M) - np.log(2 * M))/2)])
        dEnew = np.array([0.])
    else:
        Enew = np.array([np.exp(np.log(2 * M) + i * ((np.log(37.48 * M) - np.log(2 * M)) / (nNew - 1))) for i in range(nNew)])
        dEnew = np.array([np.exp(np.log(2*M)+(i+0.5)*((np.log(37.48*M)-np.log(2*M))/(nNew-1)))-np.exp(np.log(2*M)+(i-0.5)*((np.log(37.48*M)-np.log(2*M))/(nNew-1))) for i in range(nNew)])
    dEnew/=M
    Enew /=M
    EnewCgs=Enew*MeV
    ECgs =E*MeV #erg
    dECgs = dE*MeV #erg
    dEnewCgs = dEnew*MeV #erg
    return Enew,dEnew,EnewCgs,ECgs,dECgs,dEnewCgs
############################################################################
#   Yonglin Zhu 2017
############################################################################
def test(nf,nNew,whichFileP,whichFileN,note):
    #test normalization: neutrino capture files are good for different energy bins.
    #open density file, different nE files, get the minimum line number, generate random number of line for checking
    #if raduis is the same, check the total neutrino density normalization
    #if some thing is wrong, double the test case number, if still error, print number of nE, line, and the not normalized values.
    infileNC=open(dir + '/' + str(nNew) + 'NE/'  + str(nf) + note + 'nucap.txt','r')
    nlineNC=sum(1 for line in open(dir + '/' + str(nNew) + 'NE/'  + str(nf) + note + 'nucap.txt'))
    linenc = infileNC.readlines()
    inFileDensity = open(dir +'/' + str(nNew) + 'NE/' + whichFileN+ str(nf)  +  note + '.txt','r')
    nlineD=sum(1 for line in open(dir + '/' + str(nNew) + 'NE/'  +whichFileP+ str(nf) + note + '.txt'))
    linedensity = inFileDensity.readlines()  # yzhu: START FROM 0
    inFileDensity.close()
    nline=min(nlineD,nlineNC)
    ncheck=10

######################################## ve ###########################################################################################
#   Yonglin Zhu 2017
#   input: nf - file #, Ye/rho files, nNew - new # of energy bins,
#   output: ve
############################################################################################################################################
def ve(nf,nNew,whichFileP,note):
    inFileYe = open(dir +'/' + str(nNew) + 'NE/'  + 'Ye_' + whichFileP + str(nf) + note + '.txt', 'r')
    inFilerho = open(dir + '/' + str(nNew) + 'NE/'  +'rho_' + whichFileP + str(nf) + note + '.txt', 'r')
    outFileve=open(dir + '/' + str(nNew) + 'NE/'  + 'Ve_' + whichFileP + str(nf) + note + '.txt', 'w')
    lineYe=inFileYe.readlines()
    linerho=inFilerho.readlines()
    ntestpointsrho = sum(1 for line in open(dir +'/' + str(nNew) + 'NE/'  + 'rho_' + whichFileP + str(nf) + note + '.txt'))
    ntestpointsYe = sum(1 for line in open(dir +'/' + str(nNew) + 'NE/'  + 'Ye_' + whichFileP + str(nf) + note + '.txt'))
    npoints =min(ntestpointsrho,ntestpointsYe)
    for k in range(0,npoints):
        yepoints=lineYe[k]
        yeline=yepoints.split()
        rhopoints=linerho[k]
        rholine=rhopoints.split()
        outFileve.write(str(rholine[0])+"\t"+str(float(rholine[1])*float(yeline[1])*makeAPotential*NA)+"\n")
    inFilerho.close()
    inFileYe.close()
    outFileve.close()

def onebin(dir,nf,nNew,whichFileP,whichFileN,note,ntestpoints):
    nucap=0
    Enew,dEnew,EnewCgs,ECgs,dECgs,dEnewCgs=getNewE(nNew)
    inFileP = open(str(dir) + '/' + str(whichFileP) + str(nf) + '.txt', 'r')
    inFileN = open(str(dir) + '/' + str(whichFileN) + str(nf) + '.txt', 'r')
    linesP = inFileP.readlines()  # yzhu: START FROM 0
    linesN = inFileN.readlines()  # yzhu: START FROM 0
    inFileN.close()
    inFileP.close()
    outFileYe = open(dir + '/' + str(nNew) + 'NE/'  +'Ye_' + whichFileP + str(nf) + note + '.txt', 'w')
    outFilerho = open(dir +  '/' + str(nNew) + 'NE/' + 'rho_' + whichFileP + str(nf) + note + '.txt', 'w')
    outFilePotentiale = open(dir + '/' + str(nNew) + 'NE/1'  +whichFileP+ str(nf) + note + '.txt',
                            'w')  # plot linear and cubic spline
    outFileDensitye = open(dir +'/' + str(nNew) + 'NE/1' + whichFileN+ str(nf)  +  note + '.txt',
                          'w')  # plot linear and cubic spline
    outFilePotentiala = open(dir + '/' + str(nNew) + 'NE/2'  +whichFileP+ str(nf) + note + '.txt',
                         'w')  # plot linear and cubic spline
    outFileDensitya = open(dir +'/' + str(nNew) + 'NE/2' + whichFileN+ str(nf)  +  note + '.txt',
                      'w')  # plot linear and cubic spline
    outFilePotentialx = open(dir + '/' + str(nNew) + 'NE/3'  +whichFileP+ str(nf) + note + '.txt',
                        'w')  # plot linear and cubic spline
    outFileDensityx = open(dir +'/' + str(nNew) + 'NE/3' + whichFileN+ str(nf)  +  note + '.txt',
                      'w')  # plot linear and cubic spline
    outFileErr = open(dir + '/' + str(nNew) + 'NE/' +  'dir' + whichFileP + str(nf) + note + 'err.txt', 'w')
    outFileIntpCheck = open(dir  + '/' + str(nNew) + 'NE/'  + str(nf) + 'intpCheckffcub' + whichFileP+note + '.txt',
                            'w')  # plot linear and cubic spline
    outfileNC=open(dir + '/' + str(nNew) + 'NE/'  + str(nf) + note + 'nucap.txt','w')
    p = [0] * nEnergyBins
    pnew = [0]*numberofneutrinotype*nNew
    d = [0] * nEnergyBins
    dnew = [0]*numberofneutrinotype*nNew
    for n in range(0,ntestpoints):  # all points in the test trajectory
        Rcore=0. #from the core, use this for neutrino caputure rate extension
        Rs=0. #from the stating point, usually used for neutrino oscillation
        totalPotential=[0.]*numberofneutrinotype
        totalDensity=[0.]*numberofneutrinotype
        newtotalPotential=[0.]*numberofneutrinotype
        newtotalDensity=[0.]*numberofneutrinotype
        coordinatesLinesP = linesP[n]
        coordinatesLinesN = linesN[n]
        pointsCoordinatesP = coordinatesLinesP.split()
        pointsCoordinatesN = coordinatesLinesN.split()
        ################# neutrino caputure rate ##############
        if nucap==1:
            noOscillationRateMatter     = 0.
            currentRateMatter     = 0.
            noOscillationRateAntiMatter = 0.
            currentRateAntiMatter = 0.
            temp = 0.
            eDensity = np.array([0.0 for i in range(nNew)])
            eBarDensity = np.array([0.0 for i in range(nNew)])
            heavyDensity = np.array([0.0 for i in range(nNew)])
        ################# neutrino caputure rate ##############
        #get initial points of trajectory
        if n == 0:
            x0 = float(pointsCoordinatesP[0])
            y0 = float(pointsCoordinatesP[1])
            z0 = float(pointsCoordinatesP[2])
            #use R from the core
            #print(x0,y0,z0)
        Rcore = km * math.sqrt((float(pointsCoordinatesP[0])) ** 2 +
                               (float(pointsCoordinatesP[1])) ** 2 +
                               (float(pointsCoordinatesP[2])) ** 2)
        Rs = km * math.sqrt((float(pointsCoordinatesP[0]) - x0) ** 2 +
                            (float(pointsCoordinatesP[1])-y0) ** 2 +
                            (float(pointsCoordinatesP[2])-z0) ** 2)
        if n==0:
            print("Trajectory",nf, "Starting at(km):",float(pointsCoordinatesP[0]),float(pointsCoordinatesP[1]),float(pointsCoordinatesP[2]))
        outFileYe.write(str(Rs) + "\t" + str(pointsCoordinatesP[5]) + "\n")
        outFilerho.write(str(Rs) + "\t" + str(pointsCoordinatesP[3]) + "\n")
        checkP=[0.]*numberofneutrinotype
        checkD=[0.]*numberofneutrinotype
        for it in range(0, numberofneutrinotype):
            for j in range(0, nEnergyBins):
                p[j] = makeAPotential * float(pointsCoordinatesP[6 + j + it * nEnergyBins])/dECgs[j]
                d[j] = float(pointsCoordinatesN[6 + j + it * nEnergyBins])/dECgs[j]
                totalPotential[it]=totalPotential[it]+p[j]*dECgs[j] #/dE
                totalDensity[it]=totalDensity[it]+d[j]*dECgs[j]
        outFileDensitye.write(str(Rs) + "\t" +str(float(totalDensity[0]))+"\n")
        outFilePotentiale.write(str(Rs) + "\t" +str(float(totalPotential[0]))+"\n")
        outFileDensitya.write(str(Rs) + "\t" +str(float(totalDensity[1]))+"\n")
        outFilePotentiala.write(str(Rs) + "\t" +str(float(totalPotential[1]))+"\n")
        outFileDensityx.write(str(Rs) + "\t" +str(float(totalDensity[2]))+"\n")
        outFilePotentialx.write(str(Rs) + "\t" +str(float(totalPotential[2]))+"\n")
    outFileIntpCheck.close()
    outFilePotentiale.close()
    outFileDensitye.close()
    outFilePotentiala.close()
    outFileDensitya.close()
    outFilePotentialx.close()
    outFileDensityx.close()

######################################### Morebins #################################################################
#   Yonglin Zhu 2017
#   input: nf - file#, NE - # of energy bins, ntestpoints - # of line from input files, nNew - new # of energy bins
#   output: potential/density files with NE energy bins, Ye/rho files
#   check: R_original(core or not), potential/density units, potential/density normalization, interplote algorithm
#   to do: checkmorebins()
####################################################################################################################
def Morebins(outdir,nNew,whichFileP,whichFileN,note,ntestpoints):
    nucap=1
    Enew,dEnew,EnewCgs,ECgs,dECgs,dEnewCgs=getNewE(nNew)
    inFileP = open(str(whichFileP), 'r')
    inFileN = open(str(whichFileN), 'r')
    linesP = inFileP.readlines()  # yzhu: START FROM 0
    linesN = inFileN.readlines()  # yzhu: START FROM 0
    inFileN.close()
    inFileP.close()
    outFileTemp = open(outdir + '/temp.txt', 'w')
    outFileYe = open(outdir + '/Ye.txt', 'w')
    outFilerho = open(outdir +  '/rho.txt', 'w')
    outFilePotential = open(outdir + '/potential.txt', 'w')  # plot linear and cubic spline
    outFileDensity = open(outdir +'/density.txt', 'w')  # plot linear and cubic spline
    outFileErr = open(outdir + '/potential_err.txt', 'w')
    outFileIntpCheck = open(outdir  + '/intpCheckffcub_potential.txt','w')  # plot linear and cubic spline
    outfileNC=open(outdir + '/nucap.txt','w')
    p = [0] * nEnergyBins
    pnew = [0]*numberofneutrinotype*nNew
    d = [0] * nEnergyBins
    dnew = [0]*numberofneutrinotype*nNew
    for n in range(0,ntestpoints):  # all points in the test trajectory
        Rcore=0. #from the core, use this for neutrino caputure rate extension
        Rs=0. #from the stating point, usually used for neutrino oscillation
        totalPotential=[0.]*numberofneutrinotype
        totalDensity=[0.]*numberofneutrinotype
        newtotalPotential=[0.]*numberofneutrinotype
        newtotalDensity=[0.]*numberofneutrinotype
        coordinatesLinesP = linesP[n]
        coordinatesLinesN = linesN[n]
        pointsCoordinatesP = coordinatesLinesP.split()
        pointsCoordinatesN = coordinatesLinesN.split()
        ################# neutrino caputure rate ##############
        if nucap==1:
            noOscillationRateMatter     = 0.
            currentRateMatter     = 0.
            noOscillationRateAntiMatter = 0.
            currentRateAntiMatter = 0.
            temp = 0.
            eDensity = np.array([0.0 for i in range(nNew)])
            eBarDensity = np.array([0.0 for i in range(nNew)])
            heavyDensity = np.array([0.0 for i in range(nNew)])
        ################# neutrino caputure rate ##############
        #get initial points of trajectory
        if n == 0:
            x0 = float(pointsCoordinatesP[0])
            y0 = float(pointsCoordinatesP[1])
            z0 = float(pointsCoordinatesP[2])
            #use R from the core
            #print(x0,y0,z0)
        Rcore = km * math.sqrt((float(pointsCoordinatesP[0])) ** 2 +
                               (float(pointsCoordinatesP[1])) ** 2 +
                               (float(pointsCoordinatesP[2])) ** 2)
        Rs = km * math.sqrt((float(pointsCoordinatesP[0]) - x0) ** 2 +
                               (float(pointsCoordinatesP[1])-y0) ** 2 +
                               (float(pointsCoordinatesP[2])-z0) ** 2)
        if n==0:
            print("Trajectory Starting at(km):",float(pointsCoordinatesP[0]),float(pointsCoordinatesP[1]),float(pointsCoordinatesP[2]))
        outFileYe.write(str(Rs) + "\t" + str(pointsCoordinatesP[5]) + "\n")
        outFilerho.write(str(Rs) + "\t" + str(pointsCoordinatesP[3]) + "\n")
        outFileTemp.write(str(Rs) + "\t" + str(pointsCoordinatesP[4]) + "\n")
        outFilePotential.write(str(Rs) + "\t")
        outFileDensity.write(str(Rs)+"\t")
        checkP=[0.]*numberofneutrinotype
        checkD=[0.]*numberofneutrinotype
        for it in range(0, numberofneutrinotype):
            for j in range(0, nEnergyBins):
                p[j] = makeAPotential * float(pointsCoordinatesP[6 + j + it * nEnergyBins])/dECgs[j]
                d[j] = float(pointsCoordinatesN[6 + j + it * nEnergyBins])/dECgs[j]
                totalPotential[it]=totalPotential[it]+p[j]*dECgs[j] #/dE
                totalDensity[it]=totalDensity[it]+d[j]*dECgs[j]
                #print(it,j,d[j],totalDensity[it])
            fd=interp1d(E,d)
            fp=interp1d(E,p)
            for k in range(0, nNew):
                if k == 0:
                    pnew[k]=p[k]*dEnewCgs[k] #outFilePotential.write(str(p[k]) + "\t")  # first use the uninterpolated value
                    dnew[k]=d[k]*dEnewCgs[k]
                elif k == nNew - 1:
                    pnew[k]=p[nEnergyBins-1]*dEnewCgs[k]#outFilePotential.write(str(p[nEnergyBins - 1]) + "\t")
                    dnew[k]=d[nEnergyBins-1]*dEnewCgs[k]
                else:
                    pnew[k]=fp(Enew[k])*dEnewCgs[k] #outFilePotential.write(str(fp(Enew[k])) + "\t")
                    dnew[k]=fd(Enew[k])*dEnewCgs[k]
                newtotalPotential[it]=newtotalPotential[it]+ pnew[k]
                newtotalDensity[it]=newtotalDensity[it]+dnew[k]
                #print("new",k,dnew[k])
                #####################
            if newtotalDensity[it]==0 or newtotalPotential[it]==0:
                outFileErr.write("0 newtotal\t"+str(Rs)+"\t"+str(it)+"\n")
                non0=0 #record the 0 potential and out put after that
                for k in range(0,nNew):
                    outFilePotential.write(str(pnew[k]) + "\t")
                    checkP[it]+=float(pnew[k])
                    outFileDensity.write(str(dnew[k])+ "\t")
                    checkD[it]+=float(dnew[k])
            else:
                for k in range(0,nNew):
                    pnew[k]*=totalPotential[it]/newtotalPotential[it]
                    dnew[k]*=totalDensity[it]/newtotalDensity[it]
                    outFilePotential.write(str(pnew[k]) + "\t")
                    checkP[it]+=float(pnew[k])
                    outFileDensity.write(str(dnew[k])+ "\t")
                    checkD[it]+=float(dnew[k])
        ################# neutrino caputure rate ##############
                    if it ==0:
                        eDensity[k]= dnew[k]
                    elif it ==1:
                        eBarDensity[k]=dnew[k]
                    else:
                        heavyDensity[k]=dnew[k]
        #print(str(float(totalDensity[0]))+"\t"+str(float(newtotalDensity[0]))+"\t"+str(checkD[0])+"\t")
        if nucap==1:
          #  eDensity[k]=float(pnew[1 + k + 0 * nNew])
          #  eBarDensity[k]=float(pointsDensity[1 + k + 1 * nNew])
          #  heavyDensity[k]=float(pointsDensity[1 + k + 2 * nNew])
            for k in range(0,nNew):
                temp = eDensity[k]*(EnewCgs[k] + nucleonMassDifference)*(EnewCgs[k] + nucleonMassDifference)
                nm   = 1.0 #norm(Sf[e][e])
                nmSfemu=0.0
                noOscillationRateMatter = temp +noOscillationRateMatter
                currentRateMatter       = nm*temp+currentRateMatter+ (nmSfemu)*(heavyDensity[k])*(EnewCgs[k] + nucleonMassDifference)*(EnewCgs[k] + nucleonMassDifference)
                if(EnewCgs[k] > thresholdEnergy):
                    temp = eBarDensity[k]*(EnewCgs[k] - nucleonMassDifference)*(EnewCgs[k] - nucleonMassDifference)
                    nm   = 1.0 #norm(Sfbar[e][e])
                    nmSfbaremu = 0.0 #norm(Sfbar[e][mu])
                    noOscillationRateAntiMatter += temp
                    currentRateAntiMatter       += nm*temp + (nmSfbaremu)*(heavyDensity[k])*(EnewCgs[k] - nucleonMassDifference)*(EnewCgs[k] - nucleonMassDifference)
        ################# neutrino caputure rate ##############
        #print(str(float(totalDensity[2]))+"\t"+str(float(newtotalDensity[2]))+"\t"+str(checkD[2])+"\t")
        if nucap==1:
            #print(str(Rcore),str(pointsCoordinatesN[2]))
            outfileNC.write(str(Rcore)+"\t"+str(noOscillationRateMatter*Cnunu)+'\t'+str(currentRateMatter*Cnunu)+"\t"+str(noOscillationRateAntiMatter*Cnunu)
                            +"\t"+str(currentRateAntiMatter*Cnunu)+"\n")
        ################# neutrino caputure rate ##############

        outFileDensity.write(str(float(totalDensity[0]))+"\t"+str(float(newtotalDensity[0]))+"\t"+str(checkD[0])+"\t"+
                             str(float(totalDensity[1]))+"\t"+str(float(newtotalDensity[1]))+"\t"+str(checkD[1])+"\t"+
                             str(float(totalDensity[2]))+"\t"+str(float(newtotalDensity[2]))+"\t"+str(checkD[2])+"\t"+"\n")
        outFilePotential.write(str(float(totalPotential[0]))+"\t"+str(float(newtotalPotential[0]))+"\t"+str(checkP[0])+"\t"+
                               str(float(totalPotential[1]))+"\t"+str(float(newtotalPotential[1]))+"\t"+str(checkP[1])+"\t"+
                               str(float(totalPotential[2]))+"\t"+str(float(newtotalPotential[2]))+"\t"+str(checkP[2])+"\t"+"\n")
    outFileIntpCheck.close()
    outFilePotential.close()
    outFileDensity.close()
    outFileTemp.close()
######################################## neucap ###########################################################################################
#   Yonglin Zhu 2017
#   input: nf - file #, Density files, nNew - new # of energy bins, NE - # of energy bins,
#   output: neutrino caputure rate files,
############################################################################################################################################
def neucap(dir,nf,nNew,whichFileP,whichFileN,note):
    Enew,dEnew,EnewCgs,ECgs,dECgs,dEnewCgs=getNewE(nNew)
    outfileNC=open(dir + '/' + str(nNew) + 'NE/'  + str(nf) + note + 'nucap.txt','w')
    inFileDensity = open(dir +'/' + str(nNew) + 'NE/' + whichFileN+ str(nf)  +  note + '.txt','r')
    nline=sum(1 for line in open(dir + '/' + str(nNew) + 'NE/'  +whichFileP+ str(nf) + note + '.txt'))
    linedensity = inFileDensity.readlines()  # yzhu: START FROM 0
    inFileDensity.close()
    inFileN = open(dir + '/' + whichFileN + str(nf) + '.txt', 'r')
    linesN = inFileN.readlines()  # yzhu: START FROM 0
    inFileN.close()
    eDensity = np.array([0.0 for i in range(nNew)])
    eBarDensity = np.array([0.0 for i in range(nNew)])
    heavyDensity = np.array([0.0 for i in range(nNew)])
    for n in range(0,nline):
        coordinatesLinesN = linesN[n]
        pointsCoordinatesN = coordinatesLinesN.split()
############################
        noOscillationRateMatter     = 0.
        currentRateMatter     = 0.
        noOscillationRateAntiMatter = 0.
        currentRateAntiMatter = 0.
        temp = 0.
        coordinatesdensity = linedensity[n]
        pointsDensity = coordinatesdensity.split()
        for i in range(0, nNew):
            ################check the order of nue and nuebar
            eDensity[i]=float(pointsDensity[1 + i + 0 * nNew])
            eBarDensity[i]=float(pointsDensity[1 + i + 1 * nNew])
            heavyDensity[i]=float(pointsDensity[1 + i + 2 * nNew])
            #print(eDensity)
            temp = eDensity[i]*(EnewCgs[i] + nucleonMassDifference)*(EnewCgs[i] + nucleonMassDifference)
            nm   = 1.0 #norm(Sf[e][e])
            nmSfemu=0.0
            noOscillationRateMatter = temp +noOscillationRateMatter
            currentRateMatter       = nm*temp+currentRateMatter+ (nmSfemu)*(heavyDensity[i])*(EnewCgs[i] + nucleonMassDifference)*(EnewCgs[i] + nucleonMassDifference)
            if(EnewCgs[i] > thresholdEnergy):
                temp = eBarDensity[i]*(EnewCgs[i] - nucleonMassDifference)*(EnewCgs[i] - nucleonMassDifference)
                nm   = 1.0 #norm(Sfbar[e][e])
                nmSfbaremu = 0.0 #norm(Sfbar[e][mu])
                noOscillationRateAntiMatter += temp
                currentRateAntiMatter       += nm*temp + (nmSfbaremu)*(heavyDensity[i])*(EnewCgs[i] - nucleonMassDifference)*(EnewCgs[i] - nucleonMassDifference)
        print(str(pointsDensity[0]),str(pointsCoordinatesN[2]))
        outfileNC.write(str(pointsCoordinatesN[2])+"\t"+str(noOscillationRateMatter*Cnunu)+'\t'+str(currentRateMatter*Cnunu)+"\t"+str(noOscillationRateAntiMatter*Cnunu)
                        +"\t"+str(currentRateAntiMatter*Cnunu)+"\n")

######################################## input ###########################################################################################
#   Yonglin Zhu 2017
#   input: nf - file #, new Density files, nNew - new # of energy bins,
#   output: outFilePotentialBin, outFileDensityBin
############################################################################################################################################
def input(outdir,nNew,whichFileP,whichFileN,note,nCommentLines):
    #outFile = open(outdir + "/"  +whichFileP + str(nf) + note + '.txt', 'w')
    inFileNewP = open(outdir+'/potential.txt','r')  # plot linear and cubic spline
    linesNewP = inFileNewP.readlines()  # yzhu: START FROM 0
    inFileNewP.close()
    inFileNewN = open(outdir+'/density.txt', 'r')  # plot linear and cubic spline
    linesNewN = inFileNewN.readlines()  # yzhu: START FROM 0
    inFileNewN.close()
    ntestpoints = sum(1 for line in open(whichFileP))
    for it in range(0, numberofneutrinotype):
        print("it=",it)
        for k in range(0, nNew):
            print("k=",k)
            outFilePotentialBin = open(outdir + '/potential_s' + str(it + 1) + '_g' + str(
                k + 1) + '_' + note + '.txt', 'w')
            outFileDensityBin = open(outdir + '/density_s' + str(it + 1) + '_g' + str(
                k + 1) + '_' + note + '.txt', 'w')
            for n in range(nCommentLines, ntestpoints):  # all points in the test trajectory
                coordinatesLinesNewP = linesNewP[n]
                pointsCoordinatesNewP = coordinatesLinesNewP.split()
                outFilePotentialBin.write(
                    str(pointsCoordinatesNewP[0]) + "\t"
                    + str(float(pointsCoordinatesNewP[1 + k + it * nNew]))+ "\n")
                coordinatesLinesNewN = linesNewN[n]
                pointsCoordinatesNewN = coordinatesLinesNewN.split()
                outFileDensityBin.write(
                    str(pointsCoordinatesNewN[0]) + "\t"
                    + str(float(pointsCoordinatesNewN[1 + k + it * nNew]))+ "\n")
    outFilePotentialBin.close()
    #outFile.close()
############################################################
###############
