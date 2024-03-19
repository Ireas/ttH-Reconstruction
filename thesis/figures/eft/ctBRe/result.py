####  Header  #####
####################
import ROOT as root
from ROOT import TFile

def Average(lst):
    return sum(lst)/len(lst)


TFileArray = [TFile("cos_theta_starZ/output.root"), TFile("mz1/output.root"), TFile("mz1_dR/output.root"), TFile("ptz1/output.root"), TFile("ptz1_dR/output.root")]
VariableArray = ["cos_theta_starZ", "mz1", "mz1_dR", "ptz1", "ptz1_dR"]
RegionArray = ["ttZ", "ttW", "Fakes"]

FileOutput = open("Results.txt", "w")


FileOutput.write("======================================================  Ratios  ==================================================\n" )
FileOutput.write("==================================================================================================================\n" )

for i in range( len(TFileArray) ):
    CurrentVariable = VariableArray[i]

    CurrentValueC1 = [0.0, 0.0, 0.0]
    CurrentValueC2 = [0.0, 0.0, 0.0]
    CurrentValueC3 = [0.0, 0.0, 0.0]
    CurrentValueC4 = [0.0, 0.0, 0.0]
    CurrentValueDefault = [0.0, 0.0, 0.0]
    
    CurrentErrorC1 = [0.0, 0.0, 0.0]
    CurrentErrorC2 = [0.0, 0.0, 0.0]
    CurrentErrorC3 = [0.0, 0.0, 0.0]
    CurrentErrorC4 = [0.0, 0.0, 0.0]
    CurrentErrorDefault = [0.0, 0.0, 0.0]
    
    for j in range( len(RegionArray) ):
        CurrentRegionDefault = TFileArray[i].Get(RegionArray[j]+"_default")
        CurrentRegionC1 = TFileArray[i].Get(RegionArray[j]+"_c1")
        CurrentRegionC2 = TFileArray[i].Get(RegionArray[j]+"_c2")
        CurrentRegionC3 = TFileArray[i].Get(RegionArray[j]+"_c3")
        CurrentRegionC4 = TFileArray[i].Get(RegionArray[j]+"_c4")
        
        for k in range(1,CurrentRegionDefault.GetNbinsX()):
            CurrentValueC1[j]+= CurrentRegionC1.GetBinContent(k)
            CurrentValueC2[j]+= CurrentRegionC2.GetBinContent(k)
            CurrentValueC3[j]+= CurrentRegionC3.GetBinContent(k)
            CurrentValueC4[j]+= CurrentRegionC4.GetBinContent(k)
            CurrentValueDefault[j]+= CurrentRegionDefault.GetBinContent(k)

            CurrentErrorC1[j]+= CurrentRegionC1.GetBinError(k)**2
            CurrentErrorC2[j]+= CurrentRegionC2.GetBinError(k)**2
            CurrentErrorC3[j]+= CurrentRegionC3.GetBinError(k)**2
            CurrentErrorC4[j]+= CurrentRegionC4.GetBinError(k)**2
            CurrentErrorDefault[j]+= CurrentRegionDefault.GetBinError(k)**2
        
        CurrentErrorC1[j] = CurrentErrorC1[j]**(1/2)
        CurrentErrorC2[j] = CurrentErrorC2[j]**(1/2)
        CurrentErrorC3[j] = CurrentErrorC3[j]**(1/2)
        CurrentErrorC4[j] = CurrentErrorC4[j]**(1/2)
        CurrentErrorDefault[j] = CurrentErrorDefault[j]**(1/2)


        CurrentValueC1[j] = CurrentValueC1[j]/CurrentValueDefault[j]
        CurrentValueC2[j] = CurrentValueC2[j]/CurrentValueDefault[j]
        CurrentValueC3[j] = CurrentValueC3[j]/CurrentValueDefault[j]
        CurrentValueC4[j] = CurrentValueC4[j]/CurrentValueDefault[j]

        CurrentErrorC1[j] = (CurrentErrorC1[j]/CurrentValueDefault[j])**2 + (CurrentErrorDefault[j]*CurrentValueC1[j]/(CurrentValueDefault[j]**2))**2
        CurrentErrorC2[j] = (CurrentErrorC2[j]/CurrentValueDefault[j])**2 + (CurrentErrorDefault[j]*CurrentValueC2[j]/(CurrentValueDefault[j]**2))**2
        CurrentErrorC3[j] = (CurrentErrorC3[j]/CurrentValueDefault[j])**2 + (CurrentErrorDefault[j]*CurrentValueC3[j]/(CurrentValueDefault[j]**2))**2
        CurrentErrorC4[j] = (CurrentErrorC4[j]/CurrentValueDefault[j])**2 + (CurrentErrorDefault[j]*CurrentValueC4[j]/(CurrentValueDefault[j]**2))**2
        CurrentErrorC1[j] = CurrentErrorC1[j]**(1/2)
        CurrentErrorC2[j] = CurrentErrorC2[j]**(1/2)
        CurrentErrorC3[j] = CurrentErrorC3[j]**(1/2)
        CurrentErrorC4[j] = CurrentErrorC4[j]**(1/2)


    FileOutput.write("======================================================  " +  VariableArray[i] + "  ==================================================\n" )
    FileOutput.write("Coefficient C1: "+  str(CurrentValueC1) + "->\t"+ str(Average(CurrentValueC1)) +"\t+- "+ str(Average(CurrentErrorC1)) +"\n") 
    FileOutput.write("Coefficient C2: "+  str(CurrentValueC2) + "->\t"+ str(Average(CurrentValueC2)) +"\t+- "+ str(Average(CurrentErrorC2)) +"\n") 
    FileOutput.write("Coefficient C3: "+  str(CurrentValueC3) + "->\t"+ str(Average(CurrentValueC3)) +"\t+- "+ str(Average(CurrentErrorC3)) +"\n") 
    FileOutput.write("Coefficient C4: "+  str(CurrentValueC4) + "->\t"+ str(Average(CurrentValueC4)) +"\t+- "+ str(Average(CurrentErrorC4)) +"\n") 
    FileOutput.write("\n\n")




FileOutput.write("\n\n\n\n\n\n\n")
FileOutput.write("======================================================  Sep. Power  ==================================================\n" )
FileOutput.write("======================================================================================================================\n" )

for i in range( len(TFileArray) ):
    CurrentVariable = VariableArray[i]

    CurrentValueC1 = [0.0, 0.0, 0.0]
    CurrentValueC2 = [0.0, 0.0, 0.0]
    CurrentValueC3 = [0.0, 0.0, 0.0]
    CurrentValueC4 = [0.0, 0.0, 0.0]
    
    CurrentErrorC1 = [0.0, 0.0, 0.0]
    CurrentErrorC2 = [0.0, 0.0, 0.0]
    CurrentErrorC3 = [0.0, 0.0, 0.0]
    CurrentErrorC4 = [0.0, 0.0, 0.0]
    
    for j in range( len(RegionArray) ):
        CurrentRegionDefault = TFileArray[i].Get(RegionArray[j]+"_default")
        CurrentRegionC1 = TFileArray[i].Get(RegionArray[j]+"_c1")
        CurrentRegionC2 = TFileArray[i].Get(RegionArray[j]+"_c2")
        CurrentRegionC3 = TFileArray[i].Get(RegionArray[j]+"_c3")
        CurrentRegionC4 = TFileArray[i].Get(RegionArray[j]+"_c4")
        
        for k in range(1,CurrentRegionDefault.GetNbinsX()):
            CurrentValueC1[j]+= ( CurrentRegionC1.GetBinContent(k)-CurrentRegionDefault.GetBinContent(k) )**2/( CurrentRegionC1.GetBinContent(k)+CurrentRegionDefault.GetBinContent(k) )
            CurrentValueC2[j]+= ( CurrentRegionC2.GetBinContent(k)-CurrentRegionDefault.GetBinContent(k) )**2/( CurrentRegionC2.GetBinContent(k)+CurrentRegionDefault.GetBinContent(k) )
            CurrentValueC3[j]+= ( CurrentRegionC3.GetBinContent(k)-CurrentRegionDefault.GetBinContent(k) )**2/( CurrentRegionC3.GetBinContent(k)+CurrentRegionDefault.GetBinContent(k) )
            CurrentValueC4[j]+= ( CurrentRegionC4.GetBinContent(k)-CurrentRegionDefault.GetBinContent(k) )**2/( CurrentRegionC4.GetBinContent(k)+CurrentRegionDefault.GetBinContent(k) )
           
            CurrentErrorC1[j]+= CurrentRegionC1.GetBinError(k)**2 * ( (2)*(CurrentRegionC1.GetBinContent(k)-CurrentRegionDefault.GetBinContent(k))/(CurrentRegionC1.GetBinContent(k)+CurrentRegionDefault.GetBinContent(k)) + (-1)*(CurrentRegionC1.GetBinContent(k)-CurrentRegionDefault.GetBinContent(k))**2/((CurrentRegionC1.GetBinContent(k)+CurrentRegionDefault.GetBinContent(k))**2) )**2
            CurrentErrorC1[j]+= CurrentRegionDefault.GetBinError(k)**2 * ( (-2)*(CurrentRegionC1.GetBinContent(k)-CurrentRegionDefault.GetBinContent(k))/(CurrentRegionC1.GetBinContent(k)+CurrentRegionDefault.GetBinContent(k)) + (-1)*(CurrentRegionC1.GetBinContent(k)-CurrentRegionDefault.GetBinContent(k))**2/((CurrentRegionC1.GetBinContent(k)+CurrentRegionDefault.GetBinContent(k))**2) )**2
            CurrentErrorC1[j] = CurrentErrorC1[j]**(1/2)

            CurrentErrorC2[j]+= CurrentRegionC2.GetBinError(k)**2 * ( (2)*(CurrentRegionC2.GetBinContent(k)-CurrentRegionDefault.GetBinContent(k))/(CurrentRegionC2.GetBinContent(k)+CurrentRegionDefault.GetBinContent(k)) + (-1)*(CurrentRegionC2.GetBinContent(k)-CurrentRegionDefault.GetBinContent(k))**2/((CurrentRegionC2.GetBinContent(k)+CurrentRegionDefault.GetBinContent(k))**2) )**2
            CurrentErrorC2[j]+= CurrentRegionDefault.GetBinError(k)**2 * ( (-2)*(CurrentRegionC2.GetBinContent(k)-CurrentRegionDefault.GetBinContent(k))/(CurrentRegionC2.GetBinContent(k)+CurrentRegionDefault.GetBinContent(k)) + (-1)*(CurrentRegionC2.GetBinContent(k)-CurrentRegionDefault.GetBinContent(k))**2/((CurrentRegionC2.GetBinContent(k)+CurrentRegionDefault.GetBinContent(k))**2) )**2
            CurrentErrorC2[j] = CurrentErrorC2[j]**(1/2)
            
            CurrentErrorC3[j]+= CurrentRegionC3.GetBinError(k)**2 * ( (2)*(CurrentRegionC3.GetBinContent(k)-CurrentRegionDefault.GetBinContent(k))/(CurrentRegionC3.GetBinContent(k)+CurrentRegionDefault.GetBinContent(k)) + (-1)*(CurrentRegionC3.GetBinContent(k)-CurrentRegionDefault.GetBinContent(k))**2/((CurrentRegionC3.GetBinContent(k)+CurrentRegionDefault.GetBinContent(k))**2) )**2
            CurrentErrorC3[j]+= CurrentRegionDefault.GetBinError(k)**2 * ( (-2)*(CurrentRegionC3.GetBinContent(k)-CurrentRegionDefault.GetBinContent(k))/(CurrentRegionC3.GetBinContent(k)+CurrentRegionDefault.GetBinContent(k)) + (-1)*(CurrentRegionC3.GetBinContent(k)-CurrentRegionDefault.GetBinContent(k))**2/((CurrentRegionC3.GetBinContent(k)+CurrentRegionDefault.GetBinContent(k))**2) )**2
            CurrentErrorC3[j] = CurrentErrorC3[j]**(1/2)
            
            CurrentErrorC4[j]+= CurrentRegionC4.GetBinError(k)**2 * ( (2)*(CurrentRegionC4.GetBinContent(k)-CurrentRegionDefault.GetBinContent(k))/(CurrentRegionC4.GetBinContent(k)+CurrentRegionDefault.GetBinContent(k)) + (-1)*(CurrentRegionC4.GetBinContent(k)-CurrentRegionDefault.GetBinContent(k))**2/((CurrentRegionC4.GetBinContent(k)+CurrentRegionDefault.GetBinContent(k))**2) )**2
            CurrentErrorC4[j]+= CurrentRegionDefault.GetBinError(k)**2 * ( (-2)*(CurrentRegionC4.GetBinContent(k)-CurrentRegionDefault.GetBinContent(k))/(CurrentRegionC4.GetBinContent(k)+CurrentRegionDefault.GetBinContent(k)) + (-1)*(CurrentRegionC4.GetBinContent(k)-CurrentRegionDefault.GetBinContent(k))**2/((CurrentRegionC4.GetBinContent(k)+CurrentRegionDefault.GetBinContent(k))**2) )**2
            CurrentErrorC4[j] = CurrentErrorC4[j]**(1/2)


        CurrentValueC1[j]/= 2
        CurrentValueC2[j]/= 2
        CurrentValueC3[j]/= 2
        CurrentValueC4[j]/= 2
        
        CurrentErrorC1[j]/= 2
        CurrentErrorC2[j]/= 2
        CurrentErrorC3[j]/= 2
        CurrentErrorC4[j]/= 2

    FileOutput.write("======================================================  " +  VariableArray[i] + "  ==================================================\n" )
    FileOutput.write("Separation Power C1: "+  str(CurrentValueC1) + "->\t"+ str(Average(CurrentValueC1)) +"\t+- " + str(Average(CurrentErrorC1)) +"\n") 
    FileOutput.write("Separation Power C2: "+  str(CurrentValueC2) + "->\t"+ str(Average(CurrentValueC2)) +"\t+- " + str(Average(CurrentErrorC2)) +"\n") 
    FileOutput.write("Separation Power C3: "+  str(CurrentValueC3) + "->\t"+ str(Average(CurrentValueC3)) +"\t+- " + str(Average(CurrentErrorC3)) +"\n") 
    FileOutput.write("Separation Power C4: "+  str(CurrentValueC4) + "->\t"+ str(Average(CurrentValueC4)) +"\t+- " + str(Average(CurrentErrorC4)) +"\n") 
    FileOutput.write("\n\n")






FileOutput.close()
