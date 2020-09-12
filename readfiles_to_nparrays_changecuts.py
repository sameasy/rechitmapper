#   python3 readfiles_to_nparrays.py -I file.txt -T photon
#   python3 readfiles_to_nparrays.py -I file.txt -T beamhalo
import matplotlib
matplotlib.use('Agg')
import uproot
import numpy as np
import pandas as pd
import root_numpy as rp
import matplotlib.pyplot as plt
import ROOT
from tqdm import tqdm
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-I", "--inputfilelist",help="input data list",default="file.txt")
parser.add_argument("-T", "--stype",help="photon / beamhalo",default="photon")


flist = parser.parse_args().inputfilelist
sampletype = parser.parse_args().stype

savename = sampletype


def make_nparray(nev):
    histo_2d = ROOT.TH2F("histo_2d", "X vs Y", 11, -5, 6, 11, -5, 6)
    imax = np.argmax(dfs['eerhE'][nev])
#    print ("max energy",dfs['eerhE'][nev][imax])
    maxval = dfs['eerhE'][nev][imax]
#    print ("max index=",imax," ieta:",dfs['eerhiEta'][nev][imax]," iphi:",dfs['eerhiPhi'][nev][imax]) 
    ietamax = int(dfs['eerhiEta'][nev][imax])
    iphimax = int(dfs['eerhiPhi'][nev][imax])
    zsideindex = (dfs['eerhEta'][nev][imax])

    for i in range(len(dfs['eerhE'][nev])):
        if (np.sign(zsideindex) != np.sign(dfs['eerhEta'][nev][i])): 
#            print ("no sign match")
            continue
#for checks print e,eta,phi for selected event        
#        if (nev==1):
#            print ("energy :",dfs['eerhE'][nev][i],"ieta :",dfs['eerhiEta'][nev][i],"iphi :",dfs['eerhiPhi'][nev][i])
#        print ("sign match")
        binh = histo_2d.GetBin(histo_2d.GetXaxis().FindBin(int(dfs['eerhiEta'][nev][i]-ietamax)),histo_2d.GetYaxis().FindBin(int(dfs['eerhiPhi'][nev][i]-iphimax)))
        histo_2d.SetBinContent(binh,dfs['eerhE'][nev][i]/maxval) #if normalized by max
#        histo_2d.SetBinContent(binh,dfs['eerhE'][nev][i])
    return rp.hist2array(histo_2d, include_overflow=False, copy=True, return_edges=False)

def save_arr(listarr,filename):
    np_arr = np.array(listarr)
    print ("shape pf np arrays:",np_arr.shape)
    np.save(filename, np_arr)        


def plot_images(listarr,filename,nimages):
    array = np.array(listarr[:nimages])
    fig = plt.figure()
    for i in range(nimages):
        subfigr = fig.add_subplot()
        subfigr.imshow(array[i])
        subfigr.colorbar()
#        plt.show()
    fig.savefig(filename)


testfile = open(flist, 'r')
fname_list = testfile.readlines()
testfile.close()
print ("finished reading file list")
print (fname_list)

print ("starting processing arrays")
listarr = []


# In[2]:


#fName='root://cmseos.fnal.gov//store/user/sghosh/FROMSHILPIDI/anTGCtree_data_photon.root'
#fName='root://cmseos.fnal.gov//store/user/sghosh/FROMSHILPIDI/anTGCtree_MC_beamhalo.root'
#tree = uproot.open(fName,xrootdsource=dict(chunkbytes=1024**3, limitbytes=1024**3))['ggNtuplizer/EventTree']


for nfile in tqdm(fname_list):
    print ("processing file :",str(nfile))
    fName=nfile.rstrip('\n')
    tree = uproot.open(fName,xrootdsource=dict(chunkbytes=1024**3, limitbytes=1024**3))['ggNtuplizer/EventTree']
    branches = ['eerhE','eerhiEta','eerhiPhi','eerhEta','pfMET','phoE','phoEta','phoPhi','phoHoverE','phoPFChWorstIso','HLTPho']
    dictn = tree.arrays(branches=branches)
    dfs = pd.DataFrame.from_dict(dictn)
    dfs.columns=branches
    phoindex = 0 #0 for the highest enegetic photon
    dfs['HLTdecision']=dfs['HLTPho'].apply(lambda x : (x>>9)&1)
    dfs['phoE']=dfs['phoE'].apply(pd.Series)[phoindex]
    dfs['phoEta']=dfs['phoEta'].apply(pd.Series)[phoindex]
    dfs['phoPhi']=dfs['phoPhi'].apply(pd.Series)[phoindex]
    dfs['phoHoverE']=dfs['phoHoverE'].apply(pd.Series)[phoindex]
    dfs['phoPFChWorstIso']=dfs['phoPFChWorstIso'].apply(pd.Series)[phoindex]

#    dfs['eerhEmax']=dfs['eerhE'].apply(lambda x : np.max(x))
#    photoncut = (dfs['pfMET'] > 80.0)&(dfs['phoHoverE'] < 0.0260)&(dfs['phoPFChWorstIso'] < 1.146)&(abs(dfs['phoEta']) > 1.566)&(((abs(dfs['phoPhi'])<2.9)&(abs(dfs['phoPhi'])>3.0))|((abs(dfs['phoPhi'])<0)&(abs(dfs['phoPhi'])>0.2)))
    photoncut = (dfs['pfMET'] > 80.0)&(dfs['phoHoverE'] < 0.0260)&(dfs['phoPFChWorstIso'] < 1.146)&(abs(dfs['phoEta']) > 1.566)&(((abs(dfs['phoPhi'])<2.9))&((abs(dfs['phoPhi'])>0.2)))
    halocut = (dfs['pfMET'] > 150.0)&(dfs['HLTdecision'] ==1)&(abs(dfs['phoEta']) > 1.566)&(((abs(dfs['phoPhi'])>2.9)&(abs(dfs['phoPhi'])<3.0))|((abs(dfs['phoPhi'])>0)&(abs(dfs['phoPhi'])<0.2))) 
    #photoncut = (dfs['pfMET'] > 150.0)&(dfs['HLTdecision'] ==1)    
    if (sampletype == 'photon'):
        cut = photoncut
    else:
        cut = halocut
    dfs = dfs.loc[cut]
    dfs = dfs.reset_index(drop=True)
#    print ("columns in dataframe:",dfs.columns)
    nevents = dfs.shape[0]
    print ("total events: ",nevents)
    if (nevents==0): continue
    for nev in tqdm(range(nevents)):
    #    for nev in range(10):
        if (nev%1000 == 0): print (nev,":events processed") 
        if (dfs['eerhE'][nev].size == 0):  continue
        listarr.append(make_nparray(nev))



if (len(listarr) > 0):
    print ("starting saving arrays")
    save_arr(listarr,savename+'.npy')
#    print ("starting plotting arrays")
#    plot_images(listarr,savename+".pdf",len(listarr) if (len(listarr) < 10) else 10)
else:
    print ("no event selected, skipping saving")
