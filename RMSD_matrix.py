#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import mdtraj as md
import seaborn as sns # can be removed- I only use it for colormap


### set plot properties#####
plt.rcParams['svg.fonttype'] = 'none'
mpl.rcParams['font.sans-serif']="Arial"
plt.rcParams['figure.figsize'] = 3.5, 3  #width 1/3 of 7in and length 1/4 of 9.5in
params = {'legend.fontsize': 6,
          'legend.handlelength': 2}

cm=sns.color_palette("flare", as_cmap=True)

plt.rcParams['legend.frameon'] = True
plt.rcParams['legend.fancybox'] = True
mpl.rcParams['axes.linewidth'] = 0.3

mpl.rcParams['axes.labelsize'] = 7
mpl.rcParams['xtick.labelsize'] = 5
mpl.rcParams['ytick.labelsize'] = 5
mpl.rcParams['xtick.major.width'] = 0.3
mpl.rcParams['ytick.major.width'] = 0.3
mpl.rcParams['xtick.minor.width'] = 0.3
mpl.rcParams['ytick.minor.width'] = 0.3
mpl.rcParams['xtick.major.size'] = 2
mpl.rcParams['ytick.major.size'] = 2
mpl.rcParams['xtick.minor.size'] = 1 # half of the major ticks length
mpl.rcParams['ytick.minor.size'] = 1

plt.rcParams.update(params)

colors={"FF99SB-disp":"#AA3377","FF14SB":"#4477AA"} 

def main():
    
    
    #################### input parameters and data paths ###############

    ff1 = "FF14SB" 
    ff2 = "FF99SB-disp"
    
    base_path = 'C:/user/folder/'

    pdbPath = base_path+'topology.pdb'    
    
    trajPathList = [base_path+f'{ff1}_trajectory.dcd', base_path+f'{ff2}_trajectory.dcd'] 
    
    weightPathList =  [base_path+f'{ff1}_weights.dat', base_path+f'{ff2}_weights.dat']  
     

        
    cutoff = 50  # cutoff for RMSD comparison. here: top 50 most populated ensemble members will be compared
    
    
    #####################################################################
    
    
    Wpost1 = np.loadtxt(weightPathList[0])
    assert np.all(Wpost1[:,0]==np.arange(Wpost1.shape[0]))
    Wpost1 = Wpost1[Wpost1[:,1].argsort()[::-1]] 
    Wpost1_ = Wpost1[:cutoff,1]
    
    
    Wpost2 = np.loadtxt(weightPathList[1])
    assert np.all(Wpost2[:,0]==np.arange(Wpost2.shape[0]))
    Wpost2 = Wpost2[Wpost2[:,1].argsort()[::-1]] 
    Wpost2_ = Wpost2[:cutoff,1]    
   
   
    traj1 = md.load(trajPathList[0], top=pdbPath)               
    traj1_sorted = traj1[Wpost1[:cutoff,0].astype(int)]
   

    traj2 = md.load(trajPathList[1], top=pdbPath)               
    traj2_sorted = traj2[Wpost2[:cutoff,0].astype(int)]

    
        
    compare_ensembles(ff1, ff2, traj1_sorted, traj2_sorted, Wpost1_, Wpost2_)
        
    
    
def compare_ensembles(ffName1,ffName2,t1, t2, w1, w2):
    idxs = t1.topology.select('name CA')
    rmsd2d = np.zeros((t1.n_frames, t2.n_frames))
    for cl1 in range(t1.n_frames):
        rmsd2d[cl1]=md.rmsd(t2, t1, cl1, idxs)*10.0 #nm to A



    fig = plt.figure() 
    gs = fig.add_gridspec(2, 2,  height_ratios=(1,4), width_ratios=(1,4), hspace=0.08, wspace=0.08) 
       

    ax_main = fig.add_subplot(gs[1,1])  # main panel containing RMSD data 
    im=ax_main.imshow(rmsd2d, interpolation='none',cmap=cm, vmin=0, vmax=40)
    ax_main.xaxis.set_major_locator(MultipleLocator(5))
    ax_main.yaxis.set_major_locator(MultipleLocator(5))    
    ax_main.set_xlabel('cluster number')
    ax_main.tick_params(axis="y", labelleft=False)
    ax_main.set_aspect('auto')     
     
   
    ax_t = fig.add_subplot(gs[0,1], sharex=ax_main) # top panel
    ax_t.bar(np.arange(t2.n_frames),w2,color=colors[f'{ffName2}'], alpha=0.5, width=0.9)
    ax_t.tick_params(axis="x", labelbottom=False)
    ax_t.set_ylim(0,0.03)
    ax_t.yaxis.set_major_locator(MultipleLocator(0.02))
    ax_t.set_ylabel('weight')
    ax_t.text(13,0.022,'top 50 clusters in ensemble2', fontsize=5)
    
    
    ax_l = fig.add_subplot(gs[1,0], sharey=ax_main) # left panel
    ax_l.barh(np.arange(t1.n_frames),w1,color=colors[f'{ffName1}'], height=0.9)    
    ax_l.invert_xaxis()
    ax_l.set_xlim(0.03,None)
    ax_l.xaxis.set_major_locator(MultipleLocator(0.02))
    ax_l.set_xlabel('weight')
    ax_l.text(0.025,35,'top 50 clusters in ensemble1', rotation=90, fontsize=5)
    


    fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8)
    cb_ax = fig.add_axes([0.83, 0.1, 0.04, 0.616]) #color bar axis
    plt.colorbar(im, cax=cb_ax, orientation="vertical",fraction=0.046, pad=0.04).set_label(label=r'C$_{\alpha}$-RMSD [Å]',size=7)

    
    plt.savefig(f'RMSD_{ffName1}_vs_{ffName2}_dpi300.png', dpi=300, transparent=True)
    plt.savefig(f'RMSD_{ffName1}_vs_{ffName2}.svg', transparent=True)
    plt.close(fig)
  
  
    

if __name__ == '__main__':
    main()
