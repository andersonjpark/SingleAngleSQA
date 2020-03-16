import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from matplotlib.ticker import AutoMinorLocator
import matplotlib.colors as colors
import h5py
from scipy.interpolate import interp1d

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap
orig_cmap = plt.get_cmap('gist_heat')
cmap = truncate_colormap(orig_cmap, 0, 0.85, 200)

clight = 3e10 # cm/s

ns = 2
ng=25
every=100

mpl.rcParams['font.size'] = 22
mpl.rcParams['font.family'] = 'serif'
mpl.rc('text', usetex=True)
mpl.rcParams['xtick.major.size'] = 7
mpl.rcParams['xtick.major.width'] = 2
mpl.rcParams['xtick.major.pad'] = 8
mpl.rcParams['xtick.minor.size'] = 4
mpl.rcParams['xtick.minor.width'] = 2
mpl.rcParams['ytick.major.size'] = 7
mpl.rcParams['ytick.major.width'] = 2
mpl.rcParams['ytick.minor.size'] = 4
mpl.rcParams['ytick.minor.width'] = 2
mpl.rcParams['axes.linewidth'] = 2

E = 4.*(np.arange(0,25,1)+1)
noosc_color="green"

def line_color(ig,ng, cmap):
    return cmap(float(E[ig])/E[-1])

def magnitude(r,i):
    return np.sqrt(r**2 + i**2)

fig, axes = plt.subplots(3, 5, sharex=True, figsize=(16,8))

infile = h5py.File("/mnt/data/2019_IsotropicSQA/noosc/single/noosc_full_NE25/f.h5","r")
nt = np.array(infile["r(cm)"]) / clight
nf = np.array(infile["fmatrixf"])
nfz = .5 * (nf[:,:,:,0,0,0]-nf[:,:,:,1,1,0])
nfx = (nf[:,:,:,0,1,0])
nfy = (nf[:,:,:,0,1,1])
nfz0 = nfz[0,:,:]
nfx0 = nfx[0,:,:]
nl = np.sqrt(nfx**2 + nfy**2 + nfz**2)
nl0 = nl[0,:,:]
nfoffdiag = np.sqrt(nfx**2 + nfy**2)

infile = h5py.File("output.h5","r")
t = np.array(infile["r(cm)"]) / clight
f = np.array(infile["fmatrixf"])
fz = .5 * (f[:,:,:,0,0,0]-f[:,:,:,1,1,0])
fx = (f[:,:,:,0,1,0])
fy = (f[:,:,:,0,1,1])
fz0 = fz[0,:,:]
fx0 = fx[0,:,:]
l = np.sqrt(fx**2 + fy**2 + fz**2)
l0 = l[0,:,:]
foffdiag = np.sqrt(fx**2 + fy**2)
for ip in range(5):
    ig = ip*5+4
    nig = ig #*2
    color=line_color(ig,ng,cmap)
    if(ip==0):
        acolor = "0.5"
        aalpha = 1.0
    else:
        acolor=color
        aalpha=0.5
    p1, = axes[0][ip].semilogx(t,  fz[:,0,ig] /np.abs(fz0[0,ig]  ), color=color, linewidth=1.5)
    p2, = axes[0][ip].semilogx(t,  fz[:,1,ig] /np.abs(fz0[1,ig]  ), color=acolor, linewidth=1.5,alpha=aalpha)
    p3, = axes[0][ip].semilogx(nt,nfz[:,0,nig]/np.abs(nfz0[0,nig]), color=noosc_color, linewidth=2)
    p4, = axes[0][ip].semilogx(nt,nfz[:,1,nig]/np.abs(nfz0[1,nig]), color=noosc_color, linewidth=2,linestyle="--")
    
    axes[1][ip].semilogx(t,  foffdiag[:,0,ig] /np.abs(fz0[0,ig]  ), color=color, linewidth=1.5)
    axes[1][ip].semilogx(t,  foffdiag[:,1,ig] /np.abs(fz0[1,ig]  ), color=acolor, linewidth=1.5,alpha=aalpha)
    axes[1][ip].semilogx(nt,nfoffdiag[:,0,nig]/np.abs(nfz0[0,nig]), color=noosc_color, linewidth=2)
    axes[1][ip].semilogx(nt,nfoffdiag[:,1,nig]/np.abs(nfz0[1,nig]), color=noosc_color, linewidth=2,linestyle="--")
#/fx0[0,ig]  
#/fx0[1,ig]  
#/nfx0[0,nig]
#/nfx0[1,nig]

    axes[2][ip].semilogx(t,   l[:,0,ig] /np.abs(fz0[0,ig]), color=color, linewidth=1.5)
    axes[2][ip].semilogx(t,   l[:,1,ig] /np.abs(fz0[1,ig]), color=acolor, linewidth=1.5, alpha=aalpha)
    axes[2][ip].semilogx(nt, nl[:,0,nig]/np.abs(nfz0[0,nig]), color=noosc_color, linewidth=2)
    axes[2][ip].semilogx(nt, nl[:,1,nig]/np.abs(nfz0[1,nig]), color=noosc_color, linewidth=2, linestyle="--")

    axes[0][ip].text(4e-6,2.2, r"$%d\,\mathrm{MeV}$"%E[ig], color=color,fontsize=18,horizontalalignment='center',verticalalignment='center')
    if(ip==0):
        axes[1][0].legend(
            [p1,p2],
            [r"$\nu$",r"$\bar{\nu}$"],fontsize=18, ncol=1,loc="upper right", frameon=False, handletextpad=0)
        axes[1][1].legend(
            [p3,p4],
            [r"$\nu\,(\mathrm{no\,\,osc.})$",r"$\bar{\nu}\,(\mathrm{no\,\,osc.})$"],fontsize=18, ncol=1,loc="upper right", frameon=False, handletextpad=0)


# colorbar
#cax = fig.add_axes([.91, 0.1, 0.02, 0.8])
#norm = mpl.colors.Normalize(vmin=0, vmax=E[-1])
#cb1 = mpl.colorbar.ColorbarBase(cax, cmap=cmap,norm=norm,orientation='vertical')
#cb1.set_label(r'$h\nu\,(\mathrm{MeV})$')
#cax.tick_params(which='both',direction='in')

ylabelx = -.35
plt.text(ylabelx,0.5,r"$f_{(z)}(t)/L_\mathrm{eq}$",horizontalalignment='center',verticalalignment='center', transform=axes[0][0].transAxes,rotation=90,fontsize=22)
plt.text(ylabelx,0.5,r"$f_{(\perp)}(t)/L_\mathrm{eq}$",horizontalalignment='center',verticalalignment='center', transform=axes[1][0].transAxes,rotation=90,fontsize=22)
plt.text(ylabelx,0.5,r"$L(t)/L_\mathrm{eq}$",horizontalalignment='center',verticalalignment='center', transform=axes[2][0].transAxes,rotation=90,fontsize=22)

plt.subplots_adjust(wspace=0, hspace=0)
plt.minorticks_on()
axes[2][2].set_xlabel(r"$t\,(\mathrm{s})$",fontsize=22)
for ax in axes[:,1:].flatten():
    for label in ax.get_yticklabels():
        plt.setp(label,visible=False)
    #plt.setp(ax.get_xticklabels()[1], visible=False)
#axes[0].yaxis.set_major_locator(plt.MaxNLocator(5))
#axes[1].yaxis.set_major_locator(plt.MaxNLocator(6))
for ax in axes[0]:
    ax.set_ylim(-1.4,2.8)
    ax.yaxis.set_major_locator(plt.MultipleLocator(1))
for ax in axes[1]:
    ax.set_ylim(0,2.8)
    ax.yaxis.set_major_locator(plt.MultipleLocator(.5))
for ax in axes[2]:
    ax.set_ylim(.5,2.8)
    ax.yaxis.set_major_locator(plt.MultipleLocator(.5))
for ax in axes.flatten():
    #ax.grid()
    ax.set_xlim(3e-8,20e-6)
    #ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(2))
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.tick_params(which='both',direction='in')
    #ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    #ax.ticklabel_format(useOffset=False)
    
plt.savefig("f_osc.png",bbox_inches="tight")
