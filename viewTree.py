from __future__ import print_function
from __future__ import division

import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.mplot3d import Axes3D
from ROOT import TFile, TTree, TLorentzVector, TVector3, TRandom3

import numpy as np

# Load tree from file
inName = "4l_fullSM.root"
inFile = TFile(inName, "READ")
t      = inFile.Get("tree")


# Initialize RNG
rng = TRandom3()
rng.SetSeed()
nEntries = t.GetEntries()


# Create charge and flavor dictionaries
markerDict = {-1:'_',  +1:'+', 0:''}
colorDict  = {11:'#a2142f',  13:'#0072bd', 0:'#7e2f8e'}   # red = electron; blue = muon
angleDict  = {11:'#e9506d',  13:'#3eb3ff'}
lineDict   = {'p':'-', 'k':'--'} # solid = high-M; dashed = low-M
flavorDict = {11:r"$\vec{p}_{\mathrm{e}^+}$", -11:r"$\vec{p}_{\mathrm{e}^-}$", \
                      13:r"$\vec{p}_{\mu^+}$", -13:r"$\vec{p}_{\mu^-}$"}
chargeDict = {-1:'-',  +1:'+'}

rcParams['xtick.color'] = (0.75,0.75,0.75,0)
rcParams['ytick.color'] = (0.75,0.75,0.75,0)

# Create axes
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.autoscale(enable=False)
ax.xaxis.set_pane_color((1,1,1))
ax.yaxis.set_pane_color((1,1,1))
ax.zaxis.set_pane_color((1,1,1))
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_zticklabels([])
ax.w_xaxis.line.set_color((1,1,1,0))
ax.w_yaxis.line.set_color((1,1,1,0))
ax.w_zaxis.line.set_color((1,1,1,0))
ax.xaxis._axinfo['tick']['inward_factor'] = 0
ax.xaxis._axinfo['tick']['outward_factor'] = 0
ax.yaxis._axinfo['tick']['inward_factor'] = 0
ax.yaxis._axinfo['tick']['outward_factor'] = 0
ax.zaxis._axinfo['tick']['inward_factor'] = 0
ax.zaxis._axinfo['tick']['outward_factor'] = 0
ax.zaxis._axinfo['tick']['inward_factor'] = 0
#ax.xaxis.pane.fill = False

plt.rc('text', usetex=True)

rankDict = {1:r"$\vec{p}_{\ell_1}$", 2:r"", 3:r"", 4:r"", 100:r"$\vec{p}_{\ell_{2,3,4}}$"}
#rankDict = {1:r"", 2:r"$\vec{p}_{\mathrm{e}_2}$", 3:r"", 4:r"", 100:r"$\vec{p}_{\mathrm{Z}_2}$"}


# Get a random event
evtNum = 98495 #rng.Integer(nEntries)
t.GetEntry(evtNum)


# Get leptons
trio = t.k_plus + t.k_minus + t.p_minus
kk_pair, pp_pair = t.k_plus + t.k_minus, t.p_plus + t.p_minus
p4   = (t.p_plus,    t.p_minus,   t.k_plus,   t.k_minus)#kk_pair)#, trio)   # Momentum
q    = (+1,          -1,          +1,         -1, 0)          # Charge
pdg  = (t.p_pdg,     t.p_pdg,     t.k_pdg,    t.k_pdg)#t.k_pdg)#,0)     # Flavor
pair = ('p',         'p',         'p',        'p', 'k')
rank = (t.p_plus_i,  t.p_minus_i, t.k_plus_i, t.k_minus_i, 99) # Momentum order

st = (1.1, 1.5, 1.7, 1.7, 1.3)#1.4)#1.3)

# Draw lines and +/- markers
for _p4, _q, _pdg, _pair, _st, _r in zip(p4, q, pdg, pair, st, rank):
    ax.quiver(0, 0, 0, _p4.Px(), _p4.Py(), _p4.Pz(), pivot='tail', length=_p4.P(), color=colorDict[_pdg], \
            arrow_length_ratio=4/_p4.P(), linestyle=lineDict[_pair])
    ax.text(_st*_p4.Px(), _st*_p4.Py(), _st*_p4.Pz(), flavorDict[_pdg*_q], None, fontsize=21)
    #ax.text(_st*_p4.Px(), _st*_p4.Py(), _st*_p4.Pz(), rankDict[_r+1], None, fontsize=21)

name = "untitled"
#name = "zframe"

#name = "kinematics"


'''
# beta
p_2, kk_pair = t.p_minus.Vect(), kk_pair.Vect()
n = p_2.Cross(kk_pair).Unit()
u = 10 * p_2.Unit()

angle, dangle = 0, 0.01
dx, dy, dz = [], [], []
while angle < p_2.Angle(kk_pair):
    dx.append(u.X())
    dy.append(u.Y())
    dz.append(u.Z())
    u.Rotate(dangle, n)
    angle = angle + dangle

ax.plot(dx, dy, dz, color="#bf6ecf")
l = 20 * p_2.Unit()
l.Rotate(3*p_2.Angle(kk_pair)/4, n)
ax.text(l.X(), l.Y(), l.Z()+5, r"$\beta$", None, fontsize=21)
name = "beta"
'''

'''
# alpha_p
p_plus, p_minus = t.p_plus.Vect(), t.p_minus.Vect()
n = p_plus.Cross(p_minus).Unit()
u = 10 * p_plus.Unit()

angle, dangle = 0, 0.01
dx, dy, dz = [], [], []
while angle < t.pp_angle:
        dx.append(u.X())
        dy.append(u.Y())
        dz.append(u.Z())
        u.Rotate(dangle, n)
        angle = angle + dangle
                            
ax.plot(dx, dy, dz, color=angleDict[t.p_pdg])
l = 12 * p_plus.Unit()
l.Rotate(3*t.pp_angle/4, n)
ax.text(l.X(), l.Y(), l.Z()+5, r"$\alpha_{\mathrm{Z}_1}$", None, fontsize=21)


# alpha_k
k_plus, k_minus = t.k_plus.Vect(), t.k_minus.Vect()
n = k_plus.Cross(k_minus).Unit()
u = 10 * k_plus.Unit()

angle = 0
dx, dy, dz = [], [], []
while angle < t.kk_angle:
    dx.append(u.X())
    dy.append(u.Y())
    dz.append(u.Z())
    u.Rotate(dangle, n)
    angle = angle + dangle

ax.plot(dx, dy, dz, color=angleDict[t.k_pdg])
l = 30 * k_plus.Unit()
l.Rotate(t.kk_angle/3, n)
ax.text(l.X(), l.Y(), l.Z(), r"$\alpha_{\mathrm{Z}_2}$", None, fontsize=21)

name = "alpha"
'''


# phi
pp_pair, kk_pair = 0.25 * pp_pair.Vect(), 0.25 * kk_pair.Vect()
p_plus, p_minus = t.p_plus.Vect(), t.p_minus.Vect()
k_plus, k_minus = t.k_plus.Vect(), t.k_minus.Vect()
n_pp, n_kk = p_plus.Cross(p_minus).Unit(), 25 * k_plus.Cross(k_minus).Unit()


kk_cross = n_kk.Cross(kk_pair).Unit()
pp_cross = n_pp.Cross(pp_pair).Unit()

#u = 5 * pp_cross.Unit()
u = 10 * pp_pair.Unit()
angle, dangle = 0, 0.01
dx, dy, dz = [], [], []
while np.abs(angle) < t.phi:
    dx.append(u.X())
    dy.append(u.Y())
    dz.append(u.Z())
#   u.Rotate(dangle, kk_pair)
    u.Rotate(-dangle, pp_cross)
    angle = angle - dangle

ax.plot(dx, dy, dz, color="#77ac30")
#ax.quiver(dx[-10], dy[-10], dz[-10], dx[-1]-dx[-10], dy[-1]-dy[-10], dz[-1]-dz[-10], pivot='tail', length=1, \
#                  color="#77ac30", arrow_length_ratio=4)
ax.quiver(0, 0, 0, 10*pp_pair.Unit().X(), 10*pp_pair.Unit().Y(), 10*pp_pair.Unit().Z(), 
                pivot='tail', length=10, color="#77ac30", arrow_length_ratio=0, linestyle='--')
pp_pair.Rotate(-t.phi, pp_cross)
ax.quiver(0, 0, 0, 10*pp_pair.Unit().X(), 10*pp_pair.Unit().Y(), 10*pp_pair.Unit().Z(), 
                pivot='tail', length=10, color="#77ac30", arrow_length_ratio=0, linestyle='--')


#X = [0, t.p_minus.Px(), t.p_plus.Px(), -t.p_minus.Px(), -t.p_plus.Px()]
#Y = [0, t.p_minus.Py(), t.p_plus.Py(), -t.p_minus.Py(), -t.p_plus.Py()]
#Z = [0, t.p_minus.Pz(), t.p_plus.Pz(), -t.p_minus.Pz(), -t.p_plus.Pz()]

X = [0, t.p_minus.Px(), t.p_plus.Px(), t.p_plus.Px() - t.p_minus.Px()]
Y = [0, t.p_minus.Py(), t.p_plus.Py(), t.p_plus.Py() - t.p_minus.Py()]
Z = [0, t.p_minus.Pz(), t.p_plus.Pz(), t.p_plus.Pz() - t.p_minus.Pz()]
ax.plot_trisurf(X, Y, Z, color=colorDict[t.p_pdg], edgecolor='none', alpha=0.1)


#X = [0, 2*t.k_minus.Px(), 2*t.k_plus.Px(), -2*t.k_minus.Px(), -2*t.k_plus.Px()]
#Y = [0, 2*t.k_minus.Py(), 2*t.k_plus.Py(), -2*t.k_minus.Py(), -2*t.k_plus.Py()]
#Z = [0, 2*t.k_minus.Pz(), 2*t.k_plus.Pz(), -2*t.k_minus.Pz(), -2*t.k_plus.Pz()]

X = [t.k_minus.Px() - t.k_plus.Px(), t.k_minus.Px(), t.k_plus.Px(), 0]
Y = [t.k_minus.Py() - t.k_plus.Py(), t.k_minus.Py(), t.k_plus.Py(), 0]
Z = [t.k_minus.Pz() - t.k_plus.Pz(), t.k_minus.Pz(), t.k_plus.Pz(), 0]

X = [4*(t.k_minus.Px() - t.k_plus.Px()), 4*t.k_minus.Px() - 3*t.k_plus.Px(), t.k_plus.Px(), 0]
Y = [4*(t.k_minus.Py() - t.k_plus.Py()), 4*t.k_minus.Py() - 3*t.k_plus.Py(), t.k_plus.Py(), 0]
Z = [4*(t.k_minus.Pz() - t.k_plus.Pz()), 4*t.k_minus.Pz() - 3*t.k_plus.Pz(), t.k_plus.Pz(), 0]
ax.plot_trisurf(X, Y, Z, color=colorDict[t.k_pdg], edgecolor='none', alpha=0.1)


ax.plot(dx, dy, dz, color="#77ac30")

l = 15 * pp_pair.Unit()
l.Rotate(t.phi*.5, pp_cross)
ax.text(l.X(), l.Y(), l.Z(), r"$\phi$", None, fontsize=21)

name = "phi"



lim=25
ax.set_xlim(left=-lim, right=lim)
ax.set_ylim(bottom=-lim, top=lim)
ax.set_zlim(bottom=-lim, top=lim)

plt.savefig(name + ".pdf", transparent=True)
