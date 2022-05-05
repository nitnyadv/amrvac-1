import yt 
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import pylab as plt
plt.rc('font', family='serif', size=20)
plt.rcParams['mathtext.fontset'] = "stix"
plt.ion()
# Load the dataset.




units={}
units["m1"]=r"$kg/m^2/s$"
units["rho1"]=r"$kg/m^3$"
units["rho"]=r"$kg/m^3$"
units["v1"]=r"$m/s$"
units["bx1"]=r"$T$"

labels={}
labels["m1"]=r"$\rho v_z$"
labels["rho1"]=r"$\rho_1$"
labels["rho"]=r"$\rho$"
labels["v1"]=r"$v_z$"
labels["bx1"]=r"$B_{\rm x1}$"

from os.path import join

#datasets=["vn1", "vc1", "b1"]
#datasets=["vn3", "vc3", "b3", "p_n", "p_c"]
datasets=["vn3", "vc3", "b3"]
fig, ax1 = plt.subplots(nrows=len(datasets), ncols=1, figsize=(8,6*len(datasets)))
if(len(datasets)==1):
  ax1=[ax1]

#listFiles = range(0,10000,10)
listFiles = range(0,10000)
#listFiles = [16]

dir1="."

alpha = 1e4
alpha = 1e-1

#dir1="alpha-1e4"
#alpha = 1e4

nn = 5
ampl = 1e-3

gamma = 5.0/3


pn0 = 2e1    
pc0 = 1e1 
rhon0 = 2e1   
rhoc0 = 1e1   
b0 = 0.4


def getAnFunc(x, time, alpha, nn, ampl):
  from math import pi,sqrt

  k = 2*pi*nn/(x[-1]-x[0])
  va0 = b0/sqrt(rhoc0)

  #vals2
  if(alpha == 1e-3):
    omega = 2.483586707670307 + 1j * 0.009999837878226345
  elif(alpha == 1e-1):
    omega = 1.8476978931291896 + 1j*0.7139437294531664
  elif(alpha == 1e0):
    omega = 1.435586475340743 + 1j* 0.06869586285089642
  elif(alpha == 1e3):
    omega = 1.4339343253916603  + 1j * 0.00006853892165121452
  elif(alpha == 1e4):
    omega = 1.4339343237700348 + 1j * 6.853892149619427e-6
  elif(alpha == 0e0):
    omega = k * va0

  Vc = ampl * va0
  Vn = Vc * (alpha * rhon0 *rhoc0)/(1j * rhon0 * omega + alpha * rhon0 *rhoc0)
  B1 = - k * Vc * b0 / omega


  print("Vn Vc B")
  print(Vn, Vc, B1)

  wave = np.exp(1j *  (omega * time - k*(x - x[0]))) 
  vx_n1 = np.real(Vn * wave)
  vx_c1 = np.real(Vc * wave)
  bx1 = np.real(B1 * wave)
  return (vx_n1,vx_c1,bx1)

ylim = {}
for dset in datasets:
  ylim[dset] = None
  

for k in listFiles: 
  ds = yt.load(join(dir1,("w2f%04d.dat" % k)), unit_system="mks")
  data = ds.all_data()

  xx = data["x"].to_ndarray() * 1e2 ##The multiplication is because of the normalization crap in yt 
  NN=len(xx)

  xlim2=2.1
  #xlim2=1.7
  ind = np.where(xx<=xlim2)

  time = float(ds["time"])
  vx_n1, vx_c1, bx1 = getAnFunc(xx, time, alpha, nn, ampl)

  for i in range(len(datasets)):
    ax=ax1[i]
    dset=datasets[i]
    ax.cla()
    ax.grid(True)
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
    ax.tick_params(axis='both', labelsize=20)
    #ax.set_ylim((-1e-3,1e-3))
    #ax.set_ylabel(r"%s [%s]" % (labels[dset], units[dset]))
    if(i<len(datasets)-1):
      ax.set_xticklabels([])
    else:
      ax.set_xlabel("z")

    if(dset=="vn3"):
      arr = data["m_n3"].to_ndarray()/(data["rho_n"].to_ndarray() + rhon0) 
      arr2 = vx_n1
    elif(dset=="vc3"):
      arr = data["m_c3"].to_ndarray()/(data["rho_c"].to_ndarray() + rhoc0)
      arr2 = vx_c1
    elif(dset=="b3"):
      arr = data["b3"].to_ndarray()
      arr2 = bx1
    elif(dset=="p_n"):
      arr = (gamma-1)*(data["p_n"].to_ndarray() -  0.5 * (data["m_n3"].to_ndarray()**2/(data["rho_n"]+rhon0))) 
      arr2 = np.zeros(arr.shape)
    elif(dset=="rho_n"):
      arr = data["rho_n"].to_ndarray() - rhon0
      arr2 = np.zeros(arr.shape)
    elif(dset=="p_c"):
      arr = (gamma-1)*(data["p_c"].to_ndarray() - 0.5 * (data["b3"].to_ndarray()**2 + data["b1"].to_ndarray()**2) -  0.5 * (data["m_c3"].to_ndarray()**2/(data["rho_c"] + rhoc0)))
      arr2 = np.zeros(arr.shape)
    elif(dset=="rho_c"):
      arr = data["rho_c"].to_ndarray() - rhoc0
      arr2 = np.zeros(arr.shape)
    else:
      print("NOT IMPL")
      import sys
      sys.exit(0)
    ax.yaxis.offsetText.set_fontsize(20) 
    #ax.plot(xx[:NP],arr[:NP],ls='-')
    #ax.plot(xx[ind],arr[ind],ls=lss[d],label=lls[d],lw=3,color=colors[d])
    ax.plot(xx[ind],arr[ind],label=dset,lw=3,color="black")
    ax.plot(xx[ind],arr2,label=dset + " AN",lw=3, color="red", ls="--")

    ax.legend(fontsize=10)

    if(ylim[dset] is None):
      ylim[dset] = ax.get_ylim()
    else:
      ax.set_ylim(ylim[dset])
    
  fig.suptitle(r"$\alpha=%.1e$, $t=%.2f$" % (alpha,time))
  plt.subplots_adjust(hspace=0.27, left=0.15)  
  fig.canvas.draw()
  plt.draw()
  plt.savefig(join(dir1,"w2f%04d.png" % k ))

  
  
  plt.show()


