import yt 
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import pylab as plt
plt.rc('font', family='serif', size=20)
plt.rcParams['mathtext.fontset'] = "stix"
plt.ion()
# Load the dataset.

gamma=5.0/3


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

datasets=["vn1", "vc1", "p_n", "p_c", "b2"]
datasets=["vnc1", "b2", "vnc2", "p_nc"]

cls=["dodgerblue", "violet", "gray", "tab:pink"]
lss=["-", "--", "-.", ":"]
cls2=["red", "gold", "orange", "brown"]
lss2=["-", "--", "-.", ":"]
cls4=["cyan", "green", "lime", "yellow"]
lss4=["-", "--", "-.", ":"]


fig, ax1 = plt.subplots(nrows=len(datasets)//2, ncols=2, figsize=(8*2,6*len(datasets)/2))
if(len(datasets)==1):
  ax1=[ax1]

else:
  ax = []
  if(len(ax1.shape) == 1):
    for i in range(ax1.shape[0]):
        ax.append(ax1[i])
  else: 
    for i in range(ax1.shape[0]):
      for j in range(ax1.shape[1]):
        ax.append(ax1[i,j])
  ax1=ax


listFiles = [4]

dir1=["."]
ell=["2048p"]


ylim = {}
for dset in datasets:
  ylim[dset] = None
  

for k in listFiles: 


  for i in range(len(datasets)):
    dset=datasets[i]
    xx=[]
    arr =[]
    arr2=None
    arr4=None
    label1=None
    label2=None
    time = None
    if(dset in ["vnc1","vnc2", "p_nc"]):
      arr2=[]
    if(dset in ["p_nc"]):
      arr4=[]


    for d1 in dir1:
      ds = yt.load(join(d1,("shock%04d.dat" % k)), unit_system="mks")
      data = ds.all_data()
      xx1 = data["x"].to_ndarray() * 1e2 ##The multiplication is because of the normalization crap in yt 
      NN=len(xx)
      if(time is None): 
        time = float(ds["time"])
      xx.append(xx1)



      if(dset=="vn1"):
        arr.append(data["m_n1"].to_ndarray()/data["rho_n"].to_ndarray())  
      elif(dset=="vc1"):
        arr.append(data["m_c1"].to_ndarray()/data["rho_c"].to_ndarray())
      elif(dset=="vnc1"):
        arr2.append(data["m_n1"].to_ndarray()/data["rho_n"].to_ndarray())  
        arr.append(data["m_c1"].to_ndarray()/data["rho_c"].to_ndarray())
        label2=r"$vx_n$"
        label1=r"$vx_c$"
        
      elif(dset=="vn2"):
        arr.append(data["m_n2"].to_ndarray()/data["rho_n"].to_ndarray())  
      elif(dset=="vc2"):
        arr.append(data["m_c2"].to_ndarray()/data["rho_c"].to_ndarray())
      elif(dset=="vnc2"):
        arr2.append(data["m_n2"].to_ndarray()/data["rho_n"].to_ndarray())  
        arr.append(data["m_c2"].to_ndarray()/data["rho_c"].to_ndarray())
        label2=r"$vy_n$"
        label1=r"$vy_c$"
      elif(dset=="vn3"):
        arr.append(data["m_n3"].to_ndarray()/data["rho_n"].to_ndarray())  
      elif(dset=="vc3"):
        arr.append(data["m_c3"].to_ndarray()/data["rho_c"].to_ndarray())
      elif(dset=="p_n"):
        arr.append((gamma-1)*(data["p_n"].to_ndarray() -  0.5 * ((data["m_n1"].to_ndarray()**2 + data["m_n2"].to_ndarray()**2)/data["rho_n"]))) 
      elif(dset=="rho_n"):
        arr.append(data["rho_n"].to_ndarray())
      elif(dset=="p_c"):
        print(" FILE Energy = ",data["p_c"].to_ndarray())
        print("EK ", 0.5 * ((data["m_c1"].to_ndarray()**2 + data["m_c2"].to_ndarray()**2)/data["rho_c"])) 
        print("EM ", 0.5 * (data["b2"].to_ndarray()**2 + data["b1"].to_ndarray()**2))
        print("EINT ", (data["p_c"].to_ndarray() - 0.5 * (data["b2"].to_ndarray()**2 + data["b1"].to_ndarray()**2) -  0.5 * (
          (data["m_c1"].to_ndarray()**2 + data["m_c2"].to_ndarray()**2)/data["rho_c"]))) 
        arr.append((gamma-1)*(data["p_c"].to_ndarray() - 0.5 * (data["b2"].to_ndarray()**2 + data["b1"].to_ndarray()**2) -  0.5 * (
        (data["m_c1"].to_ndarray()**2 + data["m_c2"].to_ndarray()**2)/data["rho_c"]))) 
      elif(dset=="p_nc"):
        arr2.append((gamma-1)*(data["p_n"].to_ndarray() -  0.5 * ((data["m_n1"].to_ndarray()**2 + data["m_n2"].to_ndarray()**2)/data["rho_n"]))) 
        arr.append((gamma-1)*(data["p_c"].to_ndarray() - 0.5 * (data["b2"].to_ndarray()**2 + data["b1"].to_ndarray()**2) -  0.5 * (
        (data["m_c1"].to_ndarray()**2 + data["m_c2"].to_ndarray()**2)/data["rho_c"]))) 
        arr4.append((gamma-1)*(data["p_n"].to_ndarray() -  0.5 * ((data["m_n1"].to_ndarray()**2 + data["m_n2"].to_ndarray()**2)/data["rho_n"]))+
          (gamma-1)*(data["p_c"].to_ndarray() - 0.5 * (data["b2"].to_ndarray()**2 + data["b1"].to_ndarray()**2) -  0.5 * (
          (data["m_c1"].to_ndarray()**2 + data["m_c2"].to_ndarray()**2)/data["rho_c"])))
        label4=r"$p_n+p_c$"


        label2=r"$p_n$"
        label1=r"$p_c$"
      

      else:
        arr.append(data[dset].to_ndarray()) 
        if(dset=="b2"):
          label1=r"$B_y$"
        else:
          label1=dset

    ax=ax1[i]
    dset=datasets[i]
    ax.cla()
    #ax.grid(True)
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
    ax.tick_params(axis='both', labelsize=20)
    #ax.set_ylim((-1e-3,1e-3))
    #ax.set_ylabel(r"%s [%s]" % (labels[dset], units[dset]))
#    if(i<len(datasets)-1):
#      ax.set_xticklabels([])
#    else:
#      ax.set_xlabel("x/t")


    ax.set_xlabel("x/t")
    ax.set_xticks([0,0.2,0.4,0.6,0.8,1.0,1.2,1.4])

    print("TIME IS ", time)
    ax.yaxis.offsetText.set_fontsize(20) 
    #ax.plot(xx[:NP],arr[:NP],ls='-')
    #ax.plot(xx[ind],arr[ind],ls=lss[d],label=lls[d],lw=3,color=colors[d])
    xlim2=4e3
    for ll in range(len(arr)):
      ind = np.where(xx[ll]<=xlim2)
      ax.plot(xx[ll][ind]/time,arr[ll][ind],label=(r"%s,%s" % ((dset if label1 is None else label1), ell[ll])),lw=2.8,color=cls[ll], ls=lss[ll])
      if(not arr2 is None):
        ax.plot(xx[ll][ind]/time,arr2[ll][ind],label=(r"%s,%s"  %  ((dset if label2 is None else label2),ell[ll])),lw=2.8,color=cls2[ll], ls=lss2[ll])
      if(not arr4 is None):
        ax.plot(xx[ll][ind]/time,arr4[ll][ind],label=(r"%s,%s"  %  (label4,ell[ll])),lw=2.8,color=cls4[ll], ls=lss4[ll])
    if(i==3):
      ax.legend(fontsize=12,loc="upper center", ncol=2)
    else:
      ax.legend(fontsize=12)
    if(i==0):
      ax.annotate('', xy=(0.1, -0.042),  xycoords='data',
              xytext=(0.3, -0.042), textcoords='data',
              arrowprops=dict(facecolor='black', shrink=0.03),
              horizontalalignment='right', verticalalignment='bottom')



    #if(ylim[dset] is None):
    #  ylim[dset] = ax.get_ylim()
    #else:
    #  ax.set_ylim(ylim[dset])
    ax.set_xlim((0,1.5)) 
  fig.suptitle(r"$t=%.2f$" % (time))
  plt.subplots_adjust(hspace=0.27, left=0.15)  
  fig.canvas.draw()
  plt.draw()
  #plt.savefig(join(dir1,"shock%04d.png" % k ))
  plt.savefig("shock%04d.pdf" % k )

  
  
  plt.show()


