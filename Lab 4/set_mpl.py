from matplotlib import rcdefaults
from matplotlib import rcParams

def set_mpl(label_size = 16, ticklabel_size = 14):
    """Uniform plot settings."""
    
    rcdefaults()
    rcParams["text.usetex"]= True
    rcParams["font.family"]= "sans-serif"
    rcParams["font.serif"]= []                  
    rcParams["font.sans-serif"]= []
    rcParams["font.monospace"]= []
    rcParams["axes.labelsize"]= label_size               
    rcParams["font.size"]= label_size
    rcParams["axes.titlesize"]= label_size               
    rcParams["legend.fontsize"]= ticklabel_size               
    rcParams["xtick.labelsize"]= ticklabel_size
    rcParams["ytick.labelsize"]= ticklabel_size
    rcParams["legend.fontsize"]= ticklabel_size
    rcParams["axes.linewidth"]  = 2 
    rcParams["xtick.major.size"]  = 12
    rcParams["ytick.major.size"]  = 12 
    rcParams["xtick.major.width"]  = 1.5
    rcParams["ytick.major.width"]  = 1.5

