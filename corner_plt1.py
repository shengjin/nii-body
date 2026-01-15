
import numpy as np
import corner


font = {'family' : 'serif', #monospace
'weight' : 'bold', #bold
'size'   : 10,
}


n_para = 7
ndim = n_para
skipnum = 300000 #00000

fname_ch = "%s%s" % ("/home/planet/Desktop/nii-c/src/results082900/", "chain3.dat")
print(fname_ch)
# read chain 
# number of n_para
i_ch_tup = np.genfromtxt(fname_ch, skip_header=skipnum, usecols=(0, 1, 2, 3, 4, 5, 6))
#print(i_ch_tup.shape)
#### change to array
i_ch = np.asarray(i_ch_tup)

m = i_ch.shape[0]
n = i_ch.shape[1]
print(m,n)


    


figure = corner.corner(i_ch,
        labels = [
    r"$p$",
    r"$e$",
    r"$ci$",
    r"$w$",
    r"$omg$",
    r"$cm0$",
    r"$mass$"
   
            ],
        label_kwargs={"fontsize": 7},
        quantiles=[0.16, 0.5, 0.84],
        show_titles=True,
        title_kwargs={"fontsize": 7}, #smooth=True
        ) 

figure.savefig("pic_corner.png", dpi=200)


