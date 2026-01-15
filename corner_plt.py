
import numpy as np
import corner


font = {'family' : 'serif', #monospace
'weight' : 'bold', #bold
'size'   : 10,
}


n_para = 15
ndim = n_para
skipnum = 700000 #00000

fname_ch = "%s%s" % ("./results082900/", "chain7.dat")
print(fname_ch)
# read chain 
# number of n_para
i_ch_tup = np.genfromtxt(fname_ch, skip_header=skipnum, usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14))
#print(i_ch_tup.shape)
#### change to array
i_ch = np.asarray(i_ch_tup)

m = i_ch.shape[0]
n = i_ch.shape[1]
print(m,n)


figure = corner.corner(i_ch,
        labels = [
    r"$Tcont$",
    r"$Tdust$",
    r"$log(D8)$",
    r"$log(D9)$",
    r"$log(D10)$",
    r"$log(D11)$",
    r"$log(D12)$",
    r"$Tcont$",
    r"$Tdust$",
    r"$log(D8)$",
    r"$log(D9)$",
    r"$log(D10)$",
    r"$log(D11)$",
    r"$log(D11)$",
    r"$log(D12)$"
            ],
        label_kwargs={"fontsize": 14},
        quantiles=[0.16, 0.5, 0.84],
        show_titles=True,
        title_kwargs={"fontsize": 14}, #smooth=True
        ) 

figure.savefig("pic_corner.png", dpi=400)


