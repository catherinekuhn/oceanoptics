# this file extracts AC9 data written by Julia Oelker (ask, if anything is unclear, not working or wrong)

import matplotlib.pyplot as plt
import numpy as np
import os 

outdir="Lab3/"
# put all raw data files in indir and all spectra for these files will be corrected at once
indir=outdir+"raw_data/"

#a_file="ac9_274_20170711191457_a_cal1_T23_5.dat"
#a_file="ac9_274_20170711191457_a_cdom_T24.dat"
#a_file="ac9_274_20170711191457_c_cal1_T24_1.dat"
#a_file="ac9_274_20170711191457_c_cdom_T24_1.dat"



for a_file in os.listdir(indir):

    ac9_file=open(indir+"/"+a_file)

    sal=30.0
    tmp_file=a_file[-7:-4]  #######CAUTION: this assumes a file format where the last three characters are temperature##########
    tmp=float(tmp_file[:2]+"."+tmp_file[2:3])

    caltmp=19.2
    lambda_list=[412,440,488,510,532,555,650,676,715]

    phitsw=[0.0003 ,.0002, .0001, .0003 ,.0001, .0002, -.0001 ,-.0001, .0027]
    phisa=[.00018, .00008, .00008, .00009 ,.00004 ,.00008, .0001 ,.00007 ,-.00018]
    phisc=[.00007, -.00007, -.00007 ,-.00007, -.00008, -.00008, -.00005, -.00007, -.00032]

    j=0
    def tempcorr(adat, cdat):
        adat_tcorr=[]
        cdat_tcorr=[]
        
        count_wave=0
        for wave in lambda_list:
            acorr=adat[count_wave]-(phitsw[count_wave]*(tmp-caltmp)+phisa[count_wave]*sal)
            ccorr=cdat[count_wave]-(phitsw[count_wave]*(tmp-caltmp)+phisc[count_wave]*sal)
            adat_tcorr.append(acorr)
            cdat_tcorr.append(ccorr)
            count_wave+=1

        return adat_tcorr, cdat_tcorr

    data=np.genfromtxt(ac9_file, 'float', skip_header=31)



    print data
    print np.shape(data)


    a_sorting={412:7,440:8, 488:9, 510:13, 532:14, 555:15, 650:1, 676:2, 715:3}
    c_sorting={412:16,440:17, 488:18, 510:4, 532:5, 555:6, 650:10, 676:11, 715:12}

    a_data=[]
    c_data=[]
    for wavelength in lambda_list:
        median_a=np.median(data[:,a_sorting[wavelength]])
        median_c=np.median(data[:,c_sorting[wavelength]])

        a_data.append(median_a)
        c_data.append(median_c)
     
    data_tcorr=tempcorr(a_data, c_data)
    a_data_tcorr=data_tcorr[0]
    c_data_tcorr=data_tcorr[1]
    print data_tcorr

       
    a_complete=np.concatenate([(np.array(lambda_list))[:,np.newaxis], (np.array(a_data))[:,np.newaxis]], axis=1)
    c_complete=np.concatenate([(np.array(lambda_list))[:,np.newaxis], (np.array(c_data))[:,np.newaxis]], axis=1)

    a_tcorr_complete=np.concatenate([(np.array(lambda_list))[:,np.newaxis], (np.array(a_data_tcorr))[:,np.newaxis]], axis=1)
    c_tcorr_complete=np.concatenate([(np.array(lambda_list))[:,np.newaxis], (np.array(c_data_tcorr))[:,np.newaxis]], axis=1)


    np.savetxt(outdir+"a_data_"+a_file, a_complete)
    np.savetxt(outdir+"c_data_"+a_file, a_complete)
    np.savetxt(outdir+"a_corr_"+a_file, a_tcorr_complete)
    np.savetxt(outdir+"c_corr_"+a_file, c_tcorr_complete)
    print a_data
    print np.shape(a_data)   

    if "aside" in a_file:
        plt.plot(lambda_list, a_data, label='raw_data')
    elif "cside" in a_file:
        plt.plot(lambda_list, c_data, label='c')
    if "aside" in a_file:
        plt.plot(lambda_list, a_data_tcorr, label='t corrected')
    elif "cside" in a_file:
        plt.plot(lambda_list, c_data_tcorr, label='c_tcorr')

    plt.legend()
    plt.ylabel("absorption / 1/m")
    plt.xlabel("wavelength / nm")
    plt.savefig(outdir+"raw_and_tcorr_"+a_file[:-4]+".png")
    plt.savefig(outdir+"raw_and_tcorr_"+a_file[:-4]+".eps")

    plt.show()


