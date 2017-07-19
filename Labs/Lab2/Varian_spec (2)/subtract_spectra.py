#this file subtracts an already TS corrected MilliQ reference from an already TS corrected sample measurement. Extraction and TS correction have to be previously done with rd_ac9.py. written by Julia Oelker (ask, if anything is unclear, not working or wrong)
import numpy as np
import matplotlib.pyplot as plt

#cal_file="ac9_274_20170711191457_a_cal1_T23_5.dat"
#sample_file="ac9_274_20170711191457_a_cdom_T24.dat"
#cal_file="ac9_274_20170711191457_c_cal1_T24_1.dat"
#sample_file="ac9_274_20170711191457_c_cdom_T24_1.dat"

outdir="Lab3/"
indir=outdir
#cal_file="ac974_g1b_diw_aside_230.dat"
cal_file="ac974_g1b_diw_cside_232.dat"

sample_file="ac974_g1b_cul2_cside_243.dat"

prop="c" ######CAUTION: has to be switched from a to c depending on the file ######
cal_tcorr=np.genfromtxt(indir+prop+"_corr_"+cal_file)
sample_tcorr=np.genfromtxt(indir+prop+"_corr_"+sample_file)

data_mqcorr=(sample_tcorr[:,1]-cal_tcorr[:,1])[:,np.newaxis]

data_mqcorr_complete=np.concatenate([(cal_tcorr[:,0])[:,np.newaxis], data_mqcorr], axis=1)


np.savetxt(outdir+prop+"_mq_t_corr"+sample_file, data_mqcorr_complete)
plt.plot(cal_tcorr[:,0], data_mqcorr, label ="MQ and t corrected")
plt.plot(cal_tcorr[:,0], sample_tcorr[:,1], label="t corrected")
plt.legend()
plt.ylabel("absorption / 1/m")
plt.xlabel("wavelength / nm")
plt.savefig(outdir+prop+"_mq_t_corr"+sample_file[:-4]+".png")
plt.savefig(outdir+prop+"_mq_t_corr"+sample_file[:-4]+".eps")
plt.show()


