import time as tm
import numpy as Num
import sys
import csv
from optparse import OptionParser
from collections import Counter


start_time = tm.time()

class candidate:
	def __init__(self, DM, sigma, time, bin, downfact):
		self.DM = DM
		self.time = time
		self.sigma = sigma
		self.bin = bin
		self.downfact = downfact
	def __str__(self):
		return "%7.2f %7.2f %13.6f %10d   %3d\n"%\
			(self.DM, self.sigma, self.time, self.bin, self.downfact)
	def __cmp__(self, other):
		return cmp(self.bin, other.bin)            # Sort by time (i.e bin) by default

def cmpii(a, b):
    if (a > b):
        return 1
    elif (a == b):
        return 0
    else:
        return -1
    
def sk_true(DM,SNR,time,width,ld,ud,ls,us,lt,ut,lw,uw,cs):
    D = []
    T = []
    C = []
    W = []
    N = []
    S = []
    D_ = []
    T_ = []
    C_ = []
    W_ = []
    N_ = []
    S_ = []

    for jl, lk in enumerate(DM):
        if (ld <= DM[jl] <= ud and ls <= SNR[jl] <= us and lt <= time[jl] <= ut and lw <= width[jl] <= uw):
            D_.append(DM[jl])
            C_.append(SNR[jl])
            T_.append(time[jl])
            W_.append(width[jl])
            S_.append(4)
            N_.append(counts[Num.where(unique==lk)])			
            
    l = sorted(zip(C_, D_, T_, W_, N_), reverse = True)

    for i,ii in enumerate(l):
        D.append(ii[1])
        C.append(ii[0])
        W.append(ii[3])
        T.append(ii[2])
        N.append(ii[4])
        S.append(4)
        
    del N_,D_,C_,T_,W_,S_,l,DM,time,SNR,width        
        
    D = Num.asarray(D)
    C = Num.asarray(C)
    T = Num.asarray(T)
    W = Num.asarray(W)*dt
    N = Num.asarray(N)
    N=Num.squeeze(N)
    return D,C,W,T,N
                    
# --------------All Functions are Above --------------------------------------------

usage = "usage: %prog [options] .dat files _or_ .singlepulse files"

parser = OptionParser(usage)

parser.add_option("-l", "--tsample", type="float", dest="timesample", default= 0.32768,
		help="Interval of each time bin in millisecond. Default: %default")
parser.add_option("-n", "--nfile", type="string", dest="filename",
		help="Name of the file")
parser.add_option("--bd", "--badfile", type="string", dest="debug",
		help="Writes files with bad data")
parser.add_option("--lt", "--lowertime", dest="ltime", type="float",
		default= 0.0, help="Lower Limit of the Time. Default: %default")
parser.add_option("--ut", "--uppertime", dest="utime", type="float",
		default= 10000.0, help="Upper Limit of the Time in seconds. Default: %default")
parser.add_option("--ld", "--lowerdm", dest="ldm", type="float",
		default= 0.0, help="Lower Limit of the DM. Default: %default")
parser.add_option("--ud", "--upperdm", dest="udm", type="float",
		default= 50000.0, help="Upper Limit of the DM. Default: %default")
parser.add_option("--ls", "--lowersnr", dest="lsnr", type="float",
		default= 0.0, help="Lower Limit of the SNR. Default: %default")
parser.add_option("--us", "--uppersnr", dest="usnr", type="float",
		default= 50000.0, help="Upper Limit of the SNR. Default: %default")
parser.add_option("--lw", "--lowerwidth", dest="lwidth", type="float",
		default= 0.0, help="Lower Limit of the Width in millisecond. Default: %default")
parser.add_option("--uw", "--upperwidth", dest="uwidth", type="float",
		default= 1.0, help="Upper Limit of the Width in seconds. Default: %default")
parser.add_option("--sm", "--smedian", dest="snrmedian", action = "store_true",
		default= False , help="Whether lower limit of snr is median or not. Default: %default")
parser.add_option("--sk", "--skip", dest="askip", type = "string",
		default= "False" , help="Whether to skip association or not. Default: %default")
parser.add_option("--ef", "--expfactor", dest="efactor", type="float",
		default= 3.0 , help="Exponential factor for variation in size of symbols as snr. Default: %default") # (addition 0)
parser.add_option("--cs", "--colorsaturation", dest="csaturation", type="float",
		help="Value of SNR after which color of symbols will not change. Default: %default") # (addition 0)
parser.add_option("--auto", "--automode", dest="auto", type="string",
		default= "True", help="Set to True to apply limits on DM: 0 to 2000. Default: %default")
parser.add_option("--dmstep", "--dmstepsize", dest="dmstep", type="float",
		default= "100.0", help="Use to change DM step size for --auto option. Default: %default")
# to get buffer count, source location and time information from C++ code ::
parser.add_option("--bc", "--buf_count", dest="buf_count", type="int",
		help="Buffer count of FRB detection")
parser.add_option("--mjd", "--mjd_time", dest="mjd", type="float",
		help="MJD of when the first sample")
parser.add_option("--ra", "--src_ra", dest="ra", type="float",
		help="RA of source")
parser.add_option("--dec", "--src_dec", dest="dec", type="float",
		help="Declanation of source")
parser.add_option("--buf_t", "--buf_t", dest="buf_t", type="float",
		help="Amount of telescope time uniquely processed in each buffer")
parser.add_option("--new_file_flag", "--new_file_flag", dest="new_file_flag", type="int",
		help="Check weather the radec.dat has changed")       #Added by Rushikesh




(opts, args) = parser.parse_args()
#parser.print_help()

#print "_______________________________________"
#print opts.askip
lt = opts.ltime
ut = opts.utime
ld = opts.ldm
ud = opts.udm
ls = opts.lsnr
us = opts.usnr
lw = opts.lwidth
uw = opts.uwidth
sm = opts.snrmedian
dt = opts.timesample
ef = opts.efactor # (addition 0) 
cs = opts.csaturation # (addition 0)
dmstep=opts.dmstep

buf_count = opts.buf_count
mjd = opts.mjd
ra = opts.ra
dec = opts.dec
buf_t = opts.buf_t
new_file_flag = opts.new_file_flag #Added by Rushikesh

t_time = float(buf_t)*int(buf_count)
print("telescope time in each buffer (unique): ", buf_t)
print("telescope time (upto previous buffer): ", t_time)

# converting ra and dec to degrees:
if dec < 0:
  dec_const = -1
else:
  dec_const = 1

print ("buf_count: ",buf_count," mjd: ",mjd," ra: ",ra, " dec: ", dec)

#ra_hh = int(ra/10000)
#ra_mm = int((ra - ra_hh*10000)/100)
#ra_sec = ra%100
#print("RA (HH:MM:SS)", ra_hh, ra_mm, ra_sec)

#dec_hh = int(abs(dec)/10000)
#dec_mm = int(abs(abs(dec) - dec_hh*10000)/100)
#dec_sec = abs(dec%100)
#print("DEC (HH:MM:SS)", dec_hh, dec_mm, dec_sec)

ra_deg = ra
dec_deg = dec
print("RA and DEC in degrees: ", ra_deg, dec_deg)

# Converting uw and lw from seconds to the scale of data; multiplication by 1000 to scale down from second to millisecond:
time_fact=1000./dt
uw=uw*time_fact       
lw=lw*time_fact
#print (uw)


maxdm=2000
mindm = 10

opts.askip = "True"

if opts.auto == "True" :
    ld_range = Num.array([mindm,maxdm])  #the list goes upto 2000 and not 2100!
else:
    ld=opts.ldm
    ud=opts.udm
    ld_range = Num.array([ld,ud])
    
print ("ls: ",ls," auto: ",opts.auto," ld range: ",ld_range)
#------------------------------------------------------------------------------

if opts.askip == "True":
    f = Num.fromfile(args[0], dtype = Num.float32)

    if args[0].endswith(".singlepulse"):  # For plotting
        filenmbase = args[0][args[0].rfind("/"):args[0].rfind(".singlepulse")]
    elif args[0].endswith(".dat"):
        filenmbase = args[0][args[0].rfind("/"):args[0].rfind(".dat")]
    else:
        filenmbase = args[0]
	
# edits end------------------------------------------------------------------------------	
    l = len(f)
    DM =[]
    time =[]
    SNR = []
    width = []
    for i in range(0,l,4):
        DM.append(f[i])
        time.append(f[i+1])
        SNR.append(f[i+2])
        width.append(f[i+3])
#---------------------------------------------------------------------------------
# Following are various  ways to find the unique and counts of DM values:
# For updated numpy, use the Num.unique option
#    unique, counts = Num.unique(DM, return_counts = True)			    
#    unique = Num.unique(DM)			
# Brute force option but it is rather slow    
#    counts = []
#    DM = Num.array(DM)
#    for vaal in unique:
#        ccs = Num.where(DM == vaal)
#        counts.append(len(ccs[0]))			
#    counts = Num.array(counts)
# Use the following method while working with older versions of Numpy:
    count_ter = Counter(DM)
    counts = Num.array([v for k,v in sorted(count_ter.iteritems())])
    unique = Num.array([k for k,v in sorted(count_ter.iteritems())])    
#---------------------------------------------------------------------------------
		
    if (cs == None):
        cs = 5.*Num.median(SNR)          # 5 is just a constant here chosen to scale the median for color saturation
#	print ("yeeeppp cs: ",cs)

    max_ld_loc = len(ld_range)-1

    for i in range(0,max_ld_loc):
        func_val = []
        grid_val = []
        snr_val = []
        dm_val = []
        t_val = []
        w_val = []
        d_val = []

        print ("step = ", i+1, "out of ", max_ld_loc)
        ld = ld_range[i]
        ud = ld_range[i+1]
        D, C, W, T, N = sk_true(DM,SNR,time,width,ld,ud,ls,us,lt,ut,lw,uw,cs) ##RD:Outlier removal and sorted output:
        
        rad_dm = 5.0
        func_thresh = 1.3
        grid_dm=Num.arange(ld,ud+5,5)  ##THe DM is taken in the steps of 5 need to follow the single pulse methodology?

        for val in grid_dm:
            loc = Num.squeeze(Num.where(abs(val-D) < rad_dm))
            if (loc.size != 0):
                dist=abs(val-D[loc])
                f_e = Num.exp(-(Num.square(dist))/(2.*(rad_dm**2)))
                f_sum = Num.sum(f_e)
                w_f_sum = Num.dot(f_e,(C[loc]**3))
                ratio_val = w_f_sum/f_sum
#                print val, loc.size,w_f_sum,f_sum, ratio_val
                func_val.append(ratio_val)
                grid_val.append(val)
                temp_index = Num.where(C == Num.max(C[loc]))
                d_val.append(D[temp_index][0])
                snr_val.append(C[temp_index][0])
                t_val.append(T[temp_index][0])
#                f_snr = Num.sum(C[loc])
                w_val.append(W[temp_index][0])
#                w_snr = Num.dot(W[loc],C[loc])
#                w_val.append(w_snr/f_snr)


        func_val = Num.array(func_val)
        grid_val = Num.array(grid_val)

#        func_val = func_val/Num.max(func_val)
        frb_loc = Num.squeeze(Num.where(func_val > func_thresh*Num.mean(func_val)))
        
        num_frb = Num.size(frb_loc)
                    
        snr_val = Num.array(snr_val)
        t_val = Num.array(t_val)
        snr_val = Num.array(snr_val)
        w_val = Num.array(w_val)
        d_val = Num.array(d_val)
                
        if num_frb > 1:
            for loc in range(0,num_frb):
                frnd_loc = Num.where(abs(grid_val-grid_val[frb_loc[loc]]) <= 4.*rad_dm) ##Whats this for? Wants to check weather any detection within +- 20 DM.
                max_frnd_loc = Num.where(snr_val[frnd_loc] != max(snr_val[frnd_loc]))
            found_frb = Num.delete(frnd_loc, max_frnd_loc[0])
            print ("#FRB found! DM: ", d_val[found_frb][0], " Time: ", t_val[found_frb][0], " Width: ", w_val[found_frb][0], " SNR: ", snr_val[found_frb][0])
            frb_flag = True
	    frb_time = t_time + t_val[found_frb][0]	
	    print("frb time:: ", frb_time)
# Rushikesh Edits -----
	    fname_frb = "FRB_detection_" + str(mjd) + "_" + str(ra_deg) + "_" + str(dec_deg) + "_" + "{:.3f}".format(frb_time) + ".csv"
            csvRow = [frb_time, t_val[found_frb][0], d_val[found_frb][0], snr_val[found_frb][0], w_val[found_frb][0]]
            fname_frb = "FRB_detection_" + str(mjd) + "_" + str(ra_deg) + "_" + str(dec_deg) + "_" + ".csv"
            print("Name of FRB detection file: ", fname_frb)
       	    csvfile = fname_frb
	    if buf_count == 0 or buf_count%4 == 0 or new_file_flag == 1:
           	    fname_frb = "FRB_detection_" + str(mjd) + "_" + str(ra_deg) + "_" + str(dec_deg) + "_" +"bc_start_"+ str(buf_count) + "_" + ".csv"
            	    print("Name of FRB detection file: ", fname_frb)
       	    	    csvfile = fname_frb
            	    with open(csvfile, 'w') as fp:
    	        	writer = csv.writer(fp, delimiter='\t',lineterminator='\n')
	    		writer.writerow(csvRow)
                    fp.close()
	    else:   
                    fname_frb = "FRB_detection_" + str(mjd) + "_" + str(ra_deg) + "_" + str(dec_deg) + "_" +"bc_start_"+ str((buf_count)-(buf_count)%4) + "_" + ".csv"
            	    print("Name of FRB detection file: ", fname_frb)
 		    csvfile = fname_frb     #RD:This and the above filename should be exactly same
            	    with open(csvfile, 'a') as fp:
    	        	writer = csv.writer(fp, delimiter='\t',lineterminator='\n')
	    		writer.writerow(csvRow) 
                    fp.close()           
	    
            
# edits end ------
        elif num_frb == 1:
            print ("*FRB found! DM: ", d_val[frb_loc], " Time: ", t_val[frb_loc], " Width: ", w_val[frb_loc], " SNR: ", snr_val[frb_loc])
            frb_flag = True

	    frb_time = t_time + t_val[frb_loc]	
	    print("frb time:: ", frb_time)
            csvRow = [frb_time, t_val[frb_loc], d_val[frb_loc], snr_val[frb_loc], w_val[frb_loc]]
            fname_frb = "FRB_detection_" + str(mjd) + "_" + str(ra_deg) + "_" + str(dec_deg) + "_" + ".csv"
            print("Name of FRB detection file: ", fname_frb)
       	    csvfile = fname_frb
	    if buf_count == 0 or buf_count%4 == 0 or new_file_flag == 1:
           	    fname_frb = "FRB_detection_" + str(mjd) + "_" + str(ra_deg) + "_" + str(dec_deg) + "_" +"bc_start_"+ str(buf_count) + "_" + ".csv"
            	    print("Name of FRB detection file: ", fname_frb)
       	    	    csvfile = fname_frb
            	    with open(csvfile, 'w') as fp:
    	        	writer = csv.writer(fp, delimiter='\t',lineterminator='\n')
	    		writer.writerow(csvRow)
                    fp.close()
	    else:   
                    fname_frb = "FRB_detection_" + str(mjd) + "_" + str(ra_deg) + "_" + str(dec_deg) + "_" +"bc_start_"+ str((buf_count)-(buf_count)%4) + "_" + ".csv"
            	    print("Name of FRB detection file: ", fname_frb)
 		    csvfile = fname_frb     ##This and the above filename should be same
            	    with open(csvfile, 'a') as fp:
    	        	writer = csv.writer(fp, delimiter='\t',lineterminator='\n')
	    		writer.writerow(csvRow)
                    fp.close()
            import time as tm
            print ("--- seconds --- ", tm.time() - start_time)
            sys.exit()	    

import time as tm
print ("--- seconds --- ", tm.time() - start_time)
