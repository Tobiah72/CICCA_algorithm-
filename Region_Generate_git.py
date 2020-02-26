#######################################
#LAST MODIFIED BY Tobiah Steckel 07/2/19
#######################################
'''The purpose of this script is to read in a list of granules to inform
what files it should read in from a database of merged Geostationary (Meteosat or GOES)-
Active radar dataset(CloudSat). Once read in, relevant information from the given datasets
are sorted by a grouping methodology in order to generate "reference vectors"
defined by 1) geographical region and 2) CloudSat Profiling Radar (CPR) 
convective flag. The results are exported as a .txt file according to those criteria. '''

import numpy as np
import matplotlib
from numpy import *
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import glob
from numpy import *
from scipy.interpolate import interp1d
set_printoptions(threshold=nan)
import collections
import itertools
from numpy.ma import masked_where
from matplotlib.pyplot import figure, colorbar
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib.mlab as mlab
matplotlib.rcParams.update({'font.size': 50})
from matplotlib import ticker
from itertools import groupby
from operator import itemgetter
from collections import Counter
'''here we describe regional boundaries'''
Region = '1' ; Region_f = 'One'; lat_lo = 5; lat_hi = 15; lon_lo = -20; lon_hi = -10
Region = '2' ; Region_f = 'Two'; lat_lo = 5; lat_hi = 15; lon_lo = -30; lon_hi = -20
Region = '3' ; Region_f = 'Three'; lat_lo = 5; lat_hi = 15; lon_lo = -40; lon_hi = -30

Region = '4' ; Region_f = 'Four'; lat_lo = 5; lat_hi = 15; lon_lo = -50; lon_hi = -40
Region = '5' ; Region_f = 'Five'; lat_lo = 5; lat_hi = 15; lon_lo = -60; lon_hi = -50
Region = '6' ; Region_f = 'Six'; lat_lo = 5; lat_hi = 15; lon_lo = -70; lon_hi = -60
Region = '7' ; Region_f = 'Seven'; lat_lo = 5; lat_hi = 15; lon_lo = -80; lon_hi = -70

'''here we populate a list of granules in target 10 x 10 degree region, given .txt file'''
gran_list = []
f_name = '/home/tsteckel/text/Test_Region_' + str(Region_f) + '.txt'         
#f_name = '/home/tsteckel/text/Test_Region_' + str(Region_f) + '_Appended.txt'  #for testing appended .txt

with open (f_name, "r") as myfile:     # specific to file naming convention
	data=myfile.read().splitlines()
	for each_granule in data:
		gran_number = each_granule[-5:]	
		gran_list.append(gran_number)

def slicing(values_list):
	list_vals = [values_list[i::rec] for i,j in enumerate(values_list)]
	ctt_info = list_vals[:rec]
	return ctt_info
#Since we processed the data for every 5th/10 granule in CloudSat, we are now interpolating the data to get all points to match the small rings. This should be removable when just considering rings 1,2,3 and CSat info
#============================================================================================================================================

def Check_Equal(lst):
	return lst[1:] == lst[:-1]

file_name = []
all_convective_ctt = []; all_strat_ctt = []; all_shallow_ctt = []; all_non_ctt = []
error_grans = []

'''here we populate lists of .dat files , two different data sources so discern by that'''

for g in gran_list:
	granule = '%05d' % (int(g),)
	if Region in ['1','2','3']:                           # As to not seperate GOES/Meteosat
		filename_list = glob.glob('/../METEOSAT_CLOUDSAT_SMALL' + '*' + granule + '.dat')
	if Region in ['4','5','6','7']:                           # As to not seperate GOES/Meteosat
		filename_list = glob.glob('/../PATMOSX_CLOUDSAT_SMALL' + '*' + granule + '.dat')  
	if len(filename_list) > 0:
		file_name.append(filename_list)
        
'''status check'''        
stat_count = 0
denom = float(len(file_name))
for all_files in sorted(file_name):
	try:
		stat_count += 1
		print((stat_count)/(denom), ' complete')    # percentage complete
        ''' here we initialize read in parameter lists'''

		mask = []

		ring_one = [] # ring one average ctt
		ring_two = []
		ring_three = []
		ring_one_std = []
		ring_two_std = []
		ring_three_std = []

		ring_one_cod = [] # ring one average cod
		ring_two_cod = []
		ring_three_cod = []
		ring_one_cod_std = []
		ring_two_cod_std = []
		ring_three_cod_std = []

		csat_conv = [] #CloudSat CPR flag
		center_list = []
		center_cod = []
		center_height = []
		csatyr_jd = []
		cape = []
		rr = [] 
		o500 = []
        
        '''here we initialize bins after sorting process'''

		ring_one_conv = []; ring_two_conv = []; ring_three_conv = [] ; ring_one_std_conv = []; ring_two_std_conv = []; ring_three_std_conv = [] ; ring_one_cod_conv = []; ring_two_cod_conv = []; ring_three_cod_conv = [] ; ring_one_cod_std_conv = []; ring_two_cod_std_conv = []; ring_three_cod_std_conv = [] 
		center_ctt_conv = [] ; center_cth_conv = []; center_cod_conv = []  ; csat_cth_conv = []; csatyr_jd_conv = []; cape_conv = []; rr_conv = []; o500_conv = []

		ring_one_strat = []; ring_two_strat = []; ring_three_strat = [] ; ring_one_std_strat = []; ring_two_std_strat = []; ring_three_std_strat = [] ; ring_one_cod_strat = []; ring_two_cod_strat = []; ring_three_cod_strat = [] ; ring_one_cod_std_strat = []; ring_two_cod_std_strat = []; ring_three_cod_std_strat = [] 
		center_ctt_strat = [] ; center_cth_strat = []; center_cod_strat = []  ; csat_cth_strat = []; csatyr_jd_strat = []; cape_strat = []; rr_strat = []; o500_strat = []

		ring_one_shallow = []; ring_two_shallow = []; ring_three_shallow = [] ; ring_one_std_shallow = []; ring_two_std_shallow = []; ring_three_std_shallow = [] ; ring_one_cod_shallow = []; ring_two_cod_shallow = []; ring_three_cod_shallow = [] ; ring_one_cod_std_shallow = []; ring_two_cod_std_shallow = []; ring_three_cod_std_shallow = [] 
		center_ctt_shallow = [] ; center_cth_shallow = []; center_cod_shallow = []  ; csat_cth_shallow = []; csatyr_jd_shallow = []; cape_shallow = []; rr_shallow = []; o500_shallow = []

		ring_one_non = []; ring_two_non = []; ring_three_non = [] ; ring_one_std_non = []; ring_two_std_non = []; ring_three_std_non = [] ; ring_one_cod_non = []; ring_two_cod_non = []; ring_three_cod_non = [] ; ring_one_cod_std_non = []; ring_two_cod_std_non = []; ring_three_cod_std_non = [] 
		center_ctt_non = [] ; center_cth_non = []; center_cod_non = []  ; csat_cth_non = []; csatyr_jd_non = []; cape_non = []; rr_non = []; o500_non = []

		for ind_file in all_files:
			print(ind_file)
			lat_list = []; lon_list = []; time_list = []; conv_list = []; ref_list = []
			cape_list = []; o500_list = []; rr_list = []

			cttcent_list = []; cthcent_list = []; codcent_list = []
			one_cttmean_list = []; two_cttmean_list = []; three_cttmean_list = []
			one_cttstd_list = []; two_cttstd_list = []; three_cttstd_list = []

			one_codmean_list = []; two_codmean_list = []; three_codmean_list = []
			one_codstd_list = []; two_codstd_list = []; three_codstd_list = []

			patmos_cloudsat_file = open(ind_file, 'rb')

			load_file = np.load(patmos_cloudsat_file)
            
            '''bin is cPickled file, so we read in like a dictionary'''
			for key in load_file:
				if len(key) > 100:
					overpass_time  = key['CloudSat Overpass Time']; granule_number = key['Granule Number']; time = key['Time']
					time_list.append(time)
					csatyr_jd = int(ind_file[-23:-16])

					if 'SMALL' in ind_file:	

						ctt_center = key['CTT Center Pixel']; cth_center = key['CTH Center Pixel']; cod_center = key['COD Center Pixel']; lat = key['Latitude']; lon = key['Longitude']; reflectivity = key['Reflectivity (dBZ)']; conv_flag = key['Convective Flag']

						ctt_one_mean = key['Mean CTT 1']; ctt_two_mean = key['Mean CTT 2']; ctt_three_mean = key['Mean CTT 3']
						ctt_one_std = key['STD CTT 1']; ctt_two_std = key['STD CTT 2']; ctt_three_std = key['STD CTT 3']
						cod_one_mean = key['Mean COD 1']; cod_two_mean = key['Mean COD 2']; cod_three_mean = key['Mean COD 3']
						cod_one_std = key['STD COD 1']; cod_two_std = key['STD COD 2']; cod_three_std = key['STD COD 3']

						rr_read = key['Rain Rate (mm/hr)']; o500_read = key['OMEGA500']; cape_read = key['CAPE']

						cttcent_list.append(ctt_center); cthcent_list.append(cth_center);lat_list.append(lat); lon_list.append(lon); ref_list.append(reflectivity); conv_list.append(conv_flag); rr_list.append(rr_read); cape_list.append(cape_read); o500_list.append(o500_read)  ; codcent_list.append(cod_center)
						one_cttmean_list.append(ctt_one_mean); two_cttmean_list.append(ctt_two_mean); three_cttmean_list.append(ctt_three_mean) 
						one_cttstd_list.append(ctt_one_std); two_cttstd_list.append(ctt_two_std); three_cttstd_list.append(ctt_three_std) 
						one_codmean_list.append(cod_one_mean); two_codmean_list.append(cod_two_mean); three_codmean_list.append(cod_three_mean) 
						one_codstd_list.append(cod_one_std); two_codstd_list.append(cod_two_std); three_codstd_list.append(cod_three_std) 

			mask_lat = np.logical_and(np.less(lat_lo, lat_list), np.less(lat_list, lat_hi))
			mask_lon = np.logical_and(np.less(lon_lo, lon_list), np.less(lon_list, lon_hi))
			mask_select = np.logical_and(mask_lon, mask_lat)

            '''here we account for temporal lag, to chose best file time''' 
			if len(time_list) > 0:
				rec = len(set([x for x in time_list if time_list.count(x) > 1])) #This is giving the number of times that are included in the plus 3, overpass and minus 6
				times = []
				for i in time_list[0:rec]: #These are the times that are in the patmos/cldsat file.
					times.append(int(i))
				time_closest_to_overpass =  min(times, key=lambda f:abs(f-int(overpass_time[0:4])))   # Voodoo magic
				mask_sub = np.asarray((slicing(mask_select)))      # This one seems out of place
				ind_time = int('%04d' % (time_closest_to_overpass))
				if ((int(overpass_time[0:4])-ind_time) > 35 and (int(overpass_time[0:4])-ind_time) < 45):
					ind_time = int(('%04d' % (time_closest_to_overpass))[0:2]) + 1
				else:
					ind_time = int(('%04d' % (time_closest_to_overpass))[0:2])
				op_time = '%06d' % (int(overpass_time),)

				center = slicing(cttcent_list)
				height = slicing(cthcent_list)
				cod = slicing(codcent_list)

				one = slicing(one_cttmean_list)
				two = slicing(two_cttmean_list)
				three = slicing(three_cttmean_list)

				one_std = slicing(one_cttstd_list)
				two_std = slicing(two_cttstd_list)
				three_std = slicing(three_cttstd_list)

				one_cod = slicing(one_codmean_list)
				two_cod = slicing(two_codmean_list)
				three_cod = slicing(three_codmean_list)

				one_cod_std = slicing(one_codstd_list)
				two_cod_std = slicing(two_codstd_list)
				three_cod_std = slicing(three_codstd_list)

				cape_sub = slicing(cape_list)
				o500_sub = slicing(o500_list)

				csat_sub = np.asarray(conv_list[::rec])
				rr_sub = np.absolute(np.asarray(rr_list[::rec]))   # np.absolute because rain rate will report negatives 

				if len(one) > 0:      # append to lists at beginning of for loop
					center_list.append(center)
					center_height.append(height)
					center_cod.append(cod)

					ring_one.append(one)
					ring_two.append(two)
					ring_three.append(three)

					ring_one_std.append(one_std)
					ring_two_std.append(two_std)
					ring_three_std.append(three_std)

					ring_one_cod.append(one_cod)
					ring_two_cod.append(two_cod)
					ring_three_cod.append(three_cod)

					ring_one_cod_std.append(one_cod_std)
					ring_two_cod_std.append(two_cod_std)
					ring_three_cod_std.append(three_cod_std)

					o500.append(o500_sub)
					cape.append(cape_sub)

					mask.append(mask_sub)
					csat_conv.append(csat_sub)
					rr.append(rr_sub)
		#==========================================
		#==========================================
		if len(csat_conv) > 0 and len(ref_list) > 0 and len(ring_one[0]) >= ind_time:
			if len(ring_one) > 0:
				for iofct in range(ind_time,ind_time +1):
					index_list = []
					ref_list = np.asarray(ref_list[::rec])
					large_index =  int(len(ring_one[0][iofct]))
					mask_final = mask[0][iofct]
					if len(mask_final) < large_index:                # When Mask is smaller than patmos-x pixels 
						large_index = len(mask_final)
					mask_final = list(mask_final)[0:large_index]
					if not any(mask_final):
						print('Should revert to next granule file')
						continue
					mask_final = np.asarray(mask_final)

					center_final = list(np.asarray(center_list[0][iofct][0:large_index])[mask_final])
					height_final = list(np.asarray(center_height[0][iofct][0:large_index])[mask_final])
					cod_final = list(np.asarray(center_cod[0][iofct][0:large_index])[mask_final])
					cape_final = list(np.asarray(cape[0][iofct][0:large_index])[mask_final])
					o500_final = list(np.asarray(o500[0][iofct][0:large_index])[mask_final])

					ring_one_final = list(np.asarray(ring_one[0][iofct][0:large_index])[mask_final])
					ring_two_final = list(np.asarray(ring_two[0][iofct][0:large_index])[mask_final])
					ring_three_final = list(np.asarray(ring_three[0][iofct][0:large_index])[mask_final])

					ring_one_std_final = list(np.asarray(ring_one_std[0][iofct][0:large_index])[mask_final])
					ring_two_std_final = list(np.asarray(ring_two_std[0][iofct][0:large_index])[mask_final])
					ring_three_std_final = list(np.asarray(ring_three_std[0][iofct][0:large_index])[mask_final])

					ring_one_cod_final = list(np.asarray(ring_one_cod[0][iofct][0:large_index])[mask_final])
					ring_two_cod_final = list(np.asarray(ring_two_cod[0][iofct][0:large_index])[mask_final])
					ring_three_cod_final = list(np.asarray(ring_three_cod[0][iofct][0:large_index])[mask_final])

					ring_one_cod_std_final = list(np.asarray(ring_one_cod_std[0][iofct][0:large_index])[mask_final])
					ring_two_cod_std_final = list(np.asarray(ring_two_cod_std[0][iofct][0:large_index])[mask_final])
					ring_three_cod_std_final = list(np.asarray(ring_three_cod_std[0][iofct][0:large_index])[mask_final])

					ref_list_masked = list(np.array(ref_list[0:large_index])[mask_final])
					lat_list_masked = list(np.asarray(lat_list[:large_index])[mask_final])

					csat_conv = csat_conv[0][:large_index]
					csat_conv = [x[0] for x in csat_conv]
					csat_conv = list(np.array(csat_conv)[mask_final])

					rr_final = rr[0][:large_index]
					rr_final = [x[0] for x in rr_final]
					rr_final = list(np.array(rr_final)[mask_final])

		#=========================Csat_cth================

					ref_hght_start = -5
					ref_hght_end = 30.1
					ref_delta_hght = .28
					ref_hght_bins = np.arange(ref_hght_start, ref_hght_end, ref_delta_hght)
					ref_arr = np.asarray(ref_list_masked)
					lat_arr  = np.asarray(lat_list_masked)
					ref_arr_final   = masked_where(ref_arr <= -20, ref_arr); ref_arr_final = masked_where(ref_arr_final > 50., ref_arr_final); ref_final = np.fliplr(ref_arr_final)
                    
					'''here we populate CloudSat specific parameters, for filtering downstream ''' 
					validation_triple = []

					for i in range(len(height_final)):
						triple = []
						ref_ind = np.where((ref_final[i] > -20) & (ref_final[i] <= 50 )) 
						csat_cth = ref_hght_bins[ref_ind[0][-1]]
						triple.append(0)
						triple.append(csat_cth)
						triple.append(height_final[i])
						triple= np.asarray(triple, dtype = 'float')
						validation_triple.append(triple)
                        
                    '''here we initiate grouping of flags and decide record/no record'''

				   #======New Grouping===========
					if len(ring_one_final) > 0:
						current_group = []
						conv = []
						zero = []
						zero_ind = []
						deep_ind = []
						strat_ind = []
						shallow_ind = []
						pixel_size_conv = []
						pixel_size_strat = []
						pixel_size_shallow = []
						for i in range(len(csat_conv)):
							if csat_conv[i] == -1 or csat_conv[i] == -2 or csat_conv[i] == -3:
								csat_conv[i] = 0
							if csat_conv[i] == 0:                                          # Zero List
								current_group.append(i)
							else:
								if len(current_group) > 0:
									zero.append(current_group)
									current_group = []
							if i == (len(csat_conv) - 1):                                  # Cutting off last pixel of flag list
								if len(current_group) > 0:
									zero.append(current_group)
								 
						#=====Conv Pixels====

						if zero[0][0] != 0:                                #Start
							current_group = []
							conv_len = zero[0][0]
							for j in range(0,conv_len):
								conv_ind = j
								current_group.append(conv_ind)
							conv.append(current_group)

						for i in range(len(zero)-1):                     #Middle
							current_group = []
							conv_len = zero[i+1][0]-zero[i][-1]
							for j in range(1,conv_len):
								conv_ind = zero[i][-1] + j
								current_group.append(conv_ind)
							conv.append(current_group)

						if zero[-1][-1] != (len(csat_conv) - 1):                   #End
							current_group = []
							conv_len = len(csat_conv) - zero[-1][-1]
							for j in range(1,conv_len):
								conv_ind = zero[-1][-1] + j
								current_group.append(conv_ind)
							conv.append(current_group)

						#======Getting non indices 
						for i in zero:
							if len(i) > 100:
								ind_off = len(i)/2
								ind = i[0] + ind_off
								zero_ind.append(ind)
						#======Getting conv indices
						for i in conv:
							test_list = [csat_conv[j] for j in i]
							if Check_Equal:
								if test_list[0] == 1:
									len_list = np.repeat(len(i),len(i))
									zipped_ind = zip(i,len_list)
									deep_ind += zipped_ind
								if test_list[0] == 2:
									len_list = np.repeat(len(i),len(i))
									zipped_ind = zip(i,len_list)
									strat_ind += zipped_ind
								if test_list[0] == 3:
									len_list = np.repeat(len(i),len(i))
									zipped_ind = zip(i,len_list)
									shallow_ind += zipped_ind
						for i in range(len(deep_ind)):  

							ring_one_conv.append(ring_one_final[(deep_ind[i][0])]) #Using the flag indices to get the ctt values we want.
							ring_two_conv.append(ring_two_final[(deep_ind[i][0])])
							ring_three_conv.append(ring_three_final[(deep_ind[i][0])])
							ring_one_std_conv.append(ring_one_std_final[(deep_ind[i][0])]) #Using the flag indices to get the ctt values we want
							ring_two_std_conv.append(ring_two_std_final[(deep_ind[i][0])])
							ring_three_std_conv.append(ring_three_std_final[(deep_ind[i][0])])
							ring_one_cod_conv.append(ring_one_cod_final[(deep_ind[i][0])]) #Using the flag indices to get the ctt values we want.
							ring_two_cod_conv.append(ring_two_cod_final[(deep_ind[i][0])])
							ring_three_cod_conv.append(ring_three_cod_final[(deep_ind[i][0])])
							ring_one_cod_std_conv.append(ring_one_cod_std_final[(deep_ind[i][0])]) #Using the flag indices to get the ctt values we want.
							ring_two_cod_std_conv.append(ring_two_cod_std_final[(deep_ind[i][0])])
							ring_three_cod_std_conv.append(ring_three_cod_std_final[(deep_ind[i][0])])

							cape_conv.append(cape_final[(deep_ind[i][0])])
							o500_conv.append(o500_final[(deep_ind[i][0])])
							center_cod_conv.append(cod_final[(deep_ind[i][0])])
							center_ctt_conv.append(center_final[(deep_ind[i][0])])
							center_cth_conv.append(height_final[(deep_ind[i][0])])
							csat_cth_conv.append(validation_triple[(deep_ind[i][0])][1])
							pixel_size_conv.append(deep_ind[i][1])
							csatyr_jd_conv.append(csatyr_jd)
							rr_conv.append(rr_final[(deep_ind[i][0])])                           # might have to change, array

						for i in range(len(strat_ind)): 
							ring_one_strat.append(ring_one_final[(strat_ind[i][0])]) #Using the flag indices to get the ctt values we want.
							ring_two_strat.append(ring_two_final[(strat_ind[i][0])])
							ring_three_strat.append(ring_three_final[(strat_ind[i][0])])
							ring_one_std_strat.append(ring_one_std_final[(strat_ind[i][0])]) #Using the flag indices to get the ctt values we want
							ring_two_std_strat.append(ring_two_std_final[(strat_ind[i][0])])
							ring_three_std_strat.append(ring_three_std_final[(strat_ind[i][0])])
							ring_one_cod_strat.append(ring_one_cod_final[(strat_ind[i][0])]) #Using the flag indices to get the ctt values we want.
							ring_two_cod_strat.append(ring_two_cod_final[(strat_ind[i][0])])
							ring_three_cod_strat.append(ring_three_cod_final[(strat_ind[i][0])])
							ring_one_cod_std_strat.append(ring_one_cod_std_final[(strat_ind[i][0])]) #Using the flag indices to get the ctt values we want.
							ring_two_cod_std_strat.append(ring_two_cod_std_final[(strat_ind[i][0])])
							ring_three_cod_std_strat.append(ring_three_cod_std_final[(strat_ind[i][0])])

							cape_strat.append(cape_final[(strat_ind[i][0])])
							o500_strat.append(o500_final[(strat_ind[i][0])])
							center_cod_strat.append(cod_final[(strat_ind[i][0])])
							center_ctt_strat.append(center_final[(strat_ind[i][0])])
							center_cth_strat.append(height_final[(strat_ind[i][0])])
							csat_cth_strat.append(validation_triple[(strat_ind[i][0])][1])
							pixel_size_strat.append(strat_ind[i][1])
							csatyr_jd_strat.append(csatyr_jd)
							rr_strat.append(rr_final[(strat_ind[i][0])])     

						for i in range(len(shallow_ind)): 
							ring_one_shallow.append(ring_one_final[(shallow_ind[i][0])]) #Using the flag indices to get the ctt values we want.
							ring_two_shallow.append(ring_two_final[(shallow_ind[i][0])])
							ring_three_shallow.append(ring_three_final[(shallow_ind[i][0])])
							ring_one_std_shallow.append(ring_one_std_final[(shallow_ind[i][0])]) #Using the flag indices to get the ctt values we want
							ring_two_std_shallow.append(ring_two_std_final[(shallow_ind[i][0])])
							ring_three_std_shallow.append(ring_three_std_final[(shallow_ind[i][0])])
							ring_one_cod_shallow.append(ring_one_cod_final[(shallow_ind[i][0])]) #Using the flag indices to get the ctt values we want.
							ring_two_cod_shallow.append(ring_two_cod_final[(shallow_ind[i][0])])
							ring_three_cod_shallow.append(ring_three_cod_final[(shallow_ind[i][0])])
							ring_one_cod_std_shallow.append(ring_one_cod_std_final[(shallow_ind[i][0])]) #Using the flag indices to get the ctt values we want.
							ring_two_cod_std_shallow.append(ring_two_cod_std_final[(shallow_ind[i][0])])
							ring_three_cod_std_shallow.append(ring_three_cod_std_final[(shallow_ind[i][0])])

							cape_shallow.append(cape_final[(shallow_ind[i][0])])
							o500_shallow.append(o500_final[(shallow_ind[i][0])])
							center_cod_shallow.append(cod_final[(shallow_ind[i][0])])
							center_ctt_shallow.append(center_final[(shallow_ind[i][0])])
							center_cth_shallow.append(height_final[(shallow_ind[i][0])])
							csat_cth_shallow.append(validation_triple[(shallow_ind[i][0])][1])
							pixel_size_shallow.append(shallow_ind[i][1])
							csatyr_jd_shallow.append(csatyr_jd)
							rr_shallow.append(rr_final[(shallow_ind[i][0])])     

						for i in zero_ind: 
							ring_one_non.append(ring_one_final[i]) #Using the flag indices to get the ctt values we want.
							ring_two_non.append(ring_two_final[i])
							ring_three_non.append(ring_three_final[i])
							ring_one_std_non.append(ring_one_std_final[i]) #Using the flag indices to get the ctt values we want
							ring_two_std_non.append(ring_two_std_final[i])
							ring_three_std_non.append(ring_three_std_final[i])
							ring_one_cod_non.append(ring_one_cod_final[i]) #Using the flag indices to get the ctt values we want.
							ring_two_cod_non.append(ring_two_cod_final[i])
							ring_three_cod_non.append(ring_three_cod_final[i])
							ring_one_cod_std_non.append(ring_one_cod_std_final[i]) #Using the flag indices to get the ctt values we want.
							ring_two_cod_std_non.append(ring_two_cod_std_final[i])
							ring_three_cod_std_non.append(ring_three_cod_std_final[i])

							cape_non.append(cape_final[i])
							o500_non.append(o500_final[i])
							center_cod_non.append(cod_final[i])
							center_ctt_non.append(center_final[i])
							center_cth_non.append(height_final[i])
							csat_cth_non.append(validation_triple[i][1])
							csatyr_jd_non.append(csatyr_jd)
							rr_non.append(rr_final[i])     

					all_convective_ctt.append((ring_one_conv, ring_two_conv, ring_three_conv, center_ctt_conv, center_cth_conv, csat_cth_conv, pixel_size_conv, csatyr_jd_conv, ring_one_std_conv, ring_two_std_conv, ring_three_std_conv,ring_one_cod_conv, ring_two_cod_conv, ring_three_cod_conv,ring_one_cod_std_conv, ring_two_cod_std_conv, ring_three_cod_std_conv, center_cod_conv, o500_conv, cape_conv, rr_conv))
					all_strat_ctt.append((ring_one_strat, ring_two_strat, ring_three_strat, center_ctt_strat, center_cth_strat, csat_cth_strat, pixel_size_strat, csatyr_jd_strat, ring_one_std_strat, ring_two_std_strat, ring_three_std_strat,ring_one_cod_strat, ring_two_cod_strat, ring_three_cod_strat,ring_one_cod_std_strat, ring_two_cod_std_strat, ring_three_cod_std_strat, center_cod_strat, o500_strat, cape_strat, rr_strat))
					all_shallow_ctt.append((ring_one_shallow, ring_two_shallow, ring_three_shallow, center_ctt_shallow, center_cth_shallow, csat_cth_shallow, pixel_size_shallow, csatyr_jd_shallow, ring_one_std_shallow, ring_two_std_shallow, ring_three_std_shallow,ring_one_cod_shallow, ring_two_cod_shallow, ring_three_cod_shallow,ring_one_cod_std_shallow, ring_two_cod_std_shallow, ring_three_cod_std_shallow, center_cod_shallow, o500_shallow, cape_shallow, rr_shallow))
					all_non_ctt.append((ring_one_non, ring_two_non, ring_three_non, center_ctt_non, center_cth_non, csat_cth_non, csatyr_jd_non, ring_one_std_non, ring_two_std_non, ring_three_std_non,ring_one_cod_non, ring_two_cod_non, ring_three_cod_non,ring_one_cod_std_non, ring_two_cod_std_non, ring_three_cod_std_non, center_cod_non, o500_non, cape_non, rr_non))

	except KeyboardInterrupt:
		raise
	except: 
		error_grans.append(all_files)	
		print ' errant file...' 					
		continue 
    
'''here we export resulting reference vectors as .txt. There are four text files, one for each CPR flag'''
#Non Convective#

ring_one_out_non = []             # Creating empty lists that will contain CTT for every ring for all granules 
ring_two_out_non = []
ring_three_out_non = []
ring_one_std_out_non = []            
ring_two_std_out_non = []
ring_three_std_out_non = []
ring_one_cod_out_non = []        
ring_two_cod_out_non = []
ring_three_cod_out_non = []
ring_one_cod_std_out_non = []        
ring_two_cod_std_out_non = []
ring_three_cod_std_out_non = []
rr_out_non = []
cape_out_non = []
o500_out_non = []
center_ctt_out_non = []
center_cth_out_non = []
center_cod_out_non = []
csat_cth_out_non = []
csatyr_jd_out_non = []

for x in range(len(all_non_ctt)):           # appending the ring data of seperate granules together, maintain first 7 for for backwards compatability
	
	ring_one_out_non += all_non_ctt[x][0]
	ring_two_out_non += all_non_ctt[x][1]
	ring_three_out_non += all_non_ctt[x][2]
	center_ctt_out_non  += all_non_ctt[x][3]
	center_cth_out_non += all_non_ctt[x][4]
	csat_cth_out_non += all_non_ctt[x][5]
	csatyr_jd_out_non += all_non_ctt[x][6]
	ring_one_std_out_non += all_non_ctt[x][7]          
	ring_two_std_out_non += all_non_ctt[x][8]
	ring_three_std_out_non += all_non_ctt[x][9]
	ring_one_cod_out_non += all_non_ctt[x][10]        
	ring_two_cod_out_non += all_non_ctt[x][11]
	ring_three_cod_out_non += all_non_ctt[x][12]
	ring_one_cod_std_out_non += all_non_ctt[x][13]        
	ring_two_cod_std_out_non += all_non_ctt[x][14]
	ring_three_cod_std_out_non += all_non_ctt[x][15]
	center_cod_out_non += all_non_ctt[x][16]
	o500_out_non += all_non_ctt[x][17]
	cape_out_non += all_non_ctt[x][18]
	rr_out_non += all_non_ctt[x][19]

filename = '/../Region_Increased_Dimensionality_non' + str(Region) + '.txt'
with open(filename, 'wb') as f: 

	f.write('ring_one|' +str(ring_one_out_non) +'\n')    
	f.write('ring_two|' +str(ring_two_out_non) +'\n')
	f.write('ring_three|' +str(ring_three_out_non) +'\n')
	f.write('center_ctt|' +str(center_ctt_out_non) +'\n')
	f.write('center_cth|' +str(center_cth_out_non) +'\n')
	f.write('csat_cth|' +str(csat_cth_out_non) +'\n')
	f.write('csat_yrjd|' +str(csatyr_jd_out_non) +'\n')
	f.write('ring_one_std|' +str(ring_one_std_out_non) +'\n')      
	f.write('ring_two_std|' +str(ring_two_std_out_non) +'\n')
	f.write('ring_three_std|' +str(ring_three_std_out_non) +'\n')
	f.write('ring_one_cod|' +str(ring_one_cod_out_non) +'\n')        
	f.write('ring_two_cod|' +str(ring_two_cod_out_non) +'\n')
	f.write('ring_three_cod|' +str(ring_three_cod_out_non) +'\n')
	f.write('ring_one_cod_std|' +str(ring_one_cod_std_out_non) +'\n')      
	f.write('ring_two_cod_std|' +str(ring_two_cod_std_out_non) +'\n')
	f.write('ring_three_cod_std|' +str(ring_three_cod_std_out_non) +'\n')
	f.write('center_cod|' +str(center_cod_out_non) +'\n')
	f.write('o500|' +str(o500_out_non) +'\n')
	f.write('cape|' +str(cape_out_non) +'\n')
	f.write('rain_rate(mm/hr)|' +str(rr_out_non) +'\n')

f.close()
#=================================================================================#
#Deep#

ring_one_out_deep = []             # Creating empty lists that will contain CTT for every ring for all granules 
ring_two_out_deep = []
ring_three_out_deep = []
ring_one_std_out_deep = []            
ring_two_std_out_deep = []
ring_three_std_out_deep = []
ring_one_cod_out_deep = []        
ring_two_cod_out_deep = []
ring_three_cod_out_deep = []
ring_one_cod_std_out_deep = []        
ring_two_cod_std_out_deep = []
ring_three_cod_std_out_deep = []
rr_out_deep = []
cape_out_deep = []
o500_out_deep = []
center_ctt_out_deep = []
center_cth_out_deep = []
center_cod_out_deep = []
csat_cth_out_deep = []
pixel_size_out_deep = []
csatyr_jd_out_deep = []

for x in range(len(all_convective_ctt)):           # appending the ring data of seperate granules together
	
	ring_one_out_deep += all_convective_ctt[x][0]
	ring_two_out_deep += all_convective_ctt[x][1]
	ring_three_out_deep += all_convective_ctt[x][2]
	center_ctt_out_deep  += all_convective_ctt[x][3]
	center_cth_out_deep += all_convective_ctt[x][4]
	csat_cth_out_deep += all_convective_ctt[x][5]
	pixel_size_out_deep += all_convective_ctt[x][6]
	csatyr_jd_out_deep += all_convective_ctt[x][7]
	ring_one_std_out_deep += all_convective_ctt[x][8]          
	ring_two_std_out_deep += all_convective_ctt[x][9]
	ring_three_std_out_deep += all_convective_ctt[x][10]
	ring_one_cod_out_deep += all_convective_ctt[x][11]        
	ring_two_cod_out_deep += all_convective_ctt[x][12]
	ring_three_cod_out_deep += all_convective_ctt[x][13]
	ring_one_cod_std_out_deep += all_convective_ctt[x][14]        
	ring_two_cod_std_out_deep += all_convective_ctt[x][15]
	ring_three_cod_std_out_deep += all_convective_ctt[x][16]
	center_cod_out_deep += all_convective_ctt[x][17]
	o500_out_deep += all_convective_ctt[x][18]
	cape_out_deep += all_convective_ctt[x][19]
	rr_out_deep += all_convective_ctt[x][20]

filename = '/../Region_Increased_Dimensionality_deep' + str(Region) + '.txt'
with open(filename, 'wb') as f: 

	f.write('ring_one|' +str(ring_one_out_deep) +'\n')    
	f.write('ring_two|' +str(ring_two_out_deep) +'\n')
	f.write('ring_three|' +str(ring_three_out_deep) +'\n')
	f.write('center_ctt|' +str(center_ctt_out_deep) +'\n')
	f.write('center_cth|' +str(center_cth_out_deep) +'\n')
	f.write('csat_cth|' +str(csat_cth_out_deep) +'\n')
	f.write('pixel size|' +str(pixel_size_out_deep) +'\n')
	f.write('csat_yrjd|' +str(csatyr_jd_out_deep) +'\n')
	f.write('ring_one_std|' +str(ring_one_std_out_deep) +'\n')      
	f.write('ring_two_std|' +str(ring_two_std_out_deep) +'\n')
	f.write('ring_three_std|' +str(ring_three_std_out_deep) +'\n')
	f.write('ring_one_cod|' +str(ring_one_cod_out_deep) +'\n')        
	f.write('ring_two_cod|' +str(ring_two_cod_out_deep) +'\n')
	f.write('ring_three_cod|' +str(ring_three_cod_out_deep) +'\n')
	f.write('ring_one_cod_std|' +str(ring_one_cod_std_out_deep) +'\n')      
	f.write('ring_two_cod_std|' +str(ring_two_cod_std_out_deep) +'\n')
	f.write('ring_three_cod_std|' +str(ring_three_cod_std_out_deep) +'\n')
	f.write('center_cod|' +str(center_cod_out_deep) +'\n')
	f.write('o500|' +str(o500_out_deep) +'\n')
	f.write('cape|' +str(cape_out_deep) +'\n')
	f.write('rain_rate(mm/hr)|' +str(rr_out_deep) +'\n')
f.close()
#==================================================================================#
#Shallow#

ring_one_out_shallow = []             # Creating empty lists that will contain CTT for every ring for all granules 
ring_two_out_shallow = []
ring_three_out_shallow = []
ring_one_std_out_shallow = []            
ring_two_std_out_shallow = []
ring_three_std_out_shallow = []
ring_one_cod_out_shallow = []        
ring_two_cod_out_shallow = []
ring_three_cod_out_shallow = []
ring_one_cod_std_out_shallow = []        
ring_two_cod_std_out_shallow = []
ring_three_cod_std_out_shallow = []
rr_out_shallow = []
cape_out_shallow = []
o500_out_shallow = []
center_ctt_out_shallow = []
center_cth_out_shallow = []
center_cod_out_shallow = []
csat_cth_out_shallow = []
pixel_size_out_shallow = []
csatyr_jd_out_shallow = []

for x in range(len(all_shallow_ctt)):           # appending the ring data of seperate granules together
	
	ring_one_out_shallow += all_shallow_ctt[x][0]
	ring_two_out_shallow += all_shallow_ctt[x][1]
	ring_three_out_shallow += all_shallow_ctt[x][2]
	center_ctt_out_shallow  += all_shallow_ctt[x][3]
	center_cth_out_shallow += all_shallow_ctt[x][4]
	csat_cth_out_shallow += all_shallow_ctt[x][5]
	pixel_size_out_shallow += all_shallow_ctt[x][6]
	csatyr_jd_out_shallow += all_shallow_ctt[x][7]
	ring_one_std_out_shallow += all_shallow_ctt[x][8]          
	ring_two_std_out_shallow += all_shallow_ctt[x][9]
	ring_three_std_out_shallow += all_shallow_ctt[x][10]
	ring_one_cod_out_shallow += all_shallow_ctt[x][11]        
	ring_two_cod_out_shallow += all_shallow_ctt[x][12]
	ring_three_cod_out_shallow += all_shallow_ctt[x][13]
	ring_one_cod_std_out_shallow += all_shallow_ctt[x][14]        
	ring_two_cod_std_out_shallow += all_shallow_ctt[x][15]
	ring_three_cod_std_out_shallow += all_shallow_ctt[x][16]
	center_cod_out_shallow += all_shallow_ctt[x][17]
	o500_out_shallow += all_shallow_ctt[x][18]
	cape_out_shallow += all_shallow_ctt[x][19]
	rr_out_shallow += all_shallow_ctt[x][20]

filename = '/../Region_Increased_Dimensionality_shallow' + str(Region) + '.txt'
with open(filename, 'wb') as f: 

	f.write('ring_one|' +str(ring_one_out_shallow) +'\n') 
	f.write('ring_two|' +str(ring_two_out_shallow) +'\n')
	f.write('ring_three|' +str(ring_three_out_shallow) +'\n')
	f.write('center_ctt|' +str(center_ctt_out_shallow) +'\n')
	f.write('center_cth|' +str(center_cth_out_shallow) +'\n')
	f.write('csat_cth|' +str(csat_cth_out_shallow) +'\n')
	f.write('pixel size|' +str(pixel_size_out_shallow) +'\n')
	f.write('csat_yrjd|' +str(csatyr_jd_out_shallow) +'\n')
	f.write('ring_one_std|' +str(ring_one_std_out_shallow) +'\n')     
	f.write('ring_two_std|' +str(ring_two_std_out_shallow) +'\n')
	f.write('ring_three_std|' +str(ring_three_std_out_shallow) +'\n')
	f.write('ring_one_cod|' +str(ring_one_cod_out_shallow) +'\n')        
	f.write('ring_two_cod|' +str(ring_two_cod_out_shallow) +'\n')
	f.write('ring_three_cod|' +str(ring_three_cod_out_shallow) +'\n')
	f.write('ring_one_cod_std|' +str(ring_one_cod_std_out_shallow) +'\n') 
	f.write('ring_two_cod_std|' +str(ring_two_cod_std_out_shallow) +'\n')
	f.write('ring_three_cod_std|' +str(ring_three_cod_std_out_shallow) +'\n')
	f.write('center_cod|' +str(center_cod_out_shallow) +'\n')
	f.write('o500|' +str(o500_out_shallow) +'\n')
	f.write('cape|' +str(cape_out_shallow) +'\n')
	f.write('rain_rate(mm/hr)|' +str(rr_out_shallow) +'\n')
f.close()
#==================================================================================#
#Strat#

ring_one_out_strat = []             # Creating empty lists that will contain CTT for every ring for all granules 
ring_two_out_strat = []
ring_three_out_strat = []
ring_one_std_out_strat = []            
ring_two_std_out_strat = []
ring_three_std_out_strat = []
ring_one_cod_out_strat = []        
ring_two_cod_out_strat = []
ring_three_cod_out_strat = []
ring_one_cod_std_out_strat = []        
ring_two_cod_std_out_strat = []
ring_three_cod_std_out_strat = []
rr_out_strat = []
cape_out_strat = []
o500_out_strat = []
center_ctt_out_strat = []
center_cth_out_strat = []
center_cod_out_strat = []
csat_cth_out_strat = []
pixel_size_out_strat = []
csatyr_jd_out_strat = []

for x in range(len(all_strat_ctt)):           # appending the ring data of seperate granules together
	
	ring_one_out_strat += all_strat_ctt[x][0]
	ring_two_out_strat += all_strat_ctt[x][1]
	ring_three_out_strat += all_strat_ctt[x][2]
	center_ctt_out_strat  += all_strat_ctt[x][3]
	center_cth_out_strat += all_strat_ctt[x][4]
	csat_cth_out_strat += all_strat_ctt[x][5]
	pixel_size_out_strat += all_strat_ctt[x][6]
	csatyr_jd_out_strat += all_strat_ctt[x][7]
	ring_one_std_out_strat += all_strat_ctt[x][8]          
	ring_two_std_out_strat += all_strat_ctt[x][9]
	ring_three_std_out_strat += all_strat_ctt[x][10]
	ring_one_cod_out_strat += all_strat_ctt[x][11]        
	ring_two_cod_out_strat += all_strat_ctt[x][12]
	ring_three_cod_out_strat += all_strat_ctt[x][13]
	ring_one_cod_std_out_strat += all_strat_ctt[x][14]        
	ring_two_cod_std_out_strat += all_strat_ctt[x][15]
	ring_three_cod_std_out_strat += all_strat_ctt[x][16]
	center_cod_out_strat += all_strat_ctt[x][17]
	o500_out_strat += all_strat_ctt[x][18]
	cape_out_strat += all_strat_ctt[x][19]
	rr_out_strat += all_strat_ctt[x][20]

filename = '/../Region_Increased_Dimensionality_strat' + str(Region) + '.txt'
with open(filename, 'wb') as f: 

	f.write('ring_one|' +str(ring_one_out_strat) +'\n')   
	f.write('ring_two|' +str(ring_two_out_strat) +'\n')
	f.write('ring_three|' +str(ring_three_out_strat) +'\n')
	f.write('center_ctt|' +str(center_ctt_out_strat) +'\n')
	f.write('center_cth|' +str(center_cth_out_strat) +'\n')
	f.write('csat_cth|' +str(csat_cth_out_strat) +'\n')
	f.write('pixel size|' +str(pixel_size_out_strat) +'\n')
	f.write('csat_yrjd|' +str(csatyr_jd_out_strat) +'\n')
	f.write('ring_one_std|' +str(ring_one_std_out_strat) +'\n')        
	f.write('ring_two_std|' +str(ring_two_std_out_strat) +'\n')
	f.write('ring_three_std|' +str(ring_three_std_out_strat) +'\n')
	f.write('ring_one_cod|' +str(ring_one_cod_out_strat) +'\n')         
	f.write('ring_two_cod|' +str(ring_two_cod_out_strat) +'\n')
	f.write('ring_three_cod|' +str(ring_three_cod_out_strat) +'\n')
	f.write('ring_one_cod_std|' +str(ring_one_cod_std_out_strat) +'\n')     
	f.write('ring_two_cod_std|' +str(ring_two_cod_std_out_strat) +'\n')
	f.write('ring_three_cod_std|' +str(ring_three_cod_std_out_strat) +'\n')
	f.write('center_cod|' +str(center_cod_out_strat) +'\n')
	f.write('o500|' +str(o500_out_strat) +'\n')
	f.write('cape|' +str(cape_out_strat) +'\n')
	f.write('rain_rate(mm/hr)|' +str(rr_out_strat) +'\n')
f.close()

