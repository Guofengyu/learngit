import sys, getopt
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
plt.switch_backend('agg')

def usage():
	print 'Usage: python [options]' 
	print \
		"Options:\n" + \
		"\t-l <filelist>                       fsa file list\n" + \
		"\t-c <sample list>                    sampleID as control\n" + \
		"\t-o <outprefix>                      output file name's prefix\n" + \
		"\t-h help"


###Standard DNA fragments and probe DNA fragments
ref_size=np.array([35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500])
probe_size=pd.Series([108, 122, 132, 150, 164, 181, 191, 201], index=['SMA-'+str(i) for i in [1,2,3,4,5,6,8,10]])

def fsa2df(fsa_file):
	new=SeqIO.AbiIO.AbiIterator(open(fsa_file, "rb"))
	record=new.next()
	num_dyes=record.annotations['abif_raw']['Dye#1']

	new_df=pd.DataFrame()
	#new_df.columns=['FAM', 'LIZ', 'NED', 'PET', 'VIC']	
	for i in range(1, num_dyes+1):
		n=i
		dye_name=record.annotations['abif_raw']['DyeN'+str(n)]
		if dye_name=='6-FAM':
			dye_name='FAM'
		if n==5:
			n=105
		raw_data=record.annotations['abif_raw']['DATA'+str(n)]
		new_df[dye_name]=raw_data
	return new_df

def fsa2df_analyzed(fsa_file):
	new=SeqIO.AbiIO.AbiIterator(open(fsa_file, "rb"))
	record=new.next()
	num_dyes=record.annotations['abif_raw']['Dye#1']

	new_df=pd.DataFrame()
	#new_df.columns=['FAM', 'LIZ', 'NED', 'PET']	
	for i in [1,2,3,4]:
		n=i+8
		dye_name=record.annotations['abif_raw']['DyeN'+str(i)]
		if dye_name=='6-FAM':
			dye_name='FAM'
		raw_data=record.annotations['abif_raw']['DATA'+str(n)]
		new_df[dye_name]=raw_data
	return new_df

	
def findPeakElement(nums):
	"""
	:type nums: Series
	:rtype: int
	"""
	N = len(nums)
	peak = np.array([])
	threshold=1000
	for i in range(len(nums)):
		if i==0 or i==N-1:
			continue
		else:
			if nums[i]>=nums[i-1] and nums[i]>=nums[i+1] and nums[i]>threshold:
				peak=np.append(peak, i)
				#print [i, nums[i]]
	true_peak=np.array([])
	smallest_space=50
	
	for i in range(len(peak)-1):
		if len(true_peak)>0 and peak[i]-true_peak[-1]<=smallest_space:
			continue
		if peak[i+1]-peak[i]>smallest_space:
			true_peak=np.append(true_peak, peak[i])
			#print("step1:\t%d\t%d" % (i, peak[i]))
		else:
			n=i+1
			while(n<len(peak)):
				if peak[n]-peak[i]>smallest_space:
					true_peak=np.append(true_peak, peak[i])
					#print("step2:\t%d\t%d" % (i, peak[i]))
					break
				elif nums[peak[n]]>=nums[peak[i]]:
					#print("n:%d\ti:%d\t%f" % (n, i, peak[n]))
					break
				elif n==len(peak)-1:
						true_peak=np.append(true_peak, peak[i])
						#print("step3:\t%d\t%d" % (i, peak[i]))
						break
				n+=1
	if len(true_peak)==0 or peak[-1]-true_peak[-1]>smallest_space:
		true_peak=np.append(true_peak, peak[-1])
		#print("step4:\t%d\t%d" % (-1, peak[-1]))
	return true_peak

def Peakmap(new_df, index2size, case_index2size, outfile, case_dye='FAM', ref_dye='PET', ref_size=ref_size, probe_size=probe_size):
	fig=plt.figure(figsize=(20, 6))
	subfig=fig.add_subplot(111)
#	plt.xlim([0, 1])
#	plt.ylim([0, 2.5])
	subfig.plot(new_df[case_dye], color='blue', linestyle='-')
	subfig.plot(new_df[ref_dye], color='red', linestyle='-')
	for k,v in index2size.iteritems():
		subfig.plot([k, k], [0, new_df[ref_dye][k]], color='red', linestyle='-', linewidth=0.1)
		subfig.text(k, new_df[ref_dye][k], str(round(v, 1)), color='red')
	for k,v in case_index2size.iteritems():
		subfig.plot([k, k], [0, new_df[case_dye][k]], color='blue', linestyle='-', linewidth=0.1)
		subfig.text(k, new_df[case_dye][k], str(round(v, 1)), color='blue')
	plt.savefig(outfile)
	plt.close()	


def getInternalRef(peak_index, ref_size):
	ref_intervals=ref_size[1:]-ref_size[0:-1]
	space=((peak_index[1:]-peak_index[0:-1])[-5:]/ref_intervals[-5:]).mean()
	index2size=dict()
	ref_index=[peak_index[-1]]
	index2size[peak_index[-1]]=ref_size[-1]
	n=2
	if len(ref_size)>len(peak_index):
		sys.exit("peak number ERROR\nref_size length:%d\tpeak num:%d" % (len(ref_size), len(peak_index)))
	for i in range(1,len(peak_index)+1):
		if n>len(ref_size):
			break
		expected_pos=ref_index[-1]-space*(ref_size[-(n-1)]-ref_size[-n])
		if abs(peak_index[-i]-expected_pos)<abs(peak_index[-(i+1)]-expected_pos):
			ref_index.append(peak_index[-i])			
			index2size[peak_index[-i]]=ref_size[-n]
			#print [peak_index[-i], ref_size[-n]]
			n=n+1
	return index2size

def getsize(peak_index, index2size):
	index=index2size.keys()
	index.sort()
	case_index2size=dict()
	for i in peak_index:
		for j in range(len(index)):
			if i==index[j]:
				size=index2size[index[j]]
				case_index2size[i]=size
				#print [i, size]
			elif i<index[j] and i>index[j-1]:
				size=index2size[index[j-1]]+(i-index[j-1])*(index2size[index[j]]-index2size[index[j-1]])*1.0/(index[j]-index[j-1])
				case_index2size[i]=size
				#print [i, size]
	return case_index2size


def get_probe_height(new_df, case_dye='FAM', ref_dye='PET', ref_size=ref_size, probe_size=probe_size):
	trace_case=new_df[case_dye]
	trace_in_ref=new_df[ref_dye]
	index2size=getInternalRef(findPeakElement(trace_in_ref), ref_size)
	case_index2size=getsize(findPeakElement(trace_case), index2size)
	probe_height=pd.Series(index=probe_size.index)
	for i in probe_size.index:
		probe_index=sorted(case_index2size.keys(), key=lambda a: abs(case_index2size[a]-probe_size[i]))[0]
		if abs(case_index2size[probe_index]-probe_size[i])>5:
			probe_height[i]=0
			continue
		probe_height[i]=trace_case[probe_index]
	return probe_height, index2size, case_index2size

def DQ_compute(MLPA_df, ref_sample=[], ref_probe_name=['SMA-1', 'SMA-2', 'SMA-8', 'SMA-10']):
	if len(ref_sample)==0:
		ref_sample=MLPA_df.columns
	###Intra-normalisation
	intra_sample_df=pd.DataFrame(index=MLPA_df.index, columns=MLPA_df.columns)
	for probe in MLPA_df.index:
		intra_sample_df.loc[probe]=((1/MLPA_df.ix[ref_probe_name]).mul(MLPA_df.ix[probe], axis=1)).median()
	###Inter-normalisation
	inter_sample_df=pd.DataFrame(index=MLPA_df.index, columns=MLPA_df.columns)
	for sample in intra_sample_df.columns:
		inter_sample_df[sample]=((1/intra_sample_df[ref_sample]).mul(intra_sample_df[sample], axis=0)).mean(1)
	return inter_sample_df

def DQ_CI_compute(MLPA_df, ref_sample=[], ref_probe_name=['SMA-1', 'SMA-2', 'SMA-8', 'SMA-10']):
	if len(ref_sample)==0:
		ref_sample=MLPA_df.columns
	###Intra-normalisation
	DQ_df=pd.DataFrame(index=MLPA_df.index, columns=MLPA_df.columns)
	DQ_amplitude=pd.DataFrame(index=MLPA_df.index, columns=MLPA_df.columns)
	for sample in MLPA_df.columns:
		for probe in MLPA_df.index:
			DQ_df_ref=pd.DataFrame(index=ref_probe_name, columns=ref_sample)
			for h in ref_sample:
				for z in ref_probe_name:
					DQ_df_ref[h][z]=(MLPA_df[sample][probe]/MLPA_df[sample][z])/(MLPA_df[h][probe]/MLPA_df[h][z])
			MAD=abs(DQ_df_ref-DQ_df_ref.median()).median().mean()
			DQ=DQ_df_ref.median().mean()
			DQ_df[sample][probe]=DQ
			sigma=((DQ_df_ref.median()-DQ)**2).mean()
			amplitude=1.96*((1.4826*MAD)**2+sigma**2)**0.5
			DQ_amplitude[sample][probe]=amplitude
	###compute ref probe confidence range
	amplitude=1.96*((DQ_df[ref_sample].sub(DQ_df[ref_sample].mean(1), axis=0))**2).mean(1)
	ref_CI_range=[DQ_df[ref_sample].mean(1), amplitude]
	return DQ_df, DQ_amplitude, ref_CI_range

def Report_visualization(sample, outfile, ref_probe_name=['SMA-1', 'SMA-2', 'SMA-8', 'SMA-10']):
	fig=plt.figure(figsize=(16, 8))
	subfig=fig.add_subplot(111)
	plt.xlim([-0.01, 1.01])
	plt.ylim([0, 2.5])
	num_probe=len(probe_size)
	each_width=1.0/num_probe
	
	ref_background_color='bisque'
	case_background_color='cornflowerblue'
	box_width=each_width*0.4
	box_color='lightgreen'
	ref_marker_color='orange'
	###plot the CI of reference samples
	case_line_color='black'
	case_marker_color='red'
	for i in range(num_probe):
		###plot background
		if probe_size.index[i] in ref_probe_name:
			subfig.add_patch(plt.Rectangle((i*each_width, 0), each_width, 2.5, facecolor=ref_background_color, edgecolor='white'))
		else:
			subfig.add_patch(plt.Rectangle((i*each_width, 0), each_width, 2.5, facecolor=case_background_color, edgecolor='white'))
		###plot the CI of reference samples
		mean=ref_CI_range[0].iloc[i]
		amplitude=ref_CI_range[1].iloc[i]
		centerline=(i+0.5)*each_width
		subfig.add_patch(plt.Rectangle((centerline-box_width/2, mean-amplitude), box_width, 2*amplitude, color=box_color))
		subfig.plot(centerline, mean, marker='x', color=ref_marker_color, markersize=10)
		###plot case sample's DQ
		mean=DQ_df[sample].iloc[i]
		amplitude=DQ_amplitude[sample].iloc[i]
		subfig.plot([centerline-box_width/2, centerline+box_width/2], [mean-amplitude, mean-amplitude], color=case_line_color)
		subfig.plot([centerline-box_width/2, centerline+box_width/2], [mean+amplitude, mean+amplitude], color=case_line_color)
		subfig.plot([centerline, centerline], [mean-amplitude, mean+amplitude], color=case_line_color)
		subfig.plot(centerline, mean, marker='o', color=case_marker_color, markersize=10)
	subfig.plot([0.5*each_width, 1-0.5*each_width], [0.7, 0.7], color='red')
	subfig.plot([0.5*each_width, 1-0.5*each_width], [1.3, 1.3], color='blue')
	subfig.set_xticks([(i+0.5)*each_width for i in range(num_probe)]+[-0.01, 1.01])
	subfig.set_xticklabels(list(probe_size.index)+['', ''])
	ax=plt.gca()
	ax.spines['right'].set_color('none')
	ax.spines['top'].set_color('none')
	plt.savefig(outfile)
	plt.close()
	
if __name__=='__main__':
	opts, args = getopt.getopt(sys.argv[1:], "hl:")
	fsalist=''
	for op, value in opts:
		if op=='-l':
			fsalist=value
		if op=='-h':
			usage()
			sys.exit()
	if fsalist=='':
		usage()
		sys.exit()
	MLPA_df=pd.DataFrame(index=probe_size.index)
	ref_sample=[]
	all_sample=[]
	for line in open(fsalist, 'r'):
		line=line.strip()
		sample, tag, fsa_file=line.split('\t')
		all_sample.append(sample)
		if tag=="r":
			ref_sample.append(sample)
		1 if os.path.exists(fsa_file) else sys.exit("Can't find %s" % fsa_file)
		new_df=fsa2df(fsa_file)
		#new_df=fsa2df_analyzed(fsa_file)
		probe_height, index2size, case_index2size=get_probe_height(new_df)
		Peakmap(new_df, index2size, case_index2size, sample+'.fsa.Peakmap.png')
		probe_height.name=sample
		MLPA_df=pd.concat([MLPA_df, probe_height], axis=1)
	#print MLPA_df	
	DQ_df, DQ_amplitude, ref_CI_range=DQ_CI_compute(MLPA_df, ref_sample)
	#DQ_df.to_csv('test_DQ.csv', sep=',', index=True, header=True)
	for sample in all_sample:
		Report_visualization(sample, sample+'_Report.png')
		



