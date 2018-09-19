#import MLPA_analysis
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
					#print MLPA_df[h][z]
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
	MLPA_df=pd.DataFrame(index=probe_size.index)
	new_df=pd.read_csv('area_compare.data.txt', sep='\t', header=0, index_col=['sample', 'probe'], dtype={'commercial':float, 'own':float})
	commercial_df=pd.DataFrame(index=new_df.index.get_level_values('probe').drop_duplicates())
	own_df=pd.DataFrame(index=new_df.index.get_level_values('probe').drop_duplicates())
	for sample in new_df.index.get_level_values('sample').drop_duplicates():
		sample_df=new_df.xs(sample, level='sample')
		commercial_df[sample]=sample_df['commercial']
		own_df[sample]=sample_df['own']
	print own_df
	print commercial_df
	ref_sample=['17D3245681', '18D0055121']
	DQ_df, DQ_amplitude, ref_CI_range=DQ_CI_compute(commercial_df, ref_sample)
	DQ_df.to_csv('commercial_area_DQ.csv', sep=',', index=True, header=True)
	for sample in commercial_df.columns:
		Report_visualization(sample, sample+'_area__commercial_Report.png')
	DQ_df, DQ_amplitude, ref_CI_range=DQ_CI_compute(own_df, ref_sample)
	DQ_df.to_csv('own_area_DQ.csv', sep=',', index=True, header=True)		
	for sample in own_df.columns:
		Report_visualization(sample, sample+'_area_own_Report.png')



