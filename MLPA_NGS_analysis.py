#import MLPA_analysis
import sys, getopt
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
plt.switch_backend('agg')

def usage():
	print 'Usage: python [options]' 
	print \
		"Options:\n" + \
		"\t-l <filelist>                       fsa file list\n" + \
		"\t-c <sample list>                    sampleID as control\n" + \
		"\t-o <outprefix>                      output file name's prefix\n" + \
		"\t-h help"


def DQ_CI_compute(MLPA_df, ref_sample=[], ref_probe_name=[]):
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

def Report_visualization(sample, outfile, probe_name, ref_probe_name=[]):
	fig=plt.figure(figsize=(16, 8))
	subfig=fig.add_subplot(111)
	plt.xlim([-0.01, 1.01])
	plt.ylim([0, 2.5])
	num_probe=len(probe_name)
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
		if probe_name[i] in ref_probe_name:
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
	subfig.set_xticklabels(list(probe_name)+['', ''])
	ax=plt.gca()
	ax.spines['right'].set_color('none')
	ax.spines['top'].set_color('none')
	plt.savefig(outfile)
	plt.close()
	
if __name__=='__main__':
	filelist=sys.argv[1]
	MLPA_df=pd.DataFrame()
	ref_sample=[]
	for line in open(filelist, 'r'):
		line=line.strip('\n')
		sample, tag, file=line.split('\t')
		if tag=='r':
			ref_sample.append(sample)
		new_df=pd.read_csv(file, sep='\t', header=None, index_col=0, dtype={1:float})
		new_df.columns=[sample]
		MLPA_df=pd.concat([MLPA_df, new_df[sample]], axis=1, sort=True)
	MLPA_df.index.name='probe'
	ref_probe_name=[i for i in MLPA_df.index if 'SMN' not in i]
	DQ_df, DQ_amplitude, ref_CI_range=DQ_CI_compute(MLPA_df, ref_sample, ref_probe_name)
	DQ_df.to_csv('NGS_DQ.csv', sep=',', index=True, header=True)
	for sample in MLPA_df.columns:
		Report_visualization(sample, sample+'_NGS_Report.png', MLPA_df.index, ref_probe_name)



