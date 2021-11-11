import os
import argparse

parser = argparse.ArgumentParser(description="""aDNA BAM Mapping Plots For Bacteria and Viruses-  Meriam Guellil  -  September 2021 v2.0""", epilog="""Outputs Coverage, Edit Distance and Misincorporation plots for Each BAM Header""")
parser.add_argument('-b',metavar='BAM file', dest='bamR', required=True, type=str, help='Indexed BAM file, which should ideally not be filtered based on MAPQ (required)')
parser.add_argument('-d',metavar='Misincorporation file', dest='deam', required=False, type=str, help='mapDamage2 misincorporation.txt output (optional)')
parser.add_argument('-o',metavar='Output File', dest='out', required=True, type=str, help='Output file with extension for desired format (e.g. pdf, svg, png) (required)')
parser.add_argument('-i',metavar='Headers for output', dest='headlist', required=False, type=str, nargs='+', help='Space seperated list of bam headers to filter for (optional)')
parser.add_argument('-q',metavar='desired MQ cutoff', dest='mqf', required=False, default="30", type=str, help='MQ threshold for quality filtered coverage plot (default: 30)')
args= parser.parse_args()

import datetime
import numpy as np
import pysam
import pysamstats
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from io import StringIO
from tqdm import tqdm

date = datetime.datetime.now().strftime("%d/%m/%Y")

#Load BAM
bam = pysam.AlignmentFile(args.bamR, "rb")
head = {n["SN"]:n["LN"] for n in bam.header['SQ']}

#LOAD MISINCORPORATION
if args.deam:
    deam_df = pd.read_csv(args.deam,header=1,delimiter='\t',encoding='utf-8',skiprows=2)

#filter if available
if args.headlist:
    head = {k: v for k, v in head.items() if k in args.headlist}

#LOAD DEPTH
depth = pd.read_csv(StringIO(pysam.depth(args.bamR,"-aa")),delimiter='\t',encoding='utf-8',names=["id","pos","dp"])
depth30 = pd.read_csv(StringIO(pysam.depth(args.bamR,"-aa","-Q",args.mqf)),delimiter='\t',encoding='utf-8',names=["id","pos","dp"])
            
plt.rcParams.update({'font.size': 11, 'figure.titlesize' : 16, 'figure.titleweight' : "bold", 'legend.fontsize': 11})
fig = plt.figure(constrained_layout=False)
gs = gridspec.GridSpec(nrows=len(head), ncols=4, figure=fig)
fig.set_figheight(((3.5*len(head))))
fig.set_figwidth(25)
fig.suptitle('BAM File: ' + os.path.basename(args.bamR))

count = 0

pbar = tqdm(head.items(), desc="Progress", unit=" Plot", total=len(head))
for ID,LEN in pbar:
    #Empty Objects
    count += 1
    deamn3 = {}
    deamn5 = {}
    depthMQ30 = {}

    read_count = pysam.view("-c","-F","4",args.bamR,"-u",ID)
    read_count_mqf = pysam.view("-c","-F","4","-q",args.mqf,args.bamR,"-u",ID)

    if int(read_count) > 0:
        #Coverage Stats
        depthF = depth[depth['id'].values == ID]
        per_covered = (((depthF['dp'].values != 0).astype(int).sum(axis=0))/LEN)*100
        perc2x = (((depthF['dp'].values >= 2).astype(int).sum(axis=0))/LEN)*100
        meandepth = depthF["dp"].mean()
        ##Coverage Stats MQ30
        depthMQ30 = depth30[depth30['id'].values == ID]

        #DoC Stats
        Cov_df = depthF.drop(['id'], axis=1)
        Cov_Simpl = Cov_df.groupby(Cov_df.index // (LEN*0.002)).mean()
        ##DoC Stats MQ30
        Cov_df_MQ30 = depthMQ30.drop(['id'], axis=1)
        Cov_Simpl_MQ30 = Cov_df_MQ30.groupby(Cov_df_MQ30.index // (LEN*0.002)).mean()

        #MQ Stats            
        bamh = bam.fetch(ID)
        MQ_list = [read.mapping_quality for read in bamh]
        mean = np.average(MQ_list)
        #MQ Zero List
        mq = pysamstats.load_mapq_binned(bam, chrom=ID,window_size=(LEN*0.002))
        MQ = {i.pos:0 for i in mq if (i.rms_mapq==0 and i.reads_all!=0)}

        #ED Stats MQ0
        bamh = bam.fetch(ID)
        NM_list = [read.get_tag("NM") for read in bamh]
        ED_Dict = {x: NM_list.count(x)/len(NM_list) for x in range(0,6)}
        mean_ED = np.average(NM_list)
        ##ED Stats MQ30
        bamh = bam.fetch(ID)
        NM_list_MQ30 = [read.get_tag("NM") for read in bamh if read.mapping_quality >= int(args.mqf)]
        ED_Dict_MQ30 = {x: NM_list_MQ30.count(x)/len(NM_list) for x in range(0,6)}
        mean_ED_MQ30 = np.average(NM_list_MQ30)

        #Deamin Stats
        if args.deam:
            for i in range (1,6):
                df_5p = deam_df[(deam_df['Chr'].values == ID) & (deam_df['End'].values == "5p") & (deam_df['Pos'].values == i)]
                mean5 = df_5p['C>T'].sum()/df_5p['C'].sum()
                deamn5.update({i:mean5})
                df_3p = deam_df[(deam_df['Chr'].values == ID) & (deam_df['End'].values == "3p") & (deam_df['Pos'].values == i)]
                mean3 = df_3p['G>A'].sum()/df_3p['G'].sum()
                deamn3.update({i:mean3})

        #Subplot1
        f = fig.add_subplot(gs[(count-1), :3])  
        ##CoveragePlots
        plt.plot(Cov_Simpl_MQ30["pos"], Cov_Simpl_MQ30['dp'], color='#148F77',alpha = 1)
        f.fill_between(Cov_Simpl_MQ30["pos"], Cov_Simpl_MQ30['dp'], color='#148F77',alpha = 0.4)
        plt.plot(Cov_Simpl["pos"], Cov_Simpl['dp'], color='gray',alpha = 0.6)
        f.fill_between(Cov_Simpl["pos"], Cov_Simpl['dp'], color='gray',alpha = 0.2)
        f.axes = plt.gca()
        f.axes.set_ylim([0,10])
        f.axes.set_xlim([0,LEN])
        ##MQ Scatter
        ax2 = f.twinx()
        ax2.set_ylim([-1,1])
        ax2.get_yaxis().set_ticks([])
        plt.scatter(list(MQ.keys()),list(MQ.values()), color='#DB414D',alpha = 0.2)
        ##Legend Subplot 1
        f.text((LEN/100)*0.6, 9.8,ID, ha='left', va='top',weight="bold")
        f.text((LEN/100)*0.6, 9.2,'Mean Depth (MQ>=0): '+ str("{0:.2f}".format(meandepth)), ha='left', va='top')
        f.text((LEN/100)*0.6, 8.6,'Mapped Reads (MQ>=0): '+ str(f"{int(read_count):,}"), ha='left', va='top')
        f.text((LEN/100)*0.6, 8,'Mapped Reads (MQ>='+ args.mqf +'): '+ str(f"{int(read_count_mqf):,}"), ha='left', va='top')
        f.text((LEN/100)*0.6, 7.4,'Cov% (MQ>=0): '+ str("{0:.2f}".format(per_covered)), ha='left', va='top')
        f.text((LEN/100)*0.6, 6.8,'%2X (MQ>=0): '+ str("{0:.2f}".format(perc2x)), ha='left', va='top')
        f.text((LEN/100)*0.6, 6.2,'Mean MAPQ (>=0): '+ str("{0:.2f}".format(mean)), ha='left', va='top')
        f.axes.set_ylabel('')    
        f.axes.set_xlabel('')
        legend_elements1 = [Line2D([0], [0], marker='o', color='#148F77', label='CovMQ>='+args.mqf,markerfacecolor='#148F77', markersize=6),
                            Line2D([0], [0], marker='o', color='gray', label='CovMQ>=0',markerfacecolor='gray', markersize=6),
                            Line2D([0], [0], marker='o', color='#DB414D', label='MQ0', markerfacecolor='#DB414D', markersize=6,alpha = 0.3)]
        f.legend(handles=legend_elements1, bbox_to_anchor=(1,1), loc="upper right",fancybox=False, framealpha=0.4)

        #Subplot2
        ##ED
        e = fig.add_subplot(gs[(count-1),3]) 
        plt.bar(list(ED_Dict_MQ30.keys()), list(ED_Dict_MQ30.values()), color = "#148F77", alpha=0.5, width=0.8)
        plt.bar(list(ED_Dict.keys()), list(ED_Dict.values()),  alpha=0.3, color = "gray", width=0.5)
        ##Deamination
        if args.deam:  
            ax3 = e.twinx()
            plt.plot(list(deamn5.keys()),list(deamn5.values()), color = "#E15261", alpha = 0.7)
            plt.plot(list(deamn3.keys()), list(deamn3.values()), color = "#535F9F", alpha = 0.7)
            ##legend Subplot2
            e.text(5.6,max(list(ED_Dict.values())),'Mean ED (MQ>=0): '+ str("{0:.3f}".format(mean_ED)), ha='right', va='top')
            legend_elements2 = [Line2D([0], [0], marker='o', color="#148F77", label='EDistMQ>'+args.mqf, markerfacecolor='#148F77', markersize=6),
                                Line2D([0], [0], marker='o', color="gray", label='EDistMQ<'+args.mqf, markerfacecolor='gray', markersize=6, alpha=0.3),
                                Line2D([0], [0], marker='o', color="#E15261", label='5pCtoT', markerfacecolor='#E15261', markersize=6),
                                Line2D([0], [0], marker='o', color="#535F9F", label='3pGtoA', markerfacecolor='#535F9F', markersize=6)]
            e.legend(handles=legend_elements2,frameon=False, bbox_to_anchor=(1,0.9), loc="upper right")
        else:
            e.text(5.6,max(list(ED_Dict.values())),'Mean ED (MQ>=0): '+ str("{0:.3f}".format(mean_ED)), ha='right', va='top')
            legend_elements2 = [Line2D([0], [0], marker='o', color="#148F77", label='EDistMQ>'+args.mqf, markerfacecolor='#148F77', markersize=6),
                                Line2D([0], [0], marker='o', color="gray", label='EDistMQ<'+args.mqf, markerfacecolor='gray', markersize=6, alpha=0.3)]
            e.legend(handles=legend_elements2,frameon=False, bbox_to_anchor=(1,0.9), loc="upper right")
    else:
        f = fig.add_subplot(gs[(count-1), :3])
        f.axes = plt.gca()
        f.axes.set_ylim([0,10])
        f.axes.set_xlim([0,LEN])
        f.axes.set_ylabel('')
        f.axes.set_xlabel('')
        f.text((LEN/100)*1.5, 9.5,ID, ha='left', va='top',weight="bold")
        f.text((LEN/100)*20, 9.5,"NO MAPPING!", ha='left', va='top',weight="bold")
        e = fig.add_subplot(gs[(count-1),3])


#Date
fig.text(0.99, 0.99,date, ha='right', va='top',fontsize=9)

#Save Outfile
fig.tight_layout()
fig.savefig(args.out, bbox_inches='tight')


#Meriam Guellil 2021
