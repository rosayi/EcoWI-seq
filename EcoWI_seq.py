#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-

try:
    import os
    import sys
    import re
    import argparse
    import subprocess 
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    pd.options.mode.chained_assignment = None
except:
    print('Module error! Please check if os, sys, re, argparse, subprocess, pandas, matplotlib and seaborn are installed')

def remove_unmapped_reads(bam,chr,base):
    ''' remove unmapped reads (-F 4) and export all mapped reads to tmp.bam'''
    cmd = f'samtools index {bam}'
    subprocess.check_call(cmd,shell=True)
    cmd = f'samtools view -b -F 4 {bam} {chr} > {base}_{chr}.tmp.bam'
    subprocess.check_call(cmd,shell=True)

def fasta_parser(file):
    ''' read a fasta file, return a list of tuples (seqname,seq) 
    '''
    seq = []
    new = []
    with open(file,'r') as f:
        header = f.readline().strip()
        for line in f:
            if line.startswith('>'):
                new.append(tuple((header,''.join(seq))))
                header = line.strip()
                seq=[]
            else:
                seq.append(line.strip())
        new.append(tuple((header, ''.join(seq))))
    return new

def find_ptsites():
    ''' return a dict of dicts containing {contig_id:{pos: category(GAAC/GTTC)}} and export files (<contig>.PTsites) containing all the position of GTTC/GAAC in reference strand (separated by tab):
      position (1-based)  motif
            1             GAAC               e.g. 1 is the position of G
    '''   
    if '.fasta' in os.path.basename(args.genome_fasta):
        pass
    elif '.fa' in os.path.basename(args.genome_fasta):
        pass
    else:
        sys.exit('Invalid fasta!')

    fasta = fasta_parser(args.genome_fasta)
    dicts = {}
    for header,seq in fasta:
        prefix = header.split('>')[1].split(' ')[0]        # same as chr in bam file
        output = f'{prefix}.allPTsites'
        pattern = re.compile(r'(GAAC)|(GTTC)',re.I)
        matches = pattern.finditer(seq)
        sites = {}
        test = peek(matches)
        if test:
            f=open(output,'w')
            matches = pattern.finditer(seq)
            for each in matches:
                position = str(int(each.start()) + 1)
                cat = each.group()
                sites.setdefault(position,cat)
                print(f'{position}\t{each.group()}', file=f)
            f.close()
            dicts.setdefault(prefix,sites)
    return dicts

def peek(iterable):
    try:
        first = next(iterable)
    except StopIteration:
        return None
    else:
        return first

def depth(chr,basename):
    ''' calculate depth on each base (split + and - strand), return the filename of forward depth and a dictionary containing key: pos (1-based), int; and value: depth in reverse strand
    '''
    cmd = [f'samtools view -b -f 16 {basename}_{chr}.tmp.bam | samtools depth -a - > {basename}_{chr}.r.depth',f'samtools view -b -F 16 {basename}_{chr}.tmp.bam | samtools depth -a - > {basename}_{chr}.f.depth']
    for i in cmd:
        subprocess.check_call(i,shell=True)
    depth_minus = {}
    with open(f'{basename}_{chr}.r.depth','r') as f:
        for line in f:
            count = int(line.strip().split()[2])
            pos = int(line.strip().split()[1])
            depth_minus.setdefault(pos,count)
    return f'{basename}_{chr}.f.depth', depth_minus

def end_pos(chr,basename):
    cmd = [f'bedtools bamtobed -i {basename}_{chr}.tmp.bam > {basename}_{chr}.bed',
           f'grep "+" {basename}_{chr}.bed | cut -f2 - | sort -n | uniq -c > {basename}_{chr}.plus.count',
           f'grep "-" {basename}_{chr}.bed | cut -f3 - | sort -n | uniq -c > {basename}_{chr}.minus.count']
    for each in cmd:
        subprocess.check_call(each,shell=True)

    dict_plus = {}
    with open(f'{basename}_{chr}.plus.count','r') as f:           # read into dict_plus {pos:5'ends on plus count} e.g. {53:23,54:1,23:0, ...}, pos: 1-based, int
        for line in f:
            count = int(line.strip().split()[0])
            pos = int(line.strip().split()[1]) + 1
            dict_plus.setdefault(pos,count)

    dict_minus = {}
    with open(f'{basename}_{chr}.minus.count','r') as f:          # read into dict_minus {pos:5'ends on minus count} e.g. {53:23, 54:1,23:0, ...}, pos: 1-based, int
        for line in f:
            count = int(line.strip().split()[0])
            pos = int(line.strip().split()[1])
            dict_minus.setdefault(pos,count)
    
    dict = {}                                             
    for pos, plus in dict_plus.items():                     # combine two dictionaries into a new dictionary containing {pos: (5'ends on +, 5'ends on -)}
        minus = dict_minus.get(pos,0)
        dict.setdefault(pos, (plus, minus))
        dict_minus.pop(pos,0)

    for pos, minus in dict_minus.items():                   # don't forget leftover unique pos in dict_minus
        dict.setdefault(pos, (0, minus))

    return dict

def identify_feature(new_ptsites, pos):
    ''' input a position, identify if it is GAAC/GTTC cleavage site
        return feature, motif_pos (str, str)

        feature:
        GAAC       (motif_pos = position of the G)
        GTTC       (motif_pos = position of the G)
        GAAC;GTTC  (motif_pos = position of the G, separate by ;)
        other      (motif_pos = 0)

        EcoWI cleavage pattern:
        -8 -7 -6 -5 -4 -3 -2 -1 G A A C 4 5 6 7 8 9 10 11 12

        position -8: GAAC_l_minus
        position -6: GAAC_l_plus
        position +9: GAAC_r_minus
        position +11: GAAC_r_plus
    '''

    if str(pos - 11) in new_ptsites.keys():
        motif_pos = str(int(pos) - 11)
        feature = f'{new_ptsites[motif_pos]}'

        if str(int(pos) + 6) in new_ptsites.keys():
            new_pos = str(int(pos)+6)
            motif_pos += f';{new_pos}' 
            feature += f';{new_ptsites[new_pos]}'
        elif str(int(pos) + 8) in new_ptsites.keys():
            new_pos = str(int(pos)+8)
            motif_pos += f';{new_pos}'
            feature += f';{new_ptsites[new_pos]}'

    elif str(int(pos) - 9) in new_ptsites.keys():
        motif_pos = str(int(pos) - 9)
        feature = f'{new_ptsites[motif_pos]}'

        if str(int(pos) + 6) in new_ptsites.keys():
            new_pos = str(int(pos)+6)
            motif_pos += f';{new_pos}' 
            feature += f';{new_ptsites[new_pos]}'
        elif str(int(pos) + 8) in new_ptsites.keys():
            new_pos = str(int(pos)+8)
            motif_pos += f';{new_pos}'
            feature += f';{new_ptsites[new_pos]}'
        
    elif str(int(pos) + 6) in new_ptsites.keys():
        motif_pos = str(int(pos) + 6)
        feature = f'{new_ptsites[motif_pos]}'

    elif str(int(pos) + 8) in new_ptsites.keys():
        motif_pos = str(int(pos) + 8)
        feature = f'{new_ptsites[motif_pos]}'

    else:
        feature = 'other'
        motif_pos = '0'

    return feature, motif_pos

def writer(new_ptsites, output):
    with open(output,'w') as writer:
        with open(f_depth,'r') as f:
            for line in f:
                seq = line.strip().split()[0]
                pos = int(line.strip().split()[1])
                f_reads = int(line.strip().split()[2])
                r_reads = dict_r_depth[pos]
                feature, motif_pos = identify_feature(new_ptsites, pos)
                ends_plus, ends_minus = dict.get(pos,(0,0))
                new_line = f'{seq}\t{feature}\t{pos}\t{motif_pos}\t{ends_plus}\t{ends_minus}\t{f_reads}\t{r_reads}'
                print(new_line,file=writer)
    
def clear(chr, basename):
    cmd = f'rm {basename}_{chr}.tmp.bam {basename}_{chr}.*.depth {basename}_{chr}.bed {basename}_{chr}.plus.count {basename}_{chr}.minus.count {basename}_{chr}.tmp'
    subprocess.check_call(cmd,shell=True)

def calculate_modification(chr,base):
    '''
    modificatoin = (5' ends / median 5' ends within +/- bp window )
    median ends 0 denomicator is replaced by 1
    '''
    df = pd.read_csv(f'{base}_{chr}.tmp',sep='\t',header=None)
    df.columns = ['seq','feature','pos','motif_pos','ends_plus','ends_minus','reads_plus','reads_minus']
    df['total_reads'] = df['reads_plus'] + df['reads_minus']
    df['total_ends'] = df['ends_plus'] + df['ends_minus']
    df['group'] = 'other'
    df.loc[df['feature'].str.contains('G'),'group'] = 'GAAC/GTTC_CUT'

    df['median_plus'] = df.ends_plus.rolling(2*args.window).median().shift(-args.window)
    df.median_plus.fillna(1,inplace=True)
    df.median_plus.replace(0,1,inplace=True)
    df['median_ratio_plus'] = df['ends_plus']/df['median_plus']

    df['median_minus'] = df.ends_minus.rolling(2*args.window).median().shift(-args.window)
    df.median_minus.fillna(1,inplace=True)
    df.median_minus.replace(0,1,inplace=True)
    df['median_ratio_minus'] = df['ends_minus']/df['median_minus']

    return df

def export_score_file(chr,df):
    df['source'] = f'{base}.bam'
    df['strand'] = '+'
    df.to_csv(f'{base}_{chr}.modification',sep='\t',columns=['seq','source','feature','pos','pos','motif_pos','strand','ends_plus','ends_minus','reads_plus','reads_minus','median_ratio_plus'],header=False,index=False)

def plot(chr, df):
    '''plot modification score (y axis) vs total reads (x axis), split by category: GAAC/GTTC or other
    '''
    g = sns.FacetGrid(df,col='group',hue='group')
    g.map(plt.scatter,'total_reads','median_ratio_plus',s=3,alpha=0.3)
    g.add_legend()
    g.set_axis_labels(x_var='total_reads',y_var='modification_score')
    plt.savefig(f'{base}_{chr}.modification_score.png',dpi=300,bbox_inches='tight')
    plt.clf()

def filter_data(df,score_cutoff):
    # filter 1: total_reads > read_cutoff, default = 50
    filter_1 = df[df['total_reads']>args.read_cutoff]
    if filter_1.empty:
        return filter_1
    else:
    # filter 2: score > ratio_cutoff, default = 0
        filter_2 = filter_1[filter_1['median_ratio_plus']> score_cutoff] 
        if filter_2.empty:
            return filter_2
        else:
    # filter 3: rel_ratio at -2bp position > ratio_cutoff
            index_list = filter_2.index.values.tolist()
            for i in index_list:
                if i==0 or i==1:
                    filter_2.loc[i,'score_minus_2'] = 0
                else:
                    filter_2.loc[i,'score_minus_2'] = df.loc[i-2,f'median_ratio_minus']
            filter_3 = filter_2[filter_2['score_minus_2'] > score_cutoff]
            return filter_3

def collapse_data(new_ptsites,filter):
    ''' collapse sites with proximity +/- 2 bp
    '''
    filter.loc[filter['feature']=='other','motif_pos'] = filter['pos']
    other_index = filter[filter['feature']=='other'].index.values.tolist()
    for i in other_index:
        for j in [i-2,i-1,i+1,i+2]:
            if str(j) in new_ptsites.keys():
                new_pos = j
                new_feature = new_ptsites[str(j)]
                filter.loc[j,'motif_pos'] = new_pos
                filter.loc[j,'feature'] = new_feature
    return filter

def plot_ROC(df,each):
    '''plot score cutoff (x axis) vs number of modified sites (y axis)
    '''
    datapoint = []
    for i in [5,10,15,20,30,40,50,60,70,80,90,100,110,120,130,140]:
        tmp_filter = filter_data(df,i)
        if tmp_filter.empty:
            count = 0
        else:
            collapse = collapse_data(new_ptsites, tmp_filter)     
            uniq = collapse.drop_duplicates('motif_pos')
            count = len(uniq)
        datapoint.append((i,count))
    df = pd.DataFrame(datapoint, columns = ['cutoff','number of pt site'])
    sns.lineplot(x='cutoff',y='number of pt site',marker='o',data=df)
    plt.savefig(f'{base}_{each}.curve.png',dpi=300)       
    plt.clf()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input',dest='input_bam')
    parser.add_argument('-g','--genome',dest='genome_fasta')
    parser.add_argument('-r','--read_cutoff',dest='read_cutoff', type=int, default=50)
    parser.add_argument('-s','--score_cutoff', dest='ratio_cutoff', type=float, default=0)
    parser.add_argument('-w','--window',dest='window',type=int,default=50)
    parser.add_argument('--plotcurve',action='store_true')                       # add optional argument with no value
    args = parser.parse_args()

    dir_path = os.getcwd()
    base = os.path.basename(args.input_bam).split('.bam')[0]

    # get all pt sites in the genome
    ptsites = find_ptsites()    

    for each in ptsites.keys():
        new_ptsites = ptsites[each]
        remove_unmapped_reads(args.input_bam,each,base)
        f_depth, dict_r_depth = depth(each,base)
        dict = end_pos(each,base)
        writer(new_ptsites, f'{base}_{each}.tmp')

        # calculate modification_score: median_ratio
        df = calculate_modification(each, base)

        # export input.modification
        export_score_file(each, df)

        # plot modification score
        plot(each, df)    

        # clear tmp
        clear(each, base)

        # filter 
        filter = filter_data(df,args.ratio_cutoff)
        if filter.empty:
            print(f'contig:{each} No modification sites found!')
            continue

        # collapse +/-2bp
        collapse = collapse_data(new_ptsites, filter)

        # output
        uniq = collapse.drop_duplicates('motif_pos')
        uniq.to_csv(f'{base}_{each}.sites',sep='\t',columns=['motif_pos','feature'],header=False, index=False)

        if args.plotcurve:
            plot_ROC(df,each)