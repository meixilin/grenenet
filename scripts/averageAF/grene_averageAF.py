# -*- coding: utf-8 -*-
'''
Title: Average allele frequency in each site/plot/year combinations by the number_flowers sampled
Project: grenenet
Author: Meixi Lin
Date: 2023-01-17 11:11:05
Example usage:
python3 grene_averageAF.py
'''

################################################################################
## import packages
import os
import sys
import warnings
import pandas as pd
import numpy

################################################################################
## def functions
def find_plots(df: pd.DataFrame,ss: int,yy: int):
    outdf = df.query('site == @ss and year == @yy and usesample == True')
    outplots = outdf['plot'].unique().tolist() # cannot access without quoting
    return outplots

def find_samples(df: pd.DataFrame,ss: int,pp: int,yy: int):
    outdf = df.query('site == @ss and plot == @pp and year == @yy and usesample == True')
    outsamples = outdf.sampleid.values.tolist()
    outnflowers = outdf.flowerscollected.values.tolist()
    return outsamples, outnflowers

def find_freqfiles(mydir: str, sampleids: list, mychr: int):
    outl = []
    for sample in sampleids:
        fullpath='{0}/{1}.merged.bam.{2}.afSite'.format(mydir,sample,str(mychr))
        if os.path.exists(fullpath):
            outl.append(fullpath)
        else:
            warnings.warn("AF files for {0} not found".format(sample))
    return outl

def check_pos(afdt: pd.DataFrame):
    ncol = afdt.shape[1]
    poscheck = []
    for ii in range(0,ncol,2):
        poscheck.append(afdt.iloc[:,0].equals(afdt.iloc[:,ii]))
    return all(poscheck)

def load_freqfiles(myfiles: list):
    aflist = [pd.read_csv(ff) for ff in myfiles]
    afheaders = [os.path.basename(ff).split('.')[0] for ff in myfiles]
    afdt = pd.concat(aflist, axis=1)
    # sanity check if all the pos are the same
    if check_pos(afdt):
        afcalc = afdt.iloc[:,range(1,afdt.shape[1],2)]
        afcalc.columns = afheaders
        afpos = afdt.iloc[:,0].values.tolist() # only keep the first one since they are all the same
    else:
        raise ValueError('afSites have conflicting pos')
    return afcalc, afpos

# append the samples being written
def write_AFsamples(mysamples, ss, pp, yy, outf: str):
    myline = '\t'.join([str(ss),str(pp),str(yy),';'.join(mysamples)])+'\n'
    with open(outf, 'a') as ff:
        ff.write(myline)
    return None

# write the output table
def write_averageAF(aveafdt,ss,myplots,yy,outf):
    myheaders = ['pos']+['-'.join(['FH-'+str(ss), str(pp), str(yy)]) for pp in myplots]
    aveafdf = pd.DataFrame(data = aveafdt).transpose()
    aveafdf.columns = myheaders
    aveafdf['pos'] = aveafdf['pos'].astype('int64')
    aveafdf.to_csv(outf, index=False)

################################################################################
## def variables
samples_data=pd.read_csv('/central/groups/carnegie_poc/meixilin/grenenet/metadata/data/samples_data.csv')
afdir='/central/groups/carnegie_poc/lczech/grenephase1/hafpipe-231/frequencies'
outdir='/central/groups/carnegie_poc/meixilin/grenenet/analyses/data/averageAF'

################################################################################
## main
def main():
    # read parameters
    site=int(sys.argv[1])
    year=int(sys.argv[2])
    chrom=int(sys.argv[3]) # the chromosome to work on
    outaffile = '{0}/FH_site{1}_year{2}_chr{3}_averageAF.csv'.format(outdir, str(site), str(year), str(chrom))
    outafsamples = '{0}/stats/FH_site{1}_year{2}_chr{3}_averageAF.samplelist.tsv'.format(outdir, str(site), str(year), str(chrom))
    # get the available plots to iterate over for given site
    allplots = find_plots(samples_data,site,year)
    # if there are no matches
    if len(allplots) == 0:
        sys.stdout.write('FH_site{0}_year{1}_chr{2} has no afSite files'.format(str(site), str(year), str(chrom)))
        sys.exit(0)
    # start iteration
    aveafdata = [] # the final output file
    for myplot in allplots:
        sampleids, nflowers = find_samples(samples_data,site,myplot,year) # organized in the sample order
        flowerdict = dict(zip(sampleids,nflowers))
        # find frequency files
        afpaths = find_freqfiles(afdir,sampleids,chrom)
        # read all the afpaths and check pos
        afdata, afpos = load_freqfiles(afpaths)

        # add the pos for the first plot
        if len(aveafdata) == 0:
            aveafdata.append(afpos)
        # check the pos for the other plots
        else:
            afpos0 = numpy.array(aveafdata[0])
            if not numpy.array_equal(afpos0, numpy.array(afpos)):
                raise ValueError('afSites have conflicting pos')

        # get samples that were being used in the end
        afsamples = afdata.columns.tolist()
        afflowers = [flowerdict[sample] for sample in afsamples]
        # output the afsamples that were used
        write_AFsamples(afsamples,site,myplot,year,outafsamples)

        # IMPORTANT: calculate weighted average for each site (a list of values in myplot)
        aflist = afdata.values.tolist()
        afout = [numpy.average(row, weights=afflowers) for row in aflist]
        aveafdata.append(afout)

    # output aveafdata
    write_averageAF(aveafdata,site,allplots,year,outaffile)

if __name__ == "__main__":
    sys.exit(main())
