#!/usr/bin/env python
# convert vcf phased/unphased format to findhap type
# by Dzianis Prakapenka
from __future__ import print_function
import os
import argparse
import math
import re

def make_arg_parser():
    app_name="convert-from-vcf.py"
    description="converts vcf to gvchap format"
    parser = argparse.ArgumentParser(prog=app_name,
            description=description,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i","--input",
            default=argparse.SUPPRESS,
            nargs='+',
            #type=list,
            required=True,
            help="full path(s) to vcf chr file(s)")
    parser.add_argument("-m", "--missing",
            type=str,
            default="-9999",
            help="code for missing genotype in snp chr files")
    parser.add_argument("--map",
            #default=argparse.SUPPRESS,
            required=False,
            dest='mapfile',
            default='map_new.txt',
            help="output path for map file formatted 'SNPID Chr Position'")
    parser.add_argument("--genofolder",
            default='geno',
            help="output path for snp genotype chr files")
    parser.add_argument("--hapfolder",
            default='hap',
            help="output path for haplotype chr files")
    parser.add_argument("--genoprefix",
            default="chr",
            help="prefix for the seperated geno files")
    parser.add_argument("--happrefix",
            default="chr",
            help="prefix for the seperated hap files")
    parser.add_argument("--nosort",
            action="store_false",
            default=True,
            help="do not sort vcf file by position")
    parser.add_argument("--interval",
            type=int,
            default="1",
            help="number of passes, reduces memory usage but is slower")
    parser.add_argument("-V","--verbose",
            action="store_true",
            default=False,
            help="verbose output")
    return parser

def sort_files(tosort):
    # sorts list of string by digits in the string
    try:
        if args.nosort:
            return (sorted(tosort,key=lambda a: int(''.join(list(filter(str.isdigit, a))))))
        else:
            return tosort
    except ValueError:
        print("WARNING: Unable to sort by number in input files. Use --nosort option")
        return tosort

def decode_snp(s):
    try:
        alleles = re.split('[|]|/',s)
        if len(alleles)==2:
            if ('.' in alleles):
                return args.missing
            elif (alleles[0]==alleles[1]):
                if alleles[0]=='1':
                    return '2'
                elif (alleles[0]=='0'):
                    return '0'
                else:
                    print("ERROR", s, 'setting as missing')
                    return args.missing
            else:
                return '1'
        else:
            print("ERROR, cannot split", s, 'setting as missing')
            return args.missing
    except ValueError:
        print("cannot recongnize:", s)
        return args.missing

def decode_hap(s):
    try:
        alleles = re.split('[|]',s)
        if len(alleles)==2:
            out = ''
            for a in alleles:
                if (a=='0'):
                    out = out + '1'
                elif (a=='1'):
                    out = out + '2'
                else:
                    out = out + args.missing
                out = out + '\t'
            return out.strip()
        else:
            try:
                if ('/' in s):
                    print("ERROR", s, 'setting as missing')
                    print("found slash separator, are you sure your vcf is phased?")
                else:
                    print("cannot recongnize:", s)
                return args.missing
            except ValueError:
                print(s, "is not valid")
                return args.missing
    except ValueError:
        print("ERROR, cannot split", s, 'setting as missing')
        return args.missing


if __name__ == '__main__':
    print('')
    parser = make_arg_parser()
    args = parser.parse_args()

    genodir = args.genofolder if args.genofolder.startswith('/') else os.path.join(os.getcwd(), args.genofolder)
    hapdir = args.hapfolder if args.hapfolder.startswith('/') else os.path.join(os.getcwd(), args.hapfolder)
    genodir += "/" if not genodir.endswith("/") else ""
    hapdir += "/" if not hapdir.endswith("/") else ""
    if not os.path.exists(genodir): 
        os.makedirs(genodir)
    if not os.path.exists(hapdir): 
        os.makedirs(hapdir)

    # sort input files by number in name
    vcf_files=sort_files(args.input)
    snp_dict={}

    # start of columns in vcf
    chrom_col = 0
    pos_col = 1
    snp_col = 2
    geno_start = 9

    map_data = {}


    for v in vcf_files:
        idx_begin = geno_start

        if args.verbose:
            print("reading", v)
        try:
            v_data=open(v,'r')
        except IOError:
            print("can't open ", v)

        while True:
            line = v_data.readline()
            if line.startswith('##'):
                #content_pos = v_data.tell()
                continue
            elif line.startswith('#'):
                content_pos = v_data.tell()
                l = line.strip().split()
                inds = l[geno_start:]
                num_ind = len(inds)
                chrom = []
                pos = []
                snp = []
            else:
                break

        map_data_local = []
        v_data.seek(content_pos)
        for line in v_data:
                l = line.strip().split()
                map_data_local.append([l[snp_col],l[chrom_col],l[pos_col]])
                map_data.setdefault(l[chrom_col], {})
                map_data[l[chrom_col]][int(l[pos_col])] = l[snp_col]

        chromosomes = set(list(zip(*map_data_local))[1])
        interval = int(math.floor(num_ind/args.interval))
        if args.verbose:
            print('.')
            print('\tindividuals:', num_ind)
            print('\tchromosomes:', len(chromosomes))
            print('\tinterval set to', interval, 'individuals')
            print('.')


        hap_chr_name = args.hapfolder + '/' + args.happrefix
        hap_chr_files = {}
        for c in chromosomes:
            hap_file_name = (hap_chr_name + c)
            if args.verbose:
                print("\topening ", hap_file_name, "for writing")
            try:
                hap_chr_files[c]=open(hap_file_name,'w')
            except IOError:
                print("\tcan't open ", hap_file_name)
            #print('ID\t' + '\t'.join([(m[0]+'_1\t'+m[0]+'_2') for m in map_data_local if m[1]==c]), file=hap_chr_files[c])

        geno_chr_name = args.genofolder + '/' + args.genoprefix
        geno_chr_files = {}
        for c in chromosomes:
            geno_file_name = (geno_chr_name + c)
            if args.verbose:
                print("\topening ", geno_file_name, "for writing")
            try:
                geno_chr_files[c]=open(geno_file_name,'w')
            except IOError:
                print("\tcan't open ", geno_file_name)
            print('ID\t' + '\t'.join([m[0] for m in map_data_local if m[1]==c]), file=geno_chr_files[c])


        output = []
        counter = 0
        if args.verbose:
            print('.')
            print('converting')
            print('.')
        for i in range(geno_start,num_ind+interval,interval):
            end = min(num_ind+geno_start,i+interval)
            if (i >= end): break
            counter += 1
            if args.verbose:
                print('\tpass', counter)

            #print('\tdebug', counter, geno_start, num_ind, interval,i,end)
            v_data.seek(content_pos)
            #end = i+interval
            ind_start = max(0,i-geno_start)
            ind_end = end-geno_start
            #print(i,end,ind_start, ind_end)
            #print(inds[ind_start:ind_end])
            #tmpind = [[None] * interval]
            ind_chr_dict = {}
            for c in chromosomes:
                ind_chr_dict[c] = [inds[ind_start:ind_end]]
            #tmpind = []
            for line in v_data:
                l = line.strip().split()
                c = l[chrom_col]
                part = l[i:end]
                #print(i,end,part,len(part))
                ind_chr_dict[c].append(part)
                #tmpind = list(zip(tmpind,part))
            for c in chromosomes:
                output = list(zip(*ind_chr_dict[c]))
                for o in output:
                    print(o[0]+'\t'+'\t'.join([decode_hap(h) for h in o[1:]]), file=hap_chr_files[c])
                    print(o[0]+'\t'+'\t'.join([decode_snp(h) for h in o[1:]]), file=geno_chr_files[c])

        for c in chromosomes:
            hap_chr_files[c].close()

        v_data.close()

    if args.verbose:
        print('.')
        print('done\n')

    if args.mapfile is not None:
        # open and write map file
        if args.verbose:
            print('.')
            print("writing to", args.mapfile)
        try:
            mp=open(args.mapfile, "w")
        except IOError:
            print("can't open ", args.mapfile)

        print('\t'.join(['SNPID','Chr','Position']), file=mp)
        for m in sorted(map_data.keys()):
            for p in sorted(map_data[m].keys()):
                print('\t'.join([map_data[m][p],m,str(p)]), file=mp)
        mp.close()

    if args.verbose:
        print('.')
        print('done\n')
