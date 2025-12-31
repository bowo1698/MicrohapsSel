#!/usr/bin/env python
# generates haplotype block genotypes
# by Dzianis Prakapenka
from __future__ import print_function
import sys
import os
import argparse

def make_arg_parser():
    app_name="get-hap-geno.py"
    description="creates haplotype genotype files based on block\
                info files"
    parser = argparse.ArgumentParser(prog=app_name,
            description=description,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--hap",
            default=argparse.SUPPRESS,
            nargs='+',
            #type=list,
            required=True,
            help="path(s) to hap chr file(s) [required]")
    parser.add_argument("--hapinfo",
            default=argparse.SUPPRESS,
            nargs='+',
            #type=list,
            required=True,
            help="path(s) to hap_info file(s) [required]")
    parser.add_argument("-m","--missing",
            type=str,
            help="coding for missing alleles, blocks with the allele will be set to missing-value")
    parser.add_argument("--missing-value",
            type=str,
            default="-9999",
            help="value to use for blocks with missing alleles")
    parser.add_argument("-o","--output",
            default='hap_geno',
            help="output folder for haplotype genotypes")
    parser.add_argument("-p","--output-prefix",
            default='hap_geno',
            help="output prefix for haplotype genotype files")
    parser.add_argument("--nosort",
            action="store_true",
            default=False,
            help="do not sort input files by number in filename")
    parser.add_argument("--header",
            action="store_true",
            default=False,
            help="set if there is a header in the hap input files")
    parser.add_argument("-V","--verbose",
            action="store_true",
            default=False,
            help="verbose output")
    return parser

def sort_files(tosort):
    # sorts list of strings by digits in the string
    if args.nosort:
        return tosort
    else:
        return (sorted(tosort,key=lambda a: get_digits(a)))
        #return (sorted(tosort,key=lambda a: int(''.join(list(filter(str.isdigit, a))))))

def get_digits(s):
    # get digits from a string
    return int(''.join(list(filter(str.isdigit, s))))

if __name__ == '__main__':
    parser = make_arg_parser()
    args = parser.parse_args()

    hap_files=sort_files(args.hap)
    hap_info_files=sort_files(args.hapinfo)

    outdir = args.output if os.path.isabs(args.output) else os.path.join(os.getcwd(), args.output)
    outdir = outdir + "/" if not outdir.endswith("/") else outdir
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if (len(hap_files) != len(hap_info_files)):
        sys.exit("ERROR: number of hap files and hap_info files do not match")

    fn = 0
    for hap_file,hap_info_file in zip(hap_files,hap_info_files):
        tmp_str = hap_file.split('\\')
        tmp_str = tmp_str[-1].split('/')
        hap_cnum = get_digits(tmp_str[-1])

        tmp_str = hap_info_file.split('\\')
        tmp_str = tmp_str[-1].split('/')
        hap_info_cnum = get_digits(tmp_str[-1])
        if (hap_cnum == hap_info_cnum):
            chr_name = str(hap_cnum)
        else:
            fn = fn+1
            chr_name = '_u_'+str(fn)
            print('WARNING, not sure which chromosome is', hap_file, hap_info_file)
            print('assigning the name', chr_name, 'fix it before proceeding')

        outfile=outdir+args.output_prefix+'_'+chr_name

        if args.verbose:
            print("reading", hap_info_file)
        try:
            hi=open(hap_info_file,'r')
        except IOError:
            print("can't open ", hap_info_file)

        blocks = [line.strip().split()[1:] for line in hi]
        hi.close()

        if args.verbose:
            print("reading", hap_file)
        try:
            h=open(hap_file,'r')
        except IOError:
            print("can't open ", hap_file)

        if args.header:
            h.readline()
        haps = [line.strip().split() for line in h]
        h.close()

        indlist=[]
        codes=[]
        for blk in blocks:
            label = 1
            hap_dict={}
            for hh in haps:
                if len(hh)==2:
                    h=hh[1]
                else:
                    h=''.join(hh[1:])
                h_part=h[int(blk[0])*2:int(blk[-1])*2+2]
                parents = [h_part[::2],h_part[1::2]]
                for p in parents:
                    if args.missing is not None and args.missing in p:
                        hap_dict[p] = args.missing_value
                    elif p not in hap_dict:
                        hap_dict[p] = str(label)
                        label = label + 1
            codes.append(hap_dict)

        if args.verbose:
            print("writing", outfile)
        try:
            o=open(outfile,'w')
        except IOError:
            print("can't open ", outfile)

        print('\t'.join(["ID"] + ["hap_{0}_{1}\thap_{0}_{1}".format(chr_name,b) for b in range(1,len(blocks)+1)]), file=o)
        for hh in haps:
            if len(hh)==2:
                ind,h=hh
            else:
                ind = hh[0]
                h=''.join(hh[1:])
            output = [ind]
            for c,blk in enumerate(blocks):
                h_part=h[int(blk[0])*2:int(blk[-1])*2+2]
                parents = [h_part[::2],h_part[1::2]]
                for p in parents:
                    output.append(codes[c][p])
            print('\t'.join(output), file=o)

        o.close()
