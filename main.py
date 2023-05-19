import os
import sys
import numpy as np
import re
import argparse
import subprocess
import math
def convertToNumber (s):
    return int.from_bytes(s.encode(), 'little')
class PafEntry:
    def __init__(self, tabs, tags,fromstr):

        self.qr_name = tabs[0]
        self.qr_len = int(tabs[1])
        self.is_mapped = tabs[4] != ("*" if fromstr else None)
        self.label = tabs[-1]

        if self.is_mapped:
            self.qr_st = int(tabs[2])
            self.qr_en = int(tabs[3])
            self.is_fwd = tabs[4] == ('+' if fromstr else True)
            self.rf_name = tabs[5]
            self.rf_len = int(tabs[6])
            self.rf_st = int(tabs[7])
            self.rf_en = int(tabs[8])
            self.match_num = int(tabs[9])
            self.aln_len = int(tabs[10])
            self.qual = int(tabs[11])
        else:
            self.qr_st = 1
            self.qr_en = self.qr_len
            self.is_fwd=self.rf_name=self.rf_len=self.rf_st=self.rf_en=self.match_num=self.aln_len=self.qual=None

        self.tag = tags

    def rev(self):
        return PafEntry( [self.rf_name, self.rf_len, self.rf_st, self.rf_en, self.is_fwd,
                          self.qr_name, self.qr_len, self.qr_st, self.qr_en, self.match_num,
                          self.aln_len, self.qual], self.tags )

    def qry_loc(self):
        return (self.qr_name, self.qr_st, self.qr_en)

    def ref_loc(self):
        return (self.rf_name, self.rf_st, self.rf_en)

    def ext_ref(self, ext=1.0):
        st_shift = int(self.qr_st*ext)
        en_shift = int((self.qr_len - self.qr_en)*ext)

        if self.is_fwd:
            return (max(1, self.rf_st - st_shift),
                    min(self.rf_len, self.rf_en + en_shift))
        else:
            return (max(1, self.rf_st - en_shift),
                    min(self.rf_len, self.rf_en + st_shift))


    def __str__(self):
        tagstr = "\t".join( (":".join([k,v[1],str(v[0])]) for k,v in self.tags.items()))
        if self.is_mapped:
            s = "%s\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s" % (
                 self.qr_name, self.qr_len, self.qr_st, self.qr_en,
                 '+' if self.is_fwd else '-', self.rf_name,
                 self.rf_len, self.rf_st, self.rf_en,
                 self.match_num, self.aln_len, self.qual, tagstr)
        else:
            s = "\t".join((self.qr_name,str(self.qr_len)) + ("*",)*10 + (tagstr,))

        return s

def len_paf(infile):
    infile = open(infile)
    count = 0
    for i in infile:count+=1
    return count

def parse_paf_single(infile,idx,cls):
    infile = open(infile)
    data=infile.readlines()
    data.sort()
    past = ""
    for d in data:
        fromstr = type(d) != list
        if fromstr:
            tabs = d.split()
        else:
            tabs = d
        if tabs[0] != past:
            yield PafEntry(tabs,idx%cls,fromstr)
            past = tabs[0]

def parse_paf(infiles,cls):
    for idx,infile in enumerate(infiles):
        infile = open(infile)
        past = ""
        #print("idx : %d",idx)
        for l in infile:
            fromstr = type(l) != list
            if fromstr:
                tabs = l.split()
            else:
                tabs = l
            if tabs[0] != past:
                yield PafEntry(tabs,idx%cls,fromstr)
                past = tabs[0]

def fasta_id(fasta):
    f=open(fasta,'r')
    lines=f.readlines()

    hre=re.compile('>(\S+)')
    id = []
    for line in lines:
            outh = hre.search(line)
            if outh:
                    id.append(outh.group(1))
    return id

def evaluation_BE(qry, fasta,record):
    #print(qry)
    for idx,q in enumerate(qry):
        #print(convertToNumber(q.qr_name))
        if q.is_mapped:
            leng = abs(q.qr_en - q.qr_st)
            if record[idx][1] < leng:
                record[idx][0] = find_pos(fasta,q.rf_name)
                record[idx][1] = leng
    return record

def evaluation(qry, fasta,result):
    tp = list()
    tn = list()
    fp = list()
    fn = list()
    for q in qry:
        l = int(q.label)
        if q.is_mapped:
            if q.rf_name in fasta[l]:
                result[l][l] += 1
                tp.append(q)
            else:
                fp.append(q)
                result[l][find_pos(fasta,q.rf_name)] +=1
        else:
            if q.tag == l:
                tn.append(q)
            else:
                fn.append(q)

    return tp, tn, fp, fn, result

def find_pos(fasta,name):
    for y,row in enumerate(fasta):
        try:
            row.index(name)
            pos = y
            break
        except ValueError:
            pass
    return pos


def add_opts(parser):
    #parser.add_argument("infile", type=str, help="PAF file output by UNCALLED")
    parser.add_argument("-p","--pafpath", type=str,required=True, help="PAF file output by UNCALLED")
    parser.add_argument("-f","--fastapath", type=str,required=True, help="fasta file")
    parser.add_argument("-c","--classes", type=int,required=True, help="num of classes")
    parser.add_argument("-n", "--max-reads", required=False, type=int, default=None, help="Will only look at first n reads if specified")
    parser.add_argument("-r", "--ref-paf", required=False, type=str, default=None, help="Reference PAF file. Will output percent true/false positives/negatives with respect to reference. Reads not mapped in reference PAF will be classified as NA.")
    parser.add_argument("-a", "--annotate", action='store_true', help="Should be used with --ref-paf. Will output an annotated version of the input with T/P F/P specified in an 'rf' tag")

def run(args):
    ### Variable
    statsout = sys.stderr if args.annotate else sys.stdout
    cls = args.classes

    #### File Path Check
    statsout.write("Classification Result\n")
    if isinstance(args.pafpath, str) and isinstance(args.fastapath,str) and os.path.exists(args.pafpath) and os.path.exists(args.fastapath):
        paf_path = args.pafpath
        fasta_path = args.fastapath
    else:
        raise FileNotFoundError

    ### Loading PAF&FASTA Files' Path
    files = os.listdir(paf_path)
    files.sort()
    paf_list = [os.path.join(paf_path,f) for f in files if os.path.isfile(os.path.join(paf_path,f))]
    files = os.listdir(fasta_path)
    files.sort()
    fasta_list = [os.path.join(fasta_path,f) for f in files if os.path.isfile(os.path.join(fasta_path,f))]
    n_class = len(fasta_list)
    fasta_ids = []
    for f in fasta_list:
        fasta_ids.append(fasta_id(f))
    print(fasta_ids)

    ### Paf Entry
    locs = [p for p in parse_paf(paf_list,cls)]
    num_mapped = sum([p.is_mapped for p in locs])

    result = np.zeros((n_class,n_class),dtype=np.int16)
    n = len(locs)
    ### Break Even ###
    record = np.array([0])
    for idx, paf in enumerate(paf_list):
        #print(idx)
        #print(paf)
        qries = [p for p in parse_paf_single(paf,idx%cls,cls)]
        if idx % cls == 0:
            ids = len_paf(paf)
            record = np.zeros((ids,2))
            evaluation_BE(qries,fasta_ids,record)
        else:
            if len(record) == 1: raise NotImplementedError
            evaluation_BE(qries,fasta_ids,record)
            if idx % cls == cls-1:
                row = idx//cls
                for clm in range(len(result[row])):
                    result[row][clm] = np.count_nonzero(record == clm,axis=0)[0]

    ### NO Break Even ###
    """
    tp, tn, fp, fn, result  = evaluation(locs,fasta_ids,result)
    ntp,ntn,nfp,nfn = map(len, [tp, tn, fp, fn])
    """
    statsout.write("Summary: %d reads, %d mapped (%.2f%%)\n\n" % (len(locs), num_mapped, 100*num_mapped/len(locs)))
    #statsout.write("     P     N\n")
    #statsout.write("T %6.2f %5.2f\n" % (100*ntp/n, 100*ntn/n))
    #statsout.write("F %6.2f %5.2f\n" % (100*(nfp)/n, 100*nfn/n))

    #statsout.write("     P     N\n")
    #statsout.write("T  %d  %d\n" % (ntp, ntn))
    #statsout.write("F  %d  %d\n" % (nfp, nfn))

    print(result)
    print(np.sum(result,axis=1))

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='evaluation of alignment tool')
    add_opts(parser)
    run(parser.parse_args())

