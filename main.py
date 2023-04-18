import os
import sys
import numpy as np
import re
import argparse

class PafEntry:
    def __init__(self, line, tags=None):
        fromstr = type(line) != list
        if fromstr:
            tabs = line.split()
        else:
            tabs = line

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
        #self.tags = dict() if tags==None else tags 
        #for k,t,v in (s.split(":") for s in tabs[12:]):
            #if t == 'f':
                #v = float(v)
            #elif t == 'i':
                #v = int(v)
            #elif t not in ['A', 'Z','B', 'H']:
                #sys.stderr.write("Error: invalid tag type \"%s\"\n" % t)
                #sys.exit(1)
            #self.tags[k] = (v,t)

    def rev(self):
        return PafEntry( [self.rf_name, self.rf_len, self.rf_st, self.rf_en, self.is_fwd, 
                          self.qr_name, self.qr_len, self.qr_st, self.qr_en, self.match_num,
                          self.aln_len, self.qual], self.tags )

    def get_tag(self, k):
        return self.tags.get(k, (None,None))[0]
    
    def set_tag(self, k, v, t=None):
        if t == None:
            if isinstance(v, int):
                t = 'i'
            elif isinstance(v, float):
                t = 'f'
            else:
                t = 'Z'
        self.tags[k] = (v,t)

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

    def __lt__(self, paf2):
        return self.qr_name < self.qr_name

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


def parse_paf(infiles):
    for infile in infiles:
        infile = open(infile)
        for l in infile:
            #print(l)
            if l[0] == "#": continue
            yield PafEntry(l)

def paf_ref_compare(qry, ref, ret_qry=True, check_locs=True, ext=1.5):
    if type(ref) == dict:
        ref_locs = ref
    else:
        ref_locs = dict()
        for r in ref:
            l = ref_locs.get(r.qr_name, None)
            if l == None:
                ref_locs[r.qr_name] = [r]
            else:
                l.append(r)

    tp = list()
    tn = list()
    fp = list()
    fn = list()
    fp_unmap = list()

    for q in qry:
        rs = ref_locs.get(q.qr_name, [None])
        if q.is_mapped:
            if rs == [None] or not rs[0].is_mapped:
                fp_unmap.append(q if ret_qry else rs[0])
                continue
            match = False
            for r in rs:
                if ((check_locs and q.overlaps(r, ext)) or 
                    (not check_locs and q.rf_name == r.rf_name)):
                        match = True
                        tp.append(q if ret_qry else r)
                        break
            if not match:
                fp.append(q if ret_qry else r)
        else:
            if rs == [None] or not rs[0].is_mapped:
                tn.append(q if ret_qry else rs[0])
            else:
                fn.append(q if ret_qry else rs[0])

    return tp, tn, fp, fn, fp_unmap
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

def evaluation(qry, fasta,result):

    tp = list()
    tn = list()
    fp = list()
    fn = list()
    for q in qry:
        if q.is_mapped:
            l = int(q.label)
            if q.rf_name in fasta[l]:
                result[l][l] += 1
                tp.append(q)
            else:
                fp.append(q)
                result[l][find_pos(fasta,q.rf_name)] +=1
        else:
            if q.tag:
                fn.append(q)
            else:
                tn.append(q)

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
    parser.add_argument("infile", type=str, help="PAF file output by UNCALLED")
    parser.add_argument("-p","--pafpath", type=str, help="PAF file output by UNCALLED")
    parser.add_argument("-f","--fastapath", type=str, help="fasta file")
    parser.add_argument("-n", "--max-reads", required=False, type=int, default=None, help="Will only look at first n reads if specified")
    parser.add_argument("-r", "--ref-paf", required=False, type=str, default=None, help="Reference PAF file. Will output percent true/false positives/negatives with respect to reference. Reads not mapped in reference PAF will be classified as NA.")
    parser.add_argument("-a", "--annotate", action='store_true', help="Should be used with --ref-paf. Will output an annotated version of the input with T/P F/P specified in an 'rf' tag")

def run(args):
    # paf instance generated
    statsout = sys.stderr if args.annotate else sys.stdout

    statsout.write("Classification Result\n")
    if isinstance(args.pafpath, str) and isinstance(args.fastapath,str) and os.path.exists(args.pafpath) and os.path.exists(args.fastapath):
        paf_path = args.pafpath
        fasta_path = args.fastapath
    else:
        FileNotFoundError

    # loading paf&fasta files' path
    files = os.listdir(paf_path)
    paf_list = [os.path.join(paf_path,f) for f in files if os.path.isfile(os.path.join(paf_path,f))]
    files = os.listdir(fasta_path)
    fasta_list = [os.path.join(fasta_path,f) for f in files if os.path.isfile(os.path.join(fasta_path,f))]
    n_class = len(fasta_list)

    fasta_ids = []
    for f in fasta_list:
        fasta_ids.append(fasta_id(f)) 
    locs = [p for p in parse_paf(paf_list)]
    num_mapped = sum([p.is_mapped for p in locs])

    print(fasta_ids)
    print(len(fasta_ids))
    result = np.zeros((n_class,n_class))

    tp, tn, fp, fn, result  = evaluation(locs,fasta_ids,result)
    ntp,ntn,nfp,nfn = map(len, [tp, tn, fp, fn])
    n = len(locs) 
    statsout.write("Summary: %d reads, %d mapped (%.2f%%)\n\n" % (len(locs), num_mapped, 100*num_mapped/len(locs)))

    #statsout.write("     P     N\n")
    #statsout.write("T %6.2f %5.2f\n" % (100*ntp/n, 100*ntn/n))
    #statsout.write("F %6.2f %5.2f\n" % (100*(nfp)/n, 100*nfn/n))

    statsout.write("     P     N\n")
    statsout.write("T  %d  %d\n" % (ntp, ntn))
    statsout.write("F  %d  %d\n" % (nfp, nfn))

    print(result)

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='evaluation of alignment tool')
    add_opts(parser)
    run(parser.parse_args())

