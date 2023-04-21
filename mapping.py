import os
import glob
import re
import subprocess

txt_dir = '/z/kiku/Basecaller/Agn_Multi/tmp_out/'

def main():
    txt_files = []
    fasta_all = []
    out_all = []
    if os.path.exists(txt_dir):
        for name in glob.glob(txt_dir+'*.txt'):
            txt_files.append(os.path.split(name)[1])
    else:
        raise NotImplementedError
    #print(txt_files) 
    txt_files.sort()
    for file in txt_files:
        tmp = re.search(('(.*)fast5\.txt$'),file)
        fasta = tmp[1].split('_')
        fasta_all.append(fasta[0]+'.'+fasta[1])
        tmp = file.split('_')
        out_all.append(tmp[0]+'_'+tmp[1])

    for file,out in zip(txt_files,out_all):
        for idx,fasta in enumerate(fasta_all):
            print(f'uncalled map -t 16 {txt_dir}{fasta} {txt_dir}{file} > PAF/{out}_{idx}.paf')
            subprocess.run(f'uncalled map -t 16 {txt_dir}{fasta} {txt_dir}{file} > PAF/{out}_{idx}.paf',shell=True, text=True)
            #print(out.stdout)
        #print(tmp)
        #print(fasta)


if __name__ == "__main__":
    main()
