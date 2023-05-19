import os
import glob
import re
import subprocess

txt_dir = '/z/kiku/Basecaller/Agn_Multi/fastq/'
work_dir = '/z/kiku/Basecaller/Agn_Multi'
fasta_dir = '/z/kiku/Basecaller/Agn_Multi/fasta/'

def main():
    txt_files = []
    fasta_all = []
    out_all = []
    if os.path.exists(txt_dir) and os.path.exists(fasta_dir):
        for name in glob.glob(txt_dir+'*.fastq'):
            filename = os.path.basename(name)
            txt_files.append(filename)
            tmp = re.search(('(.*)\.fastq$'),filename).group(1)
            print(tmp)
            a = os.path.basename(glob.glob(fasta_dir+tmp+'*')[0]) if glob.glob(fasta_dir+tmp+'*') else None
            print(a)
            fasta_all.append(a)
            out_all.append(tmp)
    else:
        raise NotImplementedError
    txt_files.sort()
    fasta_all.sort()
    out_all.sort()

    os.chdir('/z/kiku/Basecaller/minimap2')
    time = 0
    for file,out in zip(txt_files,out_all):
        for idx,fasta in enumerate(fasta_all):
            output_str = subprocess.run(f'./minimap2 --secondary=no --paf-no-hit {fasta_dir}{fasta} {txt_dir}{file} > {work_dir}/paf_test/{out}_{idx}.paf',shell=True,text=True,capture_output=True)
            sec = re.search(('Real time: (.*) sec; CPU'),str(output_str)).group(1)
            print(sec)
            time += float(sec)
    
    print(time)
            
            #print(out.stdout)
        #print(tmp)
        #print(fasta)


if __name__ == "__main__":
    main()
