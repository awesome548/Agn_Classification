import os
import glob
import re
import subprocess

txt_dir = '/z/kiku/Basecaller/Agn_Multi/uncalled_misc/'

def main():
    """
    txt_files : fast5 address text
    fasta_all : Indexed file
    out_all : output file name
    """
    txt_files = []
    fasta_all = []
    out_all = []

    # fast5を取得
    if os.path.exists(txt_dir):
        for name in glob.glob(txt_dir+'*.txt'):
            txt_files.append(os.path.split(name)[1])
    else:
        raise NotImplementedError
    txt_files.sort()


    # Indexを取得 & 出力ファイル名の決定
    for file in txt_files:
        # 名前の部分を取得
        tmp = re.search(('(.*)fast5\.txt$'),file)
        fasta = tmp[1].split('_')
        # 整形
        fasta_all.append(fasta[0]+'.'+fasta[1])
        tmp = file.split('_')
        out_all.append(tmp[0]+'_'+tmp[1])

    for file,out in zip(txt_files,out_all):
        if file == "E_phage_fast5.txt":
            for idx,fasta in enumerate(fasta_all):
                print(f'uncalled map -t 16 {txt_dir}{fasta} {txt_dir}{file} > paf/{out}_{idx}.paf')
                subprocess.run(f'uncalled map -t 16 {txt_dir}{fasta} {txt_dir}{file} > paf/{out}_{idx}.paf',shell=True, text=True)
            


if __name__ == "__main__":
    main()
