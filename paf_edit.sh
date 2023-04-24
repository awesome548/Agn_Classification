#!/bin/zsh

paf="/z/kiku/Basecaller/Agn_Multi/paf"
#dir="/z/kiku/Basecaller/Agn_Multi/paf_edit"
dir="/z/kiku/Basecaller/Agn_Multi/paf_minimap2"
for f in $(ls $paf/$1*);do
    echo $f
    out=$(basename $f)
    #cat $f | awk "BEGIN{OFS=\" \"}{\$15=\$15\" $2\";print}" > $dir/$out
    cat $f | awk "BEGIN{OFS=\" \"}{\$18=\$18\" $2\";print}" > $dir/$out
done