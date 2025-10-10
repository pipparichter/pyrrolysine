DBNAME="arf1_cleanedDB"
NAME="arf1_cleaned"

rm ./data/mmseqs/databases/*
rm ./data/mmseqs/*
rm -f ./data/mmseqs/tmp/*

mmseqs createdb ./data/$NAME.fa ./data/mmseqs/databases/$DBNAME --dbtype 1

# Running with high sensitivity. 
mmseqs prefilter -s 7.5  ./data/mmseqs/databases/$DBNAME ./data/mmseqs/databases/$DBNAME ./data/mmseqs/databases/$DBNAME_prefilter

mmseqs align -a -e 100 --min-aln-len 10 --min-seq-id 0 ./data/mmseqs/databases/$DBNAME ./data/mmseqs/databases/$DBNAME ./data/mmseqs/databases/$DBNAME_prefilter ./data/mmseqs/$DBNAME_align

mmseqs convertalis --search-type 1 --format-output "query,target,evalue,pident,bits,qseq,tseq,alnlen,qstart,qend,tstart,tend,mismatch,gapopen" ./data/mmseqs/databases/$DBNAME ./data/mmseqs/databases/$DBNAME ./data/mmseqs/$DBNAME_align ./${NAME}_align.tsv
