#!/bin/bash

cat > srrs.txt <<'EOF'
SRR12808527
SRR12808528
SRR12808529
SRR12808497
SRR12808498
SRR12808499
EOF


echo "=== Downloading FASTQ files ==="
while read SRR; do
    echo "Processing $SRR..."
    curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$SRR&result=read_run&fields=run_accession,fastq_ftp,library_layout" \
    | tail -n +2 | while IFS=$'\t' read -r run_accession fastq_ftp library_layout; do
        for url in $(echo $fastq_ftp | tr ';' ' '); do
            fname=$(basename $url)
            echo "Downloading $fname ($library_layout)"
            curl -L "ftp://$url" -o "$RAW_DIR/$fname"
        done
    done
done < "$BASE_DIR/srrs.txt"
echo "=== FASTQ download completed ==="
