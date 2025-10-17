!/usr/bin/bash
# Script: raw_fastq.sh
# Description: Download FASTQ files for SRRs listed in srrs.txt

mkdir -p raw_fastq

while read SRR; do
  echo "=== Processing $SRR ==="
  curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$SRR&result=read_run&fields=run_accession,fastq_ftp,library_layout" \
  | tail -n +2 | while IFS=$'\t' read -r run_accession fastq_ftp library_layout; do
      for url in $(echo $fastq_ftp | tr ';' ' '); do
          fname=$(basename $url)
          echo "Downloading $fname ($library_layout)"
          curl -L "ftp://$url" -o "raw_fastq/$fname"
      done
  done
done < srrs.txt

echo "=== All downloads completed! ==="
