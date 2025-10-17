#!/usr/bin/bash
# download_reference.sh
# Download Arabidopsis thaliana genome and annotation (TAIR10) from Ensembl Plants
# and create STAR genome index

# === Directories ===
REF_DIR="reference"
GENOME_DIR="$REF_DIR/genome"
GTF_DIR="$REF_DIR/annotation"
STAR_INDEX_DIR="$REF_DIR/STAR_index"

mkdir -p "$GENOME_DIR" "$GTF_DIR" "$STAR_INDEX_DIR"

# === Download genome FASTA ===
echo "Downloading Arabidopsis thaliana genome (TAIR10) from Ensembl Plants..."
wget --no-check-certificate -c \
  https://ftp.ensemblgenomes.org/pub/plants/release-62/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz \
  -O "$GENOME_DIR/TAIR10.fa.gz"

# === Download annotation GFF3 ===
echo "Downloading Arabidopsis thaliana annotation (TAIR10) from Ensembl Plants..."
wget --no-check-certificate -c \
  https://ftp.ensemblgenomes.org/pub/plants/release-62/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.62.gff3.gz \
  -O "$GTF_DIR/TAIR10.gff3.gz"

# === Uncompress files ===
echo "Uncompressing genome and annotation files..."
gunzip -f "$GENOME_DIR/TAIR10.fa.gz"
gunzip -f "$GTF_DIR/TAIR10.gff3.gz"

# === Generate STAR genome index ===
echo "Generating STAR genome index..."
STAR --runThreadN 12 \
     --runMode genomeGenerate \
     --genomeDir "$STAR_INDEX_DIR" \
     --genomeFastaFiles "$GENOME_DIR/TAIR10.fa" \
     --sjdbGTFfile "$GTF_DIR/TAIR10.gff3" \
     --sjdbOverhang 99

echo "✅ Downloads and STAR genome index setup complete!"
echo "Genome: $GENOME_DIR/TAIR10.fa"
echo "Annotation: $GTF_DIR/TAIR10.gff3"
echo "STAR index: $STAR_INDEX_DIR"
