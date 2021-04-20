[![test](https://github.com/wtsi-hpag/SamHaplotag/actions/workflows/test.yml/badge.svg)](https://github.com/wtsi-hpag/SamHaplotag/actions/workflows/test.yml)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/samhaplotag/badges/installer/conda.svg)](https://conda.anaconda.org/bioconda)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/samhaplotag/badges/downloads.svg)](https://anaconda.org/bioconda/samhaplotag)
# SamHaplotag
Process haplotag barcodes in SAM format. Converts BC/QT barcode SAM tags into haplotag RX, QX and BX tags. See https://github.com/evolgenomics/haplotagging or https://www.fml.tuebingen.mpg.de/chan-group/haplotagging/?L=1 for details on haplotagging.

Also comes with a couple of tools:
* '10xSpoof' for converting haplotag barcodes into 10x compatible barcodes
* '16BaseBCGen' for converting haplotag barcodes into generic 16-base barcodes with 7-base joins (useful for passing to programs like [ema](https://github.com/arshajii/ema))

# Bioconda
SamHaplotag is available on [bioconda](https://bioconda.github.io/).<br/>
```sh
> conda install samhaplotag
```

# Example Usage
```bash
> samtools view -h@ 16 -F 0xF00 reads.cram | SamHaplotag | samtools view -@ 16 -o tagged_reads.cram
> samtools fastq -@ 16 -nT BX -0 /dev/null -s /dev/null tagged_reads.cram | 10xSpoof SamHaplotag_Clear_BC | bgzip -@ 16 >10x_spoofed_reads.fq.gz

> samtools fastq -@ 16 -nT BX -0 /dev/null -s /dev/null tagged_reads.cram | 16BaseBCGen | bgzip -@ 16 >16BaseBC_reads.fq.gz
> cut -f 2 HaploTag_to_16BaseBCs | tail -n +2 >16BaseBCs
```

# Installation
Requires:
* clang >= 11.0.0
* meson >= 0.57.1
```bash
> env CXX=clang meson setup --buildtype=release --prefix=<installation prefix> builddir
> cd builddir
> meson compile
> meson test
> meson install
```
