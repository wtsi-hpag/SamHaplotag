# SamHaplotag
Process haplotag barcodes in SAM format. Converts BC/QT barcode SAM tags into haplotag RX, QX and BX tags. See https://github.com/evolgenomics/haplotagging or https://www.fml.tuebingen.mpg.de/chan-group/haplotagging/?L=1 for details on haplotagging.

Also comes with a tool '10xSpoof' for converting haplotag barcodes into 10x compatible barcodes.

# Example Usage
```bash
samtools view -h@ 16 -F 0xF00 reads.cram | SamHaplotag | samtools view -@ 16 -o tagged.cram
samtools fastq -@ 16 -nT BX -0 /dev/null -s /dev/null tagged_reads.cram | 10xSpoof SamHaplotag_Clear_BC | bgzip -@ 16 >10x_spoofed_reads.fq.gz
```

# Installation
Requires:
* clang >= 11.0.0
* meson >= 0.57.1
```bash
env CXX=clang meson setup --buildtype=release --prefix=<installation prefix> builddir
cd builddir
meson compile
meson test
meson install
```
note that [Samtools](http://www.htslib.org/) must be on your PATH for the tests to run
