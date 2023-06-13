# nanoev
Cost-Effective Complete Genome Sequencing Using the MinION Platform for Identification of Recombinant Enteroviruses

**To run with Docker**

``git clone https://github.com/jade-nhri/nanoev.git``

``cd nanoev``

``docker build -t "nanoev:v1" ./``

``docker run --runtime=nvidia -h nanoev --name nanoev -i -t -v /:/MyData nanoev:v1 /bin/bash``

Please note that you need to run “conda activate homopolish” before running nanoev.py.

Installation
------------
**Installation from source**

``cd /opt``

``git clone https://github.com/jade-nhri/nanoev.git``

``cd nanoev``

``chmod +x *.py``

``export PATH="$PATH:/opt/nanoev"``


## Dependencies

- [pyspoa-0.0.3](https://github.com/nanoporetech/pyspoa)
- [medaka-1.5.0](https://github.com/nanoporetech/medaka)
- [minimap2](https://github.com/lh3/minimap2)
- [samtools-1.13](http://github.com/samtools/)
- [bcftools](https://github.com/samtools/bcftools)
- [seqkit-v2.2.0](https://github.com/shenwei356/seqkit)
- [blast-2.13.0](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/)
- [miniconda](https://repo.anaconda.com/miniconda)
- [homopolish](https://github.com/ythuang0522/)
