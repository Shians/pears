git clone https://github.com/suhrig/arriba
cd ./arriba && make
cd ..
git clone https://github.com/alexdobin/STAR
cd ./STAR/source && make
cd ../../
wget https://github.com/DavidsonGroup/flexiplex/releases/download/v1.02.5/flexiplex-v1.02.5.tar.gz
tar -xf flexiplex-v1.02.5.tar.gz
cd ./flexiplex && make
cd ../barcodes/ && gunzip *.gz
cd ..
wget https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2
tar -xf samtools-1.18.tar.bz2
cd samtools-1.18
make

