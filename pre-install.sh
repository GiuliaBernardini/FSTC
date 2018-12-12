#! /bin/sh

tar -xvf sdsl-lite.tar.gz
cd sdsl-lite
./install.sh "$(pwd)"/libsdsl
mv libsdsl/ ..
cd ..

unzip rmqo-master.zip
cd rmqo-master
make
make install
mv ./librmqo ..
cd ..
