#!/bin/bash
cd src
make criv_quad
./criv_quad ./../numbers/number_biprime.txt
make criv_quad_clean
cd ..