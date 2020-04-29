#!/bin/bash
cd src
make criv_quad
./criv_quad ./../numbers/number_biprime_to_factorize.txt
make criv_quad_clean
cd ..