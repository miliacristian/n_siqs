#!/bin/bash
cd src
make criv_quad_debug DEBUG=1
./criv_quad_debug ./../numbers/number_biprime_to_factorize.txt
make criv_quad_debug_clean
cd ..