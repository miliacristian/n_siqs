#!/bin/bash
cd src/
make criv_quad_with_timing
./criv_quad_with_timing ./../numbers/number_biprime_to_factorize.txt
make criv_quad_with_timing_clean
cd ..