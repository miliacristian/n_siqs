#!/bin/bash
cd src
make criv_quad_crash_test DEBUG=1 NUM_OF_N_TO_FACTORIZE=1
./criv_quad_crash_test ./../numbers/numbers_biprime.txt
make criv_quad_crash_test_clean
cd ..