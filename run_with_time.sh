#!/bin/bash
cd src
make criv_quad_with_timing
./criv_quad_crash_test ./../numbers/numbers_biprime.txt
make criv_quad_with_timing_clean
cd ..