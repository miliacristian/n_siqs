#!/bin/bash
cd src
make criv_quad_with_timing
./crash_test ./../numbers/numbers_biprime.txt
make criv_quad_with_timing_clean
cd ..