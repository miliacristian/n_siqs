#!/bin/bash
cd src
make criv_quad_crash_test
./crash_test ./../numbers/numbers_biprime.txt
make criv_quad_crash_test_clean
cd ..