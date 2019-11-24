#!/bin/bash
cd src
make criv_quad_debug
valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose --log-file=./../logs/valgrind-out.txt ./criv_quad_debug ./../numbers/number_biprime_test.txt
make criv_quad_debug_clean
cd ..