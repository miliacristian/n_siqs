#!/bin/bash
cd src
make criv_quad_valgrind DEBUG=1
valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose --log-file=./../log_files/valgrind-out.txt ./criv_quad_valgrind ./../numbers/number_biprime_test_to_factorize.txt
make criv_quad_valgrind_clean
cd ..