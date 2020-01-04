#!/bin/bash
make criv_quad DEBUG=1
./criv_quad ./../numbers/number_biprime.txt
make criv_quad_clean