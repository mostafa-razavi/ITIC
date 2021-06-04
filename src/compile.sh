#!/bin/bash

Dir=$PWD
rm -rf "$Dir/ITIC"
gfortran -o "$Dir/ITIC" "$Dir/ITIC.f90" "$Dir/LmdifEzCov2.for" "$Dir/MinTools.for"