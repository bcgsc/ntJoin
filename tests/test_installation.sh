#!/bin/bash

# Run this script to test your ntJoin installation. The script expects ntJoin to be in your PATH, as well as all the dependencies installed (see README for details)

set -eux

echo "Running ntJoin..."

ntJoin  assemble target=scaf.f-f.fa target_weight=1 references='ref.fa' reference_weights='2' prefix=f-f_test k=32 w=1000

echo "Done! Compare your output files to those in the 'expected_outputs' directory to ensure the run was successful."
