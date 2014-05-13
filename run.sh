#! /bin/bash

# ./bin/modal_analysis ./data/beam3_hex.stiffness  ./data/beam3_hex.mass 2 LARGEST

./bin/modal_analysis ./data/beam3_hex20.stiffness  ./data/beam3_hex20.mass 20 SMALLEST

# ./bin/modal_analysis ./data/beam3_tet.stiffness  ./data/beam3_tet.mass 20 ./data/beam3_tet.off 4