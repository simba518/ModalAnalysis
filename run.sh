#! /bin/bash

./bin/modal_analysis ./data/dino.stiffness  ./data/dino.mass 20

./bin/modal_analysis ./data/Alla_shelfa_bunny_0.9525_0.1379.stiffness  ./data/Alla_shelfa_bunny_0.9525_0.1379.mass 10 ./data/Alla_shelfa_bunny_0.9525_0.1379.off

./bin/modal_analysis ./data/bunny_0.8355_0.1000.stiffness  ./data/bunny_0.8355_0.1000.mass 10 ./data/bunny_0.8355_0.1000.off

./bin/modal_analysis ./data/bunny_0.8768_0.1001.stiffness  ./data/bunny_0.8768_0.1001.mass 10 ./data/bunny_0.8768_0.1001.off

./bin/modal_analysis ./data/bunny_0.92_0.100.stiffness  ./data/bunny_0.92_0.100.mass 10 ./data/bunny_0.92_0.100.off

