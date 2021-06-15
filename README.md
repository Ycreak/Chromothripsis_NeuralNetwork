# Chromothripsis_NeuralNetwork

##To create a dataset, run the create_trainingset.py code:

$ python3 create_trainingset.py

Make sure that the (extracted) CN data is located in ./dataset/cn/ 
https://dcc.icgc.org/releases/PCAWG/consensus_cnv (last visited 2020-06-15)

Make sure that the (extracted) SV data is located in ./dataset/sv/ 
https://dcc.icgc.org/releases/PCAWG/consensus_sv (last visited 2020-06-15)

This code calls the detect_chromothripsis.R file. This uses the ShatterSeek library. Via the Pacman package manager, all requirements should be installed automatically.

##To create the neural network, run the create_neuralnetwork.py code

All dependencies can be installed using:

$ pip3 install -r requirements.txt

