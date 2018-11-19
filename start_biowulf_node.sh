#!/bin/bash

if [[ ":$PATH:" == *":$HOME/bin:"* ]]; then
	echo "Your path is correctly set"
else
	echo 'export PATH=$PATH:$HOME/run' >> ~/.bash_profile
	source ~/.bash_profile
	echo "Your path was missing ~/bin, it has been added to your .bash_profile file."
fi
