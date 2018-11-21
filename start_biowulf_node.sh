#!/bin/bash

# makes an executable dir for scripts if does not exist
mkdir -p $HOME/bin

# add executable dir to path if not 
if [[ ":$PATH:" == *":$HOME/bin:"* ]]; then
echo "Your path is correctly set"
else
echo 'export PATH=$PATH:$HOME/bin' >> $HOME/.bash_profile
source ~/.bash_profile
echo "Your path was missing ~/bin, it has been added to your .bash_profile file."
fi

