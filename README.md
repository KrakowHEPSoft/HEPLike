# HEPLike 
## Authors: Jihyun Bhom, Marcin Chrząszczn, jihyun.bohm@cern.ch, marcin.chrzaszcz@cern.ch 
### Free software, GPL 3
#### Computer framework to store and evaluate likelihoods coming from High Energy Physics experiments. Due to its flexibility it can be interfaced with existing fitting codes and allows to uniform the interpretation of the experimental results among the users. The code is provided with large open database, which contains the experimental measurements. The code is of use for users who perform phenomenological studies, global fits or experimental averages. 

#### Prerequirements:
git
cmake, 2.8
yaml-cpp, 1.58.0
gsl, 2.1
Boost, 1.58.0
ROOT, 6.08

In debian based systems:
sudo apt-get install cmake libgsl-dev libboost-all-dev libyaml-cpp-dev

To install ROOT see:
https://root.cern.ch/building-root

#### Instalation
git clone  https://github.com/mchrzasz/HEPLike.git

cd HEPLike

mkdir build

cd build

cmake ..

make -jN

#### Structure of the program:

main - directory containing the executable programs
data_toy -  directory with example datasets. THe real measurements are in seperate git repo (see below). 
include - directory with header files for HEPLike C++ classes.
src - directory with function definitions for HEPLike C++ clasees
utils - directory containing the python scripts for handing database, creating bibtex files etc.


##### To run examples:
ln -s data_toy data


##### Measurments repository:
https://github.com/mchrzasz/HEPLikeData

