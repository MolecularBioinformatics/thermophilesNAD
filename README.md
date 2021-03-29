# thermophilesNAD
A coarse-grained temperature dependent model of NAD metabolism via [PNCA](https://www.uniprot.org/uniprot/A0A2C8ER35) and [NAMPT](https://www.uniprot.org/uniprot/P43490).

The computational framework is developed using Python3.

****
#### 1. Getting python
You can download the latest version of Python [here](https://www.python.org/downloads/).
Alternatively, you can use [Anaconda](https://www.anaconda.com/download/) to install python for your computer (Linux, Windows, Mac).

#### 2. Cloning the repository
via SSH:

      git clone git@github.com:MolecularBioinformatics/thermophilesNAD.git

via HTTPS:

      git clone https://github.com/MolecularBioinformatics/thermophilesNAD.git
      
via GitHub CLI:

      gh repo clone MolecularBioinformatics/thermophilesNAD
  
#### 2. Installing packages using requirement file
      cd thermophiles
      pip install -r requirements.txt

****
Once you have checked the steps above, start a python console to simulate thesteady-state metabolic concentrations of NAD biosynthesis. 

#### 3. To reproduce plots from the manuscript type the following in your python console:
    run plotresults.ipynb
