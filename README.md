# bluefin.killifish.analysis
This repository contains a number helper scripts and C++ code written during my Ph.D. work that was used to uncover the genomic basis for the red/yellow color polymorphism in bluefin killifish.

This project was a collaborative effort between the [Catchen](https://catchenlab.life.illinois.edu/) and [Fuller](https://beckyfullerlab.weebly.com/) labs at the University of Illinois Urbana-Champaign. 

The C++ code is largely rewrites of some of the python scripts to speed up some of the data processing steps. Many of these files were written to process the output files from [popoolation2](https://sourceforge.net/p/popoolation2/wiki/Main/), which can produce very large *\*.sync* files depending on the sliding window used.

Custom visualizations were also created using some of the python scripts using the [pycario library](https://pycairo.readthedocs.io/en/latest/). This allowed me to visualize things like sequencing depth and genome annotations in a single image.

The README.md will be periodically updated as the project moves forward, with a mature description available once the paper is written and submitted for review.
