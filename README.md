# TriangulaMarkovCoLocalizer (TMCL)

A-star 10-X genomics spatial transcriptomics hackathon.

Dataset: [LINK]

Approach wikipedia: 
1. Delaunay triangulation: https://gwlucastrig.github.io/TinfourDocs/DelaunayIntro/index.html
2. Hidden Markov: https://en.wikipedia.org/wiki/Hidden_Markov_model#:~:text=Hidden%20Markov%20models%20are%20generative,emission%20probabilities)%2C%20is%20modeled.

The project breakdown:
1) Wrangling data formatting - OUTPUT: (x,y coordinates of every cell + cell annotations)
2) Delaunay's Triangulation - https://rdrr.io/cran/interp/man/tri.mesh.html;INPUT: x,y coordinates of every cell OUTPUT: Frequency matrix for neighbor probability (cell-frequency of having a cell-type as its neighbor)
3) Permutation testing - INPUT: randomly assign the cell labels OUTPUT: null distribution frequency matrix
(NMF or HMM)?
