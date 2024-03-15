# TriangulaMarkovCoLocalizer (TMCL)

A-star 10-X genomics spatial transcriptomics hackathon.

Dataset: [LINK]

Approach wikipedia: 
1. Delaunay triangulation: https://gwlucastrig.github.io/TinfourDocs/DelaunayIntro/index.html
2. Hidden Markov: https://en.wikipedia.org/wiki/Hidden_Markov_model#:~:text=Hidden%20Markov%20models%20are%20generative,emission%20probabilities)%2C%20is%20modeled.

The project breakdown:
1) Wrangling data formatting - OUTPUT: (x,y coordinates of every cell + cell annotations)
2) Delaunay's Triangulation - INPUT: x,y coordinates of every cell and cell-type annotations. OUTPUT: Frequency matrix for neighbor probability (cells x cell-type neighbors)
3) n-size grid for tissue map - INPUT: x,y coordinates of every cell and cell-type annotations. OUTPUT: Grid frequency matrix of cell-types (cell-type x n-size grids)
4) Hidden Markov Model - INPUT: emission property matrix - nearest neighbor frequency matrix, transition matrix - cosyne-distance weighted, adjusted for euclidian distance between grids matrix. OUTPUT: Transition state matrix for each cell type.
5) Permutation testing - INPUT: Repeat steps 2-4 with randomly assigned cell labels. OUTPUT: null distribution transition score for each cell-type and p-value for each transition score.
