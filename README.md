# FilterAirComp
This repository contains the simulation code for my 2023 GLOBECOM paper, to be published in the workshop on "Next Generation Multiple Access (NGMA) for Future Wireless Communications". I also uploaded a preprint version of this paper to Arxiv, you can find that version here: "https://arxiv.org/search/eess?searchtype=author&query=Hellstr%C3%B6m,+H". The title of the paper is "Optimal Receive Filter Design for Misaligned Over-the-Air Computation".

In the code, you will see that the main Monte-Carlo loop has two versions "montecarlo_efficient.m" and "montecarlo_pedagogical.m". The former requires significantly fewer flops to run but the latter follows the structure of the paper more clearly. Both of these versions should yield identical results if initiated with the same seed.

There is a main file which contains the same setting as I used to generate the results in the paper.
