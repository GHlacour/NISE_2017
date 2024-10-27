# NISE_2017
This is the current development version of a quantum classical package for calculating coherent multidimensional spectra (as FTIR, SFG, 2DIR, 2DES, 2DIRraman, 2DSFG, and F-2DES). 
# General Description
The NISE3.1 code was originally developed by Thomas la Cour Jansen. Please, cite the appropriate references [1–6] when publishing work using this code. The code allows the calculation of the linear absorption,linear dichroism, sum-frequency generation, two-dimensional spectra (IR,UVvis, and SFG), population transfer, exciton diffusion and integrated anisotropy using the full nonadiabatic semi-classical numerical integration of the Schrödinger equation approach [2] and the sparse matrix optimization approach [4]. The code allows treating spectra of diverse systems involving intra- and intermolecular energy transfer [1,3,7–9], non-Gaussian dynamics [10, 11], surfaces [5], and chemical exchange [12]. This manual is not intended as an introduction to two-dimensional spectroscopy. The user is directed to the references including recent reviews [2,13–16,20] and books for more information [17–19]. The code use wavenumbers for frequencies and times are femtoseconds. The transition dipoles and transition polarizabilities may be given in any desired units. This version has MPI and OpenMP implementation for parallel use for all 2D methods [21].
Feedback on the program and the manual are welcome via e-mail: t.l.c.jansen@rug.nl or contribute an issue on the gitHub repository. Change in the code is allowed, but on own risk, and should be reported clearly in publications. Redistribution of the code must happen in accordance with the license.

Hamiltonians for the NISE code can be created with the AIM program [22,23]. An external tutorial is available [24] and a YouTube video demonstration of the installation of the programme [25].

## Official NISE developer team
Computational Spectroscopy group of Thomas la Cour Jansen, University of Groningen\
Group of Ana Cunha, University of Antwerp

# References
1. T. L. C. Jansen and J. Knoester. J. Phys. Chem. B, 110:22910–22916, (2006).
2. T. L. C. Jansen and J. Knoester. Acc. Chem. Res., 42(9):1405–1411, (2009).
3. T. L. C. Jansen, B. M. Auer, M. Yang and J. L. Skinner. J. Chem. Phys., 132:224503, (2010).
4. C. Liang and T. L. C. Jansen. J. Chem. Theory Comput., 8:1706–1713, (2012).
5. C. Liang, M. Louhivuori, S. J. Marrink, T. L. C. Jansen and J. Knoester. J. Phys.
Chem. Lett., 4:448–452, (2013).
6. C. D. N. van Hengel, K. E. van Adrichem and T. L. C. Jansen, J. Chem. Phys. 158, 064106 (2023)
7. D. Cringus, T. L. C. Jansen, M. S. Pshenichnikov and D. A. Wiersma. J. Chem.
Phys., 127:084507, (2007).
8. T. L. C. Jansen and J. Knoester. Biophys. J., 94:1818–1825, (2008).
9. A. G. Dijkstra, T. L. C. Jansen and J. Knoester. J. Phys. Chem. A, 114:7315–7320, (2010).
10. T. L. C. Jansen, D. Cringus and M. S. Pshenichnikov. J. Phys. Chem. A, 113:6260, (2009).
11. S. Roy, M. S. Pshenichnikov and T. L. C. Jansen. J. Phys. Chem. B, 115:5431–5440, (2011).
12. T. L. C. Jansen and J. Knoester. J. Chem. Phys., 127:234502, (2007).
13. P. Hamm, M. H. Lim and R. M. Hochstrasser. J. Phys. Chem. B, 102:6123–6138,
(1998).
14. R. M. Hochstrasser. Chem. Phys., 266(2-3):273–284, (2001).
15. M. Cho. Chem. Rev., 108:1331, (2008).
16. S. Mukamel. Annu. Rev. Phys. Chem., 51:691, (2000).
17. M. Cho. Two-dimensional Optical Spectroscopy. CRC Press, Boca Raton, 2009.
18. S. Mukamel. Principles of Nonlinear Optical Spectroscopy. Oxford University Press,
New York, 1995.
19. P. Hamm and M. T. Zanni. Concepts and Methods of 2D Infrared Spectroscopy. Cambridge University Press, Cambridge, 2011.
20. T. L. C. Jansen J. Chem. Phys. 155 (17): 170901, (2021)
21. A. S. Sardjan, F. P. Westerman, J. P. Ogilvie, and T. L. C. Jansen, J. Phys. Chem. B 124: 9420-9427 (2020).
22. K. E. van Adrichem and T. L. C. Jansen, J. Chem. Theory Comput. 18: 3089–3098 (2022).
23. AIM: [https://github.com/Kimvana/AIM](https://github.com/Kimvana/AIM) 
24. NISE Tutorials: [https://github.com/GHlacour/NISE_Tutorials](https://github.com/GHlacour/NISE_Tutorials)
25. NISE installation video: [YouTube](https://www.youtube.com/watch?v=npvV9UOFmDg) 
