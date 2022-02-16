Explication of import rules

Current import rules:

first line contains all the interactions, and then all of the levels of the chemicals to start with

next, the chemcial reactions for each species are created. This gives the total number of reaction terms for each species, with a reaction term being a single polynomial expression arising from mass action kinetics that has a prefactor (positive or negative), which is the rate.

We can extend this to larger systems in the following way:

first line is the total number of components
next line is the species which are undergoing phase separation (1s and zeros)
the next line is (n-1)*(n-2)/2 off diagonal components of the interaction matrix
the next line is the amount of each chemical that we begin with
then each following line will be the chemistry affecting species i in the setting described above.

