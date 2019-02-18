self.sites:  
(actual undisturbed lattice sites)
[[ 0.          0.          0.        ]
 [ 2.71055907  1.56494201  4.42632442]
...
 [29.81614973 17.21436208  8.85264884]
 [32.5267088  18.77930408 13.27897326]]

self.state: 
['S' 'S' 'S' 'M' 'M' 'M' 'A' 'A' 'A' 'A' ....'V' 'A' 'A' ]

ostr:
ostr = "".join(self.state)  = e.g. AAAAAAAMAAASAAAAMAAAAAAAAAAASASAAMAAAAVA...

ivac (constant): if 64 atoms sits & 1vac -> ivac = [63]
ivac (constant): if 64 atoms sits & 2vac -> ivac = [63 ,62]
ivac (constant): is always the index lattice site associated the vacancy (e.g. 63 when 64 sites)

svac = self.idx[ivac] # lattice site associated with this vacancy
svac: index of the current vacancy position 

self.neigh 

ridx: [ 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47
 48 49 63 51 52 53 54 55 56 57 58 59 60 61 62 50]
ridx: # reverse lookup index [i.e. ridx[i] gives the index of the atoms/vacancy at site i] (where the 'site' is the original site when the simulation started) 


idx: connects the initial positions of the atoms 
self.idx[ivac] : is sorted in such a way that last position gives alwasy the index 
to the vacancy.

svac = self.idx[ivac]  # ivac is consatnt = 63
      (self.idx is changed with every swap to connect to original lattice sites)
sneigh = self.neigh[svac] -> loop over all neighbors
		is every neighbor of vacancy self.neigh[svac]


in self.ridx tauschen immer svac <--> sneigh  (das ausgesuchte sneigh)

in self.idx  tauschen immer ivac <--> ineigh  

ineigh = self.ridx[sneigh] # gets index of atom associated with the neighboring site

%cat vacancy_pos_unwrapped_ka_0.txt
[17.21239301  9.93757974  7.02693002] svac 63 sneigh 50 ivac 63 ineigh 50
[11.47492868  1.65626329  4.68462002] svac 50 sneigh 54 ivac 63 ineigh 54
[12.90929476  4.14065823  4.68462002] svac 54 sneigh 5  ivac 63 ineigh 5
[2.86873217 3.31252658 2.34231001]    svac 5  sneigh 6  ivac 63 ineigh 6
[4.30309825 4.14065823 4.68462002]    svac 6  sneigh 2  ivac 63 ineigh 2

difference between idx und ridx:
step 0:
atom 50 geht auf vacancy position 63 -> site 50 is now vacant
svac = 50
 idx [ 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47
 48 49 63 51 52 53 54 55 56 57 58 59 60 61 62 50]
ridx [ 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47
 48 49 63 51 52 53 54 55 56 57 58 59 60 61 62 50]
 idx[63] = atom = 50 
ridx[63] = atom = 50
 idx[50] = vac  = 63
ridx[50] = vac  = 63
self.idx[ivac=63] (=50) -> gives always the index of where the vacancy is currently
self.ridx[svac] is always ivac == 63
self.idx[ivac] is always svac == 50
self.ridx[ivac=63] (=50) -> gives always the index of where the vacancy is currently

step 1:
atom 54 geht auf vacancy position 50: -> site 54 is now vacant
svac = 54
 idx [ 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47
 48 49 63 51 52 53 50 55 56 57 58 59 60 61 62 54]  -> shows 
ridx [ 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47
 48 49 54 51 52 53 63 55 56 57 58 59 60 61 62 50]
 idx[63] = 54
ridx[63] = 54
 idx[50] = 63
ridx[50] = 63
self.idx[ivac=63] (54) -> gives always the index of where the vacancy is currently
self.ridx[svac] is always ivac == 63
self.idx[ivac] is always svac == 54
self.ridx[ivac=63] (54) -> gives always the index of where the vacancy is currently

self.idx[63]  = 54 -> original site 63 is currently on site 54 
self.ridx[54] = 63 -> current site 54 is originally from site 54 
self.ridx[svac]    -> always 63 (for every step)
self.sites[54]     -> always gives array([24.3950316 ,  7.82471003,  8.85264884])
				   -> for every step
