README FOR DOM SELF ADAPT

Forward-in-time simulations for determining neutral diversity around a sweeping beneficial allele, while considering arbitrary levels of self-fertilisation and dominance. As used in the study "Selective sweeps under dominance and self-fertilisation".

Simulation uses routines found with the GNU Scientific Library (GSL) (http://www.gnu.org/software/gsl/) Since GSL is distributed under the GNU General Public License (http://www.gnu.org/copyleft/gpl.html), you must download it separately from this file.

Programs can be compiled with e.g. GCC using a command like: gcc DomSelfAdaptFWD_BIRec -lm -lgsl -lgslcblas -I/usr/local/include -L/usr/local/lib DomSelfAdaptFWD_BIRec.c (similar for DomSelfAdaptFWD_RepsRec.c).

'DomSelfAdaptFWD_BIRec.c' needs to be run first in order to generate neutral diversity (i.e. in the absence of selection). The program is executed as:

./DomSelfAdaptFWD_BIRec N self R 4Nu basereps suffix

Where:
- N is the population size
- self is the rate of self-fertilisation (related to inbreeding rate F = self/(2-self)
- R is the population-level recombination rate 2Nr across the entire sample
- 4Nu is the population-level mutation rate
- base is the number of simulations to run (i.e. number of repetitions). Must be set to one.
- suffix is the baseline numbering for indexing simulation outputs

The outputs are files "Pop_X.dat" where X is the rep+suffix number. These list the location of each polymorphism and it's state for all 2N haplotypes (0 = ancestral; 1 = derived).

Next one runs DomSelfAdaptFWD_RepsRec to simulate a sweep from these burn-in populations. The program is executed as:

./DomSelfAdaptFWD_RepsRec N s h self R x0 xfix 4Nu basereps reps-per-sim samples suffix npoly

Where:
- N is the population size
- s is the homozygote selective advantage
- h is the dominance coefficient (has to lie between 0 and 1)
- self is the rate of self-fertilisation (related to inbreeding rate F = self/(2-self))
- R is the population-level recombination rate 2Nr across the entire sample
- x0 is the frequency at which the derived allele becomes advantageous
- xfix is the fixation frequency of the sweep, at which point haplotypes are sampled
- 4Nu is the population-level mutation rate
- basereps is the number of simulations to run (i.e. number of repetitions). Must be set to one.
- reps-per-sim is the number of times haplotypes are sampled from the fixed state
- samples is the number of haplotypes to sample from the final population
- suffix is the index of the first population file to read in (number matches 'X' from 'Pop_X.dat')
- npoly is the number of polymorphisms of the inputted Pop file.

The outputs are placed in two folders: 'Polymorphisms', which is a file of polymorphisms for the entire population; and 'Mutations', which are files of sampled mutations. Note that the state of the selected allele (at location 0) is also provided.

Note that although it is possible to define more than one repetition of each simulation, these programs are designed to be run on a cluster machine, where each node only executes a single neutral population (for DomSelfAdaptFWD_BIRec) or selective sweep (for DomSelfAdaptFWD_RepsRec). Hence it is recommended to only run one rep per simulation.