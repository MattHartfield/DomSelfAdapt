/* DomSelfAdaptFWD.c

Forward-in-time simulation of an adaptive alleles with dominance and selfing, with linked neutral fragment
(A check of the coalescent simulation)

This is the BATCH version - to be used on cluster machine to produce many simulation replicates at once

Reps program: Reads in pre-calculated tables and introduces neutral mutations based on them

< Add further preamble here once the program is near release - e.g. runtime instructions, etc. >

Simulation uses routines found with the GNU Scientific Library (GSL)
(http://www.gnu.org/software/gsl/)
Since GSL is distributed under the GNU General Public License 
(http://www.gnu.org/copyleft/gpl.html), you must download it 
separately from this file.

*/

/* Preprocessor statements */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stddef.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <sys/stat.h>
#include <sys/types.h>

#define INITTS 500

/* Function prototypes */
void Wait();
void fitness(unsigned int N, double h, double s, unsigned int *inpop, double *outfit, double *outcf, double *fitsum);
double sumT_UI(unsigned int *Tin, unsigned int nrow);
unsigned int parchoose(unsigned int N, double *culfit, double *fitsum, const gsl_rng *r);
unsigned int bitswitch(unsigned int x);
void generation(unsigned int N, double self, double R, unsigned int *nums, unsigned int *selin, unsigned int *selout, unsigned int **neutin, unsigned int **neutout, double *culfit, double *fitsum, unsigned int npoly, const gsl_rng *r);
void addpoly(unsigned int N, unsigned int **neutin, double *polyloc, unsigned int *npoly, double theta, const gsl_rng *r);
void reassign(unsigned int **neutin, unsigned int **neutout, unsigned int *selin, unsigned int *selout, unsigned int npoly, unsigned int N);
void reassign2(unsigned int **neutin, unsigned int **neutout, double *posin, double *posout, unsigned int npoly, unsigned int N);
void polyprint(unsigned int **neutin, double *posin, unsigned int npoly, unsigned int N);
void ptrim(unsigned int size, unsigned int **neutin, unsigned int **neutout, double *posin, double *posout, unsigned int npoly, unsigned int *npolyT, double *avpi);
void mutsamp(unsigned int N, unsigned int samps, double *polypos, unsigned int **neutin, unsigned int *nums2, unsigned int npoly, unsigned int ei, unsigned int rep, const gsl_rng *r);
void popread(unsigned int **neutin, double *posin, unsigned int npoly, unsigned int N, unsigned int ei, unsigned int suffix);

void Wait(){
	printf("Press Enter to Continue");
	while( getchar() != '\n' );
	printf("\n");	
}

/* Fitness calculation routine */
void fitness(unsigned int N, double h, double s, unsigned int *inpop, double *outfit, double *outcf, double *fitsum){
	unsigned int i;				/* Pop counter */
	unsigned int ac = 0;		/* Allele copies in indv */
	
	/* Calculating new fitnesses */
    *fitsum = 0;
	for(i = 0; i < N; i++){
		ac = 0;
		ac = (*(inpop + 2*i)) + (*(inpop + 2*i + 1));
		if(ac == 0){
			*(outfit + i) = 1;
		}else if(ac == 1){
			*(outfit + i) = 1 + h*s;
		}else if(ac == 2){
			*(outfit + i) = 1 + s;
		}
		*fitsum += *(outfit + i);
		*(outcf + i) = *fitsum;
    }
}	/* End of fitness routine */

/* Summing entire table (unsigned int) */
double sumT_UI(unsigned int *Tin, unsigned int nrow){
	unsigned int i;
	double res = 0;
	for(i = 0; i < nrow; i++){
		res += (*(Tin + i));
	}
	return res;
}

/* Parent choosing */
unsigned int parchoose(unsigned int N, double *culfit, double *fitsum, const gsl_rng *r){
		
	double a;					/* Random number */
	unsigned int done = 0;		/* Signal to determine whether selection complete */
	unsigned int cases = 0;		/* Determines whether one should add up or decrease */
	unsigned int choose1 = 0;	/* Chromosome to be selected */
	unsigned int firstgo = 0;	/* First try at choosing progeny */
	
	done = 0;
	cases = 0;
	choose1 = 0;
	a = gsl_ran_flat(r,0.0,1.0);
	firstgo = (int)N*a;
	
	/* Choosing case */
	if (((*(culfit + firstgo)/(*(fitsum))) >= a) && firstgo == 0){
		cases = 3;
	} else if (((*(culfit + firstgo)/(*(fitsum))) >= a) && firstgo != 0){
		cases = 2;
	} else {
		cases = 1;
	}
	
	switch(cases)
	{
		case 1:
			while(done != 1){
				firstgo++;
				if((*(culfit + firstgo)/(*(fitsum))) >= a){
					choose1 = firstgo;
					done = 1;
				}
			}
			break;
		case 2:
			while(done != 1){
				firstgo--;
				if((*(culfit + firstgo)/(*(fitsum))) < a || firstgo == 0){
					choose1 = (firstgo + 1);
					done = 1;
				}
			}
		break;
		case 3:
			choose1 = 0;
			break;
	}
	
	return choose1;
}

/* Bit-switching routine */
unsigned int bitswitch(unsigned int x){
	unsigned int ret;
	switch(x)
	{
		case 0:
			ret = 1;
			break;
		case 1:
			ret = 0;
			break;		
		default:	/* If none of these cases chosen, exit with error message */
			fprintf(stderr,"Error in bitswitch: input is not 0 or 1.\n");
			exit(1);
            break;
	}
	return ret;
}

/* Reproduction routine */
void generation(unsigned int N, double self, double R, unsigned int *nums, unsigned int *selin, unsigned int *selout, unsigned int **neutin, unsigned int **neutout, double *culfit, double *fitsum, unsigned int npoly, const gsl_rng *r){
	unsigned int i, j, k;		/* Pop counter, neutral marker counter, rec counter */
	unsigned int isself = 0;	/* Is this a selfing reproduction? */
	unsigned int choose1 = 0;	/* Chromosome to be selected */
	unsigned int choose2 = 0;	/* Chromosome to be selected (2nd with outcrossing) */
	unsigned int wc1 = 0;		/* Which chrom (parent 1) */
	unsigned int wc2 = 0;		/* Which chrom (parent 2) */
	unsigned int index1 = 0;	/* Index of chosen sample 1 */
	unsigned int index2 = 0;	/* Index of chosen sample 2 */
	unsigned int nself = 0;		/* Number of selfing events */
	unsigned int nrec = 0;		/* Number of rec events */
	unsigned int selfix = 0;	/* Selfing index */
	unsigned int recix = 0;		/* Recombination index */
	unsigned int repix = 0;		/* Reproduction index */	
	
	/* First deciding number of selfing, recombination events */
	nself = gsl_ran_binomial(r, self, N);
	nrec = gsl_ran_binomial(r, R, (2*N-nself));
	
	/* Creating vector of events */
	unsigned int *selfev = calloc(nself,sizeof(unsigned int));				/* Which events are selfing */
	unsigned int *recev = calloc(nrec,sizeof(unsigned int));				/* Which events are recombinants */
	unsigned int *maxrec = calloc((2*N-nself),sizeof(unsigned int));		/* Max number of rec events (to draw from) */
	for(k = 0; k < (2*N-nself); k++){
		*(maxrec + k) = k;
	}
	
	/* Choosing events */
	gsl_ran_choose(r, selfev, nself, nums, N, sizeof(unsigned int));
	gsl_ran_choose(r, recev, nrec, maxrec, (2*N-nself), sizeof(unsigned int));
		
	for(i = 0; i < N; i++){		/* Regenerating population */
    	
    	/* Choosing first parent, and relevant chromosomes */
		choose1 = parchoose(N, culfit, fitsum, r);
		wc1 = gsl_ran_bernoulli(r,0.5);
		wc2 = gsl_ran_bernoulli(r,0.5);
		index1 = 2*choose1 + wc1;
		
		/* Copying over selected sites, first checking if selfing occurs */
		isself = 0;
		if(nself != 0){
			if(*(selfev + selfix) == i){
				isself = 1;
				selfix++;
			}
		}
		if(isself == 1){
			index2 = 2*choose1 + wc2;
		}else if(isself == 0){
			choose2 = parchoose(N, culfit, fitsum, r);
			index2 = 2*choose2 + wc2;
		}
		(*(selout + 2*i)) = (*(selin + index1));
		(*(selout + 2*i + 1)) = (*(selin + index2));
		
		/* Now copying over neutral fragment, accounting for recombination */
		if(nrec != 0){
			/*
			printf("First entry is...%d\n",*(recev + recix));
			printf("Repix is %d\n",repix);
			*/
			if(*(recev + recix) == repix){
				recix++;
				wc1 = bitswitch(wc1);
				index1 = 2*choose1 + wc1;
				if(isself == 1){
					wc2 = bitswitch(wc2);
					index2 = 2*choose1 + wc2;
				}
			}
		}
		repix++;
		
		if(nrec != 0 && isself == 0){
			if(*(recev + recix) == repix){
				recix++;
				wc2 = bitswitch(wc2);
				index2 = 2*choose2 + wc2;
			}
			repix++;
		}
		
		for(j = 0; j < npoly; j++){
			*((*(neutout + 2*i)) + j) = *((*(neutin + index1)) + j);
			*((*(neutout + 2*i + 1)) + j) = *((*(neutin + index2)) + j);
		}
		
		/*		
		if(isrec == 1){
			indext = index2;
			index2 = index1;
			index1 = indext;
		}
		*/
				
	}

	free(maxrec);
	free(recev);
	
	free(selfev);
	
}	/* End of reproduction routine */

/* Adding neutral polymorphism routine */
void addpoly(unsigned int N, unsigned int **neutin, double *polyloc, unsigned int *npoly, double theta, const gsl_rng *r){
	unsigned int i, j;
	int k;
	unsigned int newpoly = 0;			/* Number of new polymorphisms added */
	unsigned int pchr = 0;				/* Chromosome of appearance */
	unsigned int pfound = 0;			/* Flag to denote if right index found */
	double ploc = 0;					/* Location of polymorphism */
/*	FILE *ofp_poly;				 Pointer for polymorphisms output */	
	
	newpoly = gsl_ran_poisson(r,(0.5*theta));
/*	printf("Add %d new polymorphisms\n",newpoly);*/
	
	for(j = 0; j < newpoly; j++){
		ploc = gsl_ran_flat(r,0.0,1.0);
/*		printf("ploc is %lf\n",ploc);		*/
		pchr = gsl_rng_uniform_int(r,2*N);
/*		printf("pchr is %d\n",pchr);		*/
		
		/* Inserting new polymorphism */
		if(*(npoly) == 0){
/*			printf("zero entry here\n");	*/
			*(polyloc + 0) = ploc;
				for(i = 0; i < 2*N; i++){
					*((*(neutin + i)) + 0) = 0;
				}
			*((*(neutin + pchr)) + 0) = 1;
		}else if(*(npoly) > 0){
/*			printf("non-zero entry here\n");		*/
			pfound = 0;
			for(k = (*(npoly) - 1); k >= 0; k--){
/*				printf("k is %d\n",k);*/
				if(*(polyloc + k) > ploc){	
/*					printf("is here 1\n");*/
					*(polyloc + k + 1) = *(polyloc + k);
					for(i = 0; i < 2*N; i++){
						*((*(neutin + i)) + k + 1) = *((*(neutin + i)) + k);
					}
				}else if(*(polyloc + k) <= ploc){
/*					printf("is here 2\n");				*/
					pfound = 1;
					*(polyloc + k + 1) = ploc;
					for(i = 0; i < 2*N; i++){
						*((*(neutin + i)) + k + 1) = 0;
					}
					*((*(neutin + pchr)) + k + 1) = 1;
					break;
				}			
			}

			if(pfound == 0){		/* In this case polymorphism is inserted at the start */
				*(polyloc + 0) = ploc;
				for(i = 0; i < 2*N; i++){
					*((*(neutin + i)) + 0) = 0;
				}
				*((*(neutin + pchr)) + 0) = 1;
			}
			
		}
		
		*(npoly) += 1;
/*		printf("n poly is %d\n",*(npoly));			*/
		if(*(npoly) == INITTS){
			fprintf(stderr,"Too many neutral polymorphisms.\n");
			exit(1);
		}
		
	}
}

/* Reassigning new population routine */
void reassign(unsigned int **neutin, unsigned int **neutout, unsigned int *selin, unsigned int *selout, unsigned int npoly, unsigned int N){
	unsigned int i, j;		/* pop counter, poly counter */
	
	for(i = 0; i < 2*N; i++){
		*(selout + i) = *(selin + i);
		for(j = 0; j < npoly; j++){
			*((*(neutout + i)) + j) = *((*(neutin + i)) + j);
		}
	}
	
}	/* End of reassignment routine */

/* Reassigning new population routine (after trimming) */
void reassign2(unsigned int **neutin, unsigned int **neutout, double *posin, double *posout, unsigned int npoly, unsigned int N){
	unsigned int i, j;		/* pop counter, poly counter */
	
	for(i = 0; i < 2*N; i++){
		for(j = 0; j < npoly; j++){
			*((*(neutout + i)) + j) = *((*(neutin + i)) + j);
			if(i == 0){
				*(posout + j) = *(posin + j);
			}
		}
	}
	
}	/* End of reassignment routine */

/* Printing polymorphism routines */
void polyprint(unsigned int **neutin, double *posin, unsigned int npoly, unsigned int N){
	unsigned int i, j;		/* pop counter, poly counter */
	
	for(j = 0; j < npoly; j++){
		printf("%lf ",*(posin + j));
		for(i = 0; i < 2*N; i++){
			printf("%d ",*((*(neutin + i)) + j));
		}
		printf("\n");
	}
	Wait();
	
}	/* End of reassignment routine */

/* Cutting out non-polymorphic sites */
void ptrim(unsigned int size, unsigned int **neutin, unsigned int **neutout, double *posin, double *posout, unsigned int npoly, unsigned int *npolyT, double *avpi){
	unsigned int i, j;
	unsigned int newp = 0;		/* New number of polymorphisms */
	unsigned int count = 0;		/* Polymorphism count */
	double cumpi = 0;			/* Cumulative pi (diversity) */

/*	printf("Trimming routine activated\n");	*/
	for(j = 0; j < npoly; j++){
		count = 0;
		*(posout + newp) = *(posin + j);
		for(i = 0; i < size; i++){
			*((*(neutout + i)) + newp) = *((*(neutin + i)) + j);
			count += *((*(neutin + i)) + j);
		}
/*		printf("For poly %d, count is %d\n",j,count);*/
		if((count > 0) && (count < size)){
/*			printf("Keeping poly %d\n",j);	*/
			newp++;
			cumpi += ((count*(size-count))/(1.0*size*size));
/*			printf("Count is %d, freq is %lf, cumulative value is %lf\n",count,count/(1.0*size),cumpi);	*/
		}
	}
	
	*npolyT = newp;
	*avpi = (cumpi/(1.0*newp));
}

/* Producing mutational samples */
void mutsamp(unsigned int N, unsigned int samps, double *polypos, unsigned int **neutin, unsigned int *nums2, unsigned int npoly, unsigned int ei, unsigned int rep, const gsl_rng *r){

	unsigned int a, j, x;
	unsigned int currsamp;
	unsigned int npolyT;
	double avpi;
	char Mout[32];				 /* String to hold filename in (Mutations) */
	FILE *ofp_mut = NULL;		 /* Pointer for data output */
	
	for(x = 0; x < rep; x++){
	
		double *polyposF = calloc(npoly,sizeof(double));						/* Position of neutral mutations (final for polymorphism sampling) */
		unsigned int **nsamp = calloc(samps,sizeof(unsigned int *));			/* Table of neutral markers per individual (final for sampling) */
		unsigned int **nsampT = calloc(samps,sizeof(unsigned int *));			/* Table of neutral markers per individual (final for sampling) AFTER TRIMMING */
		for(a = 0; a < samps; a++){
			nsamp[a] = calloc(npoly,sizeof(unsigned int));	
			nsampT[a] = calloc(npoly,sizeof(unsigned int));				
		}
		
		npolyT = npoly;
	
		/* Choosing which haplotypes to sample */
		unsigned int *thesamp = calloc(samps,sizeof(unsigned int));
		gsl_ran_choose(r, thesamp, samps, nums2, 2*N, sizeof(unsigned int));
	
		for(a = 0; a < samps; a++){
			currsamp = *(thesamp + a);
			for(j = 0; j < npoly; j++){
				*((*(nsamp + a)) + j) = *((*(neutin + currsamp)) + j);
			}
		}
		
		/* Check, then include code to ONLY print out polymorphic sites in table? */
		ptrim(samps, nsamp, nsampT, polypos, polyposF, npoly, &npolyT, &avpi);
	
		/* Printing out sample table */
		sprintf(Mout,"Mutations/Muts_%d.dat",ei*rep + x);
		ofp_mut = fopen(Mout,"w");
		for(j = 0; j < npolyT; j++){
			fprintf(ofp_mut,"%lf ",*(polyposF + j));
			for(a = 0; a < samps; a++){
				fprintf(ofp_mut,"%d ",*((*(nsampT + a)) + j));
			}
			fprintf(ofp_mut,"\n");
		}
		fclose(ofp_mut);
	
		for(a = 0; a < samps; a++){
			free(nsampT[a]);		
			free(nsamp[a]);
		}
		free(nsampT);		
		free(nsamp);
		free(thesamp);
		free(polyposF);
	
	}
}

/* Reading population state from file */
void popread(unsigned int **neutin, double *posin, unsigned int npoly, unsigned int N, unsigned int ei, unsigned int suffix){
	unsigned int i, j;		/* pop counter, poly counter */
	char Pin[32];				 /* String to hold filename in (Mutations) */
	FILE *ifp_poly = NULL;		 /* Pointer for data output */
	
	/* Printing out polymorphism table */
	sprintf(Pin,"Pop_%d.dat",ei + suffix);
	ifp_poly = fopen(Pin,"r");
	for(j = 0; j < npoly; j++){
		/* Is this the right way of reading in the first column of file? */
		for(i = 0; i < (2*N + 1); i++){
			if(i == 0){
				fscanf(ifp_poly,"%lf",&(posin[j]));
			}else if(i > 0){
				fscanf(ifp_poly,"%d",&(neutin[(i-1)][j]));
			}
		}
		fprintf(ifp_poly,"\n");
	}
	fclose(ifp_poly);
		
}	/* End of population reading routine */

/* Main program */
int main(int argc, char *argv[]){

	/* Declare variables here */
	unsigned int a, i, x;		/* Counters */
	unsigned int N = 0;			/* Population Size */
	unsigned int initp = 0;		/* Initial location of mutant (individual) */
	unsigned int done = 0;		/* Is simulation done or not? */
	unsigned int dcopies = 0;	/* Number of copies of derived allele */
	unsigned int afix = 0;		/* Has allele fixed? */
	unsigned int npoly = 0;		/* Number of neutral alleles */
	unsigned int npolyT = 0;	/* Number of neutral alleles (after cutting out fixed sites) */	
	unsigned int nowsel = 0;	/* If derived allele selected for? */
	unsigned int samps = 0;		/* Number of samples to take */
	unsigned int suffix = 0;	/* Number of times to subsample from population */
	unsigned int base = 0;		/* Number of baseline forward-in-time simulations */
	unsigned int rps = 0;		/* Samples taken per baseline simulation */
/*	unsigned int n = 0;			sprintf counter */
	unsigned int ttot = 0;		/* Total time of simulation */
	double s = 0;				/* Strength of selection */
	double h = 0;				/* Dominance level */
	double self = 0;			/* Selfing rate */
	double Rin = 0;				/* Recombination rate (input) */
	double R = 0;				/* Recombination rate (at a time) */	
	double x0 = 0;				/* Initial allele frequency */
	double theta = 0;			/* Rate of neutral mutation */
	double fitsum = 0;			/* Summed Fitness */
	double scurr = 0;			/* Current fitness of mutant (to account for initial neutrality) */
	double avpi = 0;
	char Pout[32];				/* String to hold filename in (Mutations) */	
	char Sout[32];				/* String to hold filename in (Seed) */		
	FILE *ofp_sd;				/* Pointer for seed output */
	FILE *ofp_poly;				/* Pointer for polymorphisms output */
	
	/* GSL random number definitions */
	const gsl_rng_type * T;
	gsl_rng * r;
	
	/* Reading in data from command line */
	if(argc != 12){
		fprintf(stderr,"Not enough inputs (need: N s h self R x0 theta basereps samples suffix npoly).\n");
		exit(1);
	}
	
	N = atoi(argv[1]);
	if(N <= 0){
		fprintf(stderr,"Total Population size N is zero or negative, not allowed.\n");
		exit(1);
	}
	
	s = strtod(argv[2],NULL);
	if(s < 0){
		fprintf(stderr,"Allele strength s is negative, not allowed.\n");
		exit(1);
	}
	
	h = strtod(argv[3],NULL);
	if(h < 0 || h > 1){
		fprintf(stderr,"Dominance value must lie between 0 and 1.\n");
		exit(1);
	}
	
	self = strtod(argv[4],NULL);
	if(self < 0 || self > 1){
		fprintf(stderr,"Selfing rate must lie between 0 and 1 (inclusive).\n");
		exit(1);
	}
		
	Rin = strtod(argv[5],NULL);
	if(Rin < 0){
		fprintf(stderr,"Recombination rate must be positive or zero.\n");
		exit(1);
	}
	Rin /= 2.0*N;
	
	x0 = strtod(argv[6],NULL);
	if(x0 < (1.0/(2.0*N)) || x0 > 1.0-(1.0/(2.0*N)) ){
		fprintf(stderr,"Initial mutant frequency must lie between 1/2N and 1-1/2N.\n");
		exit(1);
	}
	
	theta = strtod(argv[7],NULL);
	if(theta <= 0){
		fprintf(stderr,"Mutation rate must be a positive value.\n");
		exit(1);
	}
	
	base = atoi(argv[8]);
	if(base <= 0){
		fprintf(stderr,"Number of baseline simulations must be greater than zero.\n");
		exit(1);
	}
	
	rps = atoi(argv[9]);
	if(rps <= 0){
		fprintf(stderr,"Samples must be > 0.\n");
		exit(1);
	}
	
	suffix = atoi(argv[10]);
	if(argv[10] < 0){
		fprintf(stderr,"File index must be greater than or equal to zero.\n");
		exit(1);
	}
	
	npoly = atoi(argv[11]);
	
	/* create a generator chosen by the 
    environment variable GSL_RNG_TYPE */
    
    mkdir("SeedsRep/", 0777); 
	gsl_rng_env_setup();
	if (!getenv("GSL_RNG_SEED")) gsl_rng_default_seed = time(0);
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	sprintf(Sout,"SeedsRep/Seed_%d.dat",suffix);
	ofp_sd = fopen(Sout,"w");
	fprintf(ofp_sd,"%lu\n",gsl_rng_default_seed);
	fclose(ofp_sd);

	mkdir("Polymorphisms/", 0777);		
	mkdir("Mutations/", 0777);
	
	/* Vector of 0 to (N-1) for sampling random numbers */
	unsigned int *nums = calloc(N,sizeof(unsigned int));
	unsigned int *nums2 = calloc(2*N,sizeof(unsigned int));
	for(a = 0; a < 2*N; a++){
		*(nums2 + a) = a;
		if(a < N){
			*(nums + a) = a;
		}
	}
	
	/* Executing simulation */
	for(i = 0; i < base; i++){

/*		printf("Starting run %d\n",i);*/
	
		afix = 0;
		R = Rin;
		while(afix != 1){
	
			/* Setting up selection, neutral tables per individual */
			double *fit = calloc(N,sizeof(double));									/* Fitness of each individual */
			double *cumfit = calloc(N,sizeof(double));								/* Cumulative fitness of each individual */
			double *polypos = calloc(INITTS,sizeof(double));						/* Position of neutral mutations */
			double *polyposF = calloc(INITTS,sizeof(double));						/* Position of neutral mutations (final for polymorphism sampling) */
			unsigned int *selindvP = calloc(2*N,sizeof(unsigned int));				/* Table of derived allele states per individual (parents) */
			unsigned int *selindvO = calloc(2*N,sizeof(unsigned int));				/* Table of derived allele states per individual (offspring) */
			unsigned int **neutindvP = calloc(2*N,sizeof(unsigned int *));			/* Table of neutral markers per individual (parents) */
			unsigned int **neutindvO = calloc(2*N,sizeof(unsigned int *));			/* Table of neutral markers per individual (offspring) */
			unsigned int **neutindvF = calloc(2*N,sizeof(unsigned int *));			/* Table of neutral markers per individual (final for sampling) */
			for(a = 0; a < 2*N; a++){
				neutindvP[a] = calloc(INITTS,sizeof(unsigned int));	
				neutindvO[a] = calloc(INITTS,sizeof(unsigned int));
				neutindvF[a] = calloc(INITTS,sizeof(unsigned int));			
			}
			ttot = 0;

			/* Reassigning ancestral polymorphism state */
			popread(neutindvP, polypos, npoly, N, i, suffix);
/*
			polyprint(neutindvP, polypos, npoly, N);
			exit(1);
*/
			/* Assigning initial mutant, calculating fitness */
			initp = gsl_rng_uniform_int(r,2*N);
			*(selindvP + initp) = 1;
		
			if(x0 == 1.0/(2.0*N)){
				scurr = s;
				nowsel = 1;
			}else if(x0 > 1.0/(2.0*N)){
				scurr = 0;
				nowsel = 0;
			}
	/*		printf("scurr is %lf\n",scurr);	*/
			fitness(N, h, scurr, selindvP, fit, cumfit, &fitsum);
	
			done = 0;
			while(done != 1){
			
				/* Generating new population */
				generation(N,self,R,nums,selindvP,selindvO,neutindvP,neutindvO,cumfit,&fitsum,npoly,r);
		
				/* Neutral Mutation phase */
				addpoly(N, neutindvO, polypos, &npoly, theta, r);
			
				/* Reassigning matrices */
				reassign(neutindvO, neutindvP, selindvO, selindvP, npoly, N);
				
				ttot++;
				if(ttot%1 == 0){
					ptrim(2*N, neutindvP, neutindvF, polypos, polyposF, npoly, &npolyT,&avpi);
					npoly = npolyT;
					reassign2(neutindvF, neutindvP, polyposF, polypos, npoly, N);
				}
			
				/* Checking derived allele copies and whether it has fixed or lost */
				dcopies = sumT_UI(selindvP, 2.0*N);
/*				printf("Dcopies are %d\n",dcopies);				*/
/*				printf("copies are %d, npoly are %d\n",dcopies,npoly);*/
				if(dcopies == 0){
					done = 1;
				}
				if(dcopies == 2.0*N){
					done = 1;
					afix = 1;
				}
			
				/* Deciding if derived allele becomes selected for, then updating fitness */
				if(((dcopies/(2.0*N)) >= x0) && nowsel == 0){
					nowsel = 1;
					scurr = s;
	/*				printf("scurr is %lf\n",scurr);	*/
				}
				fitness(N, h, scurr, selindvP, fit, cumfit, &fitsum);
			}
		
			/* Routine to print out polymorphisms if allele fixed*/
			if(afix == 1){
			
				printf("Selected allele fixed\n");
			
				/* First, remove non-polymorphic sites */
				ptrim(2*N, neutindvP, neutindvF, polypos, polyposF, npoly, &npolyT,&avpi);
				npoly = npolyT;
		
				/* Printing out polymorphism table */
				sprintf(Pout,"Polymorphisms/Poly_%d.dat",i + suffix);
				ofp_poly = fopen(Pout,"w");
				for(a = 0; a < npoly; a++){
					fprintf(ofp_poly,"%lf ",*(polyposF + a));
					for(x = 0; x < 2*N; x++){
						fprintf(ofp_poly,"%d ",*((*(neutindvF + x)) + a));
					}
					fprintf(ofp_poly,"\n");
				}
				fclose(ofp_poly);
		
				/* Then, sample from this table to produce pseudo-coalescent results */
				samps = 10;
				mutsamp(N, samps, polyposF, neutindvF, nums2, npoly, i + suffix, rps, r);
			
			}
		
			/* End of run; freeing memory */
			for(a = 0; a < 2*N; a++){
				free(neutindvF[a]);		
				free(neutindvO[a]);
				free(neutindvP[a]);
			}	
			free(neutindvF);
			free(neutindvO);
			free(selindvO);
			free(neutindvP);
			free(selindvP);
			free(polyposF);
			free(polypos);
			free(cumfit);		
			free(fit);
		}
	
	}

	free(nums2);	
	free(nums);
	gsl_rng_free(r);	
	return 0;
	
}	/* End of main program */

/* End of File */