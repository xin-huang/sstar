/*  Link in this file for random number generation with rand()
     and seeding from the clock  */
/*  Thanks to Alan Rogers for suggestion of using pid.   17 Nov 2018 */


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

	double
ran1()
{
	int rand();
	return( rand()/(RAND_MAX+1.0)  );
}


	void seedit( const char *flag)
{
	FILE *fopen(), *pfseed;
	unsigned int seed2, tempseed ;

  if( flag[0] == 's' ) {
      time_t      currtime = time(NULL);
      unsigned long pid = (unsigned long) getpid();
      tempseed = (unsigned int)currtime^pid;
	srand( seed2 = tempseed ) ;
        printf("\n%d\n", seed2 );    
	}

}

	int
commandlineseed( char **seeds)
{
	unsigned int seed2 ;
    void srand(unsigned int seed);

	seed2 = atoi( seeds[0] );

	printf("\n%d\n", seed2 );    

	srand(seed2) ; 
	return(1);
}

