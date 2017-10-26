/* nfasack.c DNA Fasta file sanity check: principally it will say 
 * if it a multisequence file */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pcre.h>

#define OVECCOUNT 300    /* should be a multiple of 3 */
#define MXFSZ (1<<16)
#define BUF 4 /* tested for a low a value as 2, and is fine */

#if defined DBG || defined DBG2
#define GBUF 4
#else
#define GBUF 128
#endif
#define SSZ 2 /* CG count, first, AT count second, third are the anomalous characters */
#define HISTBUCKETSZ 10
#define HTCOLWIDTH 120
#define NAMESTRSZ 256 // an arbitrary size for a name.

#define CONDREALLOC(x, b, c, a, t); \
	if((x)==((b)-1)) { \
		(b) += (c); \
		(a)=realloc((a), (b)*sizeof(t)); \
		memset((a)+(b)-(c), '\0', (c)*sizeof(t)); \
	}

#define CONDREALLOC2(x, b, c, a1, a2, t); \
	if((x)==((b)-1)) { \
		(b) += (c); \
		(a1)=realloc((a1), (b)*sizeof(t)); \
		memset((a1)+(b)-(c), '\0', (c)*sizeof(t)); \
		(a2)=realloc((a2), (b)*sizeof(t)); \
		memset((a2)+(b)-(c), '\0', (c)*sizeof(t)); \
	}

typedef struct /* onefa */
{
	char *id;
	char *sq;
	unsigned idz, sqz;
} onefa;

typedef struct /* i_s; sequence index and number of symbols */
{
	unsigned idx;
	size_t sylen; /* this is the precise symbol length of the sequence */
	size_t sy[SSZ]; /* used to hold counts of symbols */
	float cgp;
	unsigned ambano[2]; /* number of ambiguous symbols (first), then number of anomalous symbols */
	char *id; // the ID name
	char *sq; // the sequence itself
	unsigned ibf, sbf; // buffers for ID and SQ strings
	unsigned idz, sqz; // actual size  of the ID and SQ strings. Is almost a duplication of sylen, can be removed if memory is a consideration.
} i_s; /* sequence index and number of symbols */

typedef struct /* ou_t unique occurence vector: often uov */
{
	unsigned *ua; /* unique size array */
	unsigned *oa; /* occurence array: will be same size as unique array */
	unsigned ub; /* buffer for the unique numbers */
	unsigned usz; /* size of the two arrays */
} uo_t;

size_t fszfind(FILE *fp)
{
    fseek(fp, 0, SEEK_END);
    size_t fbytsz = (size_t)ftell(fp);
    rewind(fp);
    return fbytsz;
}

void prtuo(uo_t *uov)
{
	unsigned j;
	printf("number of different sequence lengths: %u\n", uov->usz);
	size_t totsz=0;
	for(j=0; j<uov->usz;++j)
		totsz += uov->ua[j]*uov->oa[j];
#ifdef DBG
	printf("vals: "); 
	for(j=0; j<uov->usz;++j)
		printf("%5u ", uov->ua[j]);
	printf("\n"); 
	printf("ocs: "); 
	for(j=0; j<uov->usz;++j)
		printf("%5u ", uov->oa[j]);
	printf("\n"); 
#endif
	printf("Total size of sequences added together=%zu\n", totsz); 
	return;
}

uo_t *uniquelens(i_s *sqi, unsigned numsq)
{
	unsigned char new;
	unsigned i, j;
	unsigned ai=0;
	uo_t *uov=malloc(sizeof(uo_t));
	uov->ub=GBUF;
	uov->ua=calloc(uov->ub, sizeof(unsigned));
	uov->oa=calloc(uov->ub, sizeof(unsigned));
	for(i=0; i<numsq;++i) {
		new=1;
		for(j=0; j<=ai;++j) {
			if(uov->ua[j] == sqi[i].sylen) {
				uov->oa[j]++;
				new=0;
				break;
			}
		}
		if(new) {
			CONDREALLOC2(ai, uov->ub, GBUF, uov->ua, uov->oa, unsigned);
			uov->ua[ai]=sqi[i].sylen;
			uov->oa[ai]++;
			ai++;
		}
	}

	uov->usz=ai;
#if defined DBG || defined DBG2
	prtuo(uov);
#endif

	return uov;
}

void prtfa(onefa *fac)
{
	printf(">");
	printf("%s\n", fac->id);
	printf("%s\n", fac->sq);
}

void prtfaf(onefa *fac, FILE *fp)
{
	fprintf(fp, ">");
	fprintf(fp, "%s\n", fac->id);
	fprintf(fp, "%s\n", fac->sq);
}

void prtsq(i_s *sqisz, int sz)
{
	int i;
	printf("Number of different sequences=%i\n", sz); 
#ifdef DBG
	for(i=0;i<sz;++i) {
		printf("%s\n", sqisz[i].id);
		printf("%s\n", sqisz[i].sq);
	}
#endif
	return;
}

void prtfa2(onefa *fac)
{
	int i;
	printf("SQZ=%d:", fac->sqz);
	for(i=0;i<3;++i) 
		putchar(fac->sq[i]);
	printf("\n"); 
}

void prti_s(i_s *sqisz, int sz, float *mxcg, float *mncg)
{
	int i;
	char *sqgood;
	*mxcg=.0;
	*mncg=1.;

	size_t tsz;
	for(i=0;i<sz;++i) {
		if(sqisz[i].ambano[1] != 0)
			sqgood="AnoSQ";
		else
			sqgood="SQ";
		tsz = sqisz[i].sy[0] + sqisz[i].sy[1];
		sqisz[i].cgp=(float)sqisz[i].sy[0]/tsz;
		if(sqisz[i].cgp>*mxcg)
			*mxcg=sqisz[i].cgp;
		if(sqisz[i].cgp<*mncg)
			*mncg=sqisz[i].cgp;

		printf("| %s#%i=TOT:%zu CG:%.4f ", sqgood, i, sqisz[i].sylen, sqisz[i].cgp);
	}
	printf("|\n"); 
}

i_s *procfastaf(char *fname, unsigned *nsq)
{
	/* general declarations */
	FILE *fin;
	char IGLINE, begline;
	size_t lidx, mxsylen, mnsylen;
	unsigned mxamb, mnamb;
	int i, c, sqidx;
	int gbuf;
	i_s *sqisz=NULL;
	int whatint; // a number reflecting the type of symbol read
	unsigned numsq, numano;
	int ididx0=0;

	// OK open the file
	if(!(fin=fopen(fname, "r")) ) { /*should one check the extension of the fasta file ? */
		printf("Error. Cannot open \"%s\" file.\n", fname);
		exit(EXIT_FAILURE);
	}

	IGLINE=0, begline=1;
	lidx=0, mxsylen=0, mnsylen=0XFFFFFFFFFFFFFFFF;
	mxamb=0, mnamb=0xFFFFFFFF;

	sqidx=-1; /* this is slightly dangerous, you need very much to know what you're doing */
	gbuf=GBUF;
	// sqisz=malloc(gbuf*sizeof(i_s));
	sqisz=realloc(sqisz, gbuf*sizeof(i_s));
	for(i=0;i<gbuf;++i) {
		sqisz[i].ibf=GBUF;
		sqisz[i].sbf=GBUF;
		sqisz[i].id=calloc(sqisz[i].ibf, sizeof(char));
		sqisz[i].sq=calloc(sqisz[i].sbf, sizeof(char));
	}
	for(i=gbuf-GBUF;i<gbuf;++i) {
		sqisz[i].ambano[0]=0;
		sqisz[i].ambano[1]=0;
	}
	whatint=0; /* needs explanation */
	ididx0=0;

	while( ( (c = fgetc(fin)) != EOF) ) {
		if(c =='\n') {
			IGLINE=0;
			begline=1;
			lidx++;
		} else if( (begline==1) & (c == '>') ) { /* this condition catches the beginning of a new sequence, and uses it to prepare the nextsequence.*/
			IGLINE =1;
			begline=0; 
			if(sqidx>=0) { /* chancing my arm here ... operating on the past sequence */
				if(sqisz[sqidx].sylen > mxsylen)
					mxsylen = sqisz[sqidx].sylen;
				if(sqisz[sqidx].sylen < mnsylen)
					mnsylen = sqisz[sqidx].sylen;
				if(sqisz[sqidx].ambano[0] > mxamb)
					mxamb = sqisz[sqidx].ambano[0];
				if(sqisz[sqidx].ambano[0] < mnamb)
					mnamb = sqisz[sqidx].ambano[0];

				CONDREALLOC(ididx0, sqisz[sqidx].ibf, GBUF, sqisz[sqidx].id, char);
				sqisz[sqidx].id[ididx0]='\0';
				CONDREALLOC(sqisz[sqidx].sylen, sqisz[sqidx].sbf, GBUF, sqisz[sqidx].sq, char);
				sqisz[sqidx].sq[sqisz[sqidx].sylen]='\0';
				sqisz[sqidx].idz=1+ididx0;
				sqisz[sqidx].sqz=1+sqisz[sqidx].sylen;
			}

			sqidx++;
			if(sqidx==gbuf) {
				gbuf+=GBUF;
				sqisz=realloc(sqisz, gbuf*sizeof(i_s));
				for(i=gbuf-GBUF;i<gbuf;++i) {
					sqisz[i].ibf=GBUF;
					sqisz[i].sbf=GBUF;
					sqisz[i].id=calloc(sqisz[i].ibf, sizeof(char));
					sqisz[i].sq=calloc(sqisz[i].sbf, sizeof(char));
				}
			}
			sqisz[sqidx].idx=sqidx;

			/* resetting stuff */
			sqisz[sqidx].sylen=0;
			ididx0=0;
			for(i=0;i<SSZ;++i)
				sqisz[sqidx].sy[i]=0;
			for(i=0;i<2;++i)
				sqisz[sqidx].ambano[i]=0;
		} else if (IGLINE==1) {
			CONDREALLOC(ididx0, sqisz[sqidx].ibf, GBUF, sqisz[sqidx].id, char);
			sqisz[sqidx].id[ididx0++]=c;
		} else if (IGLINE==0) {
			CONDREALLOC(sqisz[sqidx].sylen, sqisz[sqidx].sbf, GBUF, sqisz[sqidx].sq, char);
			sqisz[sqidx].sq[sqisz[sqidx].sylen]=c;
			sqisz[sqidx].sylen++;
			switch(c) {
				case 'A': case 'a':
					whatint=1; break;
				case 'C': case 'c':
					whatint=2; break;
				case 'G': case 'g':
					whatint=3; break;
				case 'T': case 't':
					whatint=4; break;
				case 'R': case 'r':
					whatint=5; break;
				case 'Y': case 'y':
					whatint=6; break;
				case 'K': case 'k': /* the ketos */
					whatint=7; break;
				case 'M': case 'm': /* the aminoids */
					whatint=8; break;
				case 'S': case 's':
					whatint=9; break;
				case 'W': case 'w':
					whatint=10; break;
				case 'B': case 'b':
					whatint=11; break;
				case 'D': case 'd':
					whatint=12; break;
				case 'H': case 'h':
					whatint=13; break;
				case 'V': case 'v':
					whatint=14; break;
				case 'N': case 'n':
					whatint=15; break;
				case '-':
					whatint=16; break;
				default:
					whatint=17; /* unknown this means your fasta file is naff. */
			}
		}
		if( (whatint == 2) || (whatint == 3) ) {
			sqisz[sqidx].sy[0]++;
			sqisz[sqidx].ambano[1]++;
		} else if (whatint < 5) {
			sqisz[sqidx].sy[1]++;
			sqisz[sqidx].ambano[1]++;
		} else 
			sqisz[sqidx].ambano[0]++;
	}
	fclose(fin);
	/* postprocessing on the final sequence */
	if(sqisz[sqidx].sylen > mxsylen)
		mxsylen = sqisz[sqidx].sylen;
	if(sqisz[sqidx].sylen < mnsylen)
		mnsylen = sqisz[sqidx].sylen;
	if(sqisz[sqidx].ambano[0] > mxamb)
		mxamb = sqisz[sqidx].ambano[0];
	if(sqisz[sqidx].ambano[0] < mnamb)
		mnamb = sqisz[sqidx].ambano[0];

	/* the last sequence */
	CONDREALLOC(ididx0, sqisz[sqidx].ibf, GBUF, sqisz[sqidx].id, char);
	sqisz[sqidx].id[ididx0]='\0';
	CONDREALLOC(sqisz[sqidx].sylen, sqisz[sqidx].sbf, GBUF, sqisz[sqidx].sq, char);
	sqisz[sqidx].sq[sqisz[sqidx].sylen]='\0';
	sqisz[sqidx].idz=1+ididx0;
	sqisz[sqidx].sqz=1+sqisz[sqidx].sylen;

	numsq=sqidx+1, numano=0;
	for(i=0;i<numsq;++i) {
		if(sqisz[i].ambano[1])
			numano++;
	}
	// return 0;
	for(i=numsq;i<gbuf;++i) {
		free(sqisz[i].id);
		free(sqisz[i].sq);
	}
	sqisz=realloc(sqisz, numsq*sizeof(i_s));

	/* check for uniform sequence size, necessary for alignments */
	uo_t *uov=uniquelens(sqisz, numsq);
	prtuo(uov);

	prtsq(sqisz, numsq);

	free(uov->ua);
	free(uov->oa);
	free(uov);
	*nsq=numsq;
	return sqisz;
}

int main(int argc, char **argv)
{
    unsigned char *name_table;
    char find_all; /* boolean to see if -g option is called or not */
    int namecount;
    int name_entry_size;
    int ovector[OVECCOUNT]; /* very important, this is for the substring ... the part of the targetxt which matches the pattern */
    /* the ovec will appear to hold the start and end byte of a match, with reference to the total targetxt. The third of the triplets is the lenght
     * qhich appears to be put in there for convenience purposes */
    int ttxtlen;
    int i;

	/* by default, find all matches by default (originally a -g option) */
    find_all = 1;

    /* After the options, we require exactly two arguments, which are the pattern, and the targetxt string. */
    if (argc != 3) {
        printf("Two arguments required: a regex and a fasta file name: for now only first sequence will be analysed.\n");
        exit(EXIT_FAILURE);
    }

    char *pattern = argv[1]; /* OK, here is our pattern */
	unsigned numsq;
    i_s *sqisz=procfastaf(argv[2], &numsq); /* and here is our target text:sqisz[0].sq is the first sequence */
    char *targetxt = sqisz[0].sq;
    ttxtlen = sqisz[0].sylen;

    /* OK, now we're going for our first pcre function call, we'll need twwo error variables. Pattern gets compiled into "re" here and errors can be captured */
    const char *error;
    int erroffset;
    pcre *re = pcre_compile(
            pattern,              /* the pattern */
            0,                    /* default options */
            &error,               /* for error message */
            &erroffset,           /* for error offset */
            NULL);                /* use default character tables */

    /* Handle error, if it occurs: compilation has failed: print the error message and exit */
    if (re == NULL) {
        printf("PCRE compilation failed at offset %d: %s\n", erroffset, error);
        exit(EXIT_FAILURE);
    }

    /* If the compilation succeeded, we can go on to execute. in order to do a     *
     * pattern match against the targetxt string. This does just ONE match. If *
     * further matching is needed, it will be done later. the returned integer
	reflects the amoutn of free space in the output vector because after the first three entries
 (which are the full pattern match itself), the subgroups are all stored. */
    int rc = pcre_exec(
            re,                   /* the compiled pattern */
            NULL,                 /* no extra data - we didn't study the pattern */
            targetxt,              /* the targetxt string */
            ttxtlen,       /* the length of the targetxt */
            0,                    /* start at offset 0 in the targetxt */
            0,                    /* default options */
            ovector,              /* output vector for substring information */
            OVECCOUNT);           /* number of elements in the output vector */

    /* As usual, a call to pcre will require we look out for errors if they have occurred, This time, "no match" is treated as an error, although it's fairly valid as a result. */
    if (rc < 0) {
        switch(rc) {
            case PCRE_ERROR_NOMATCH:
                printf("No match\n"); break;
            /* Handle other special cases if you like */
            default: printf("Matching error %d\n", rc); break;
        }
        pcre_free(re);     /* With an error we should bail out gracefully, that is release memory used for the compiled pattern */
        exit(EXIT_FAILURE);
    }

    /* Match succeded  ok we can continue */
    printf("\nMatch succeeded at offset %i, ov1= %i.", ovector[0], ovector[1]);
    printf("Return val from pcre_exec was %i\n", rc);

    /* We have found the first match within the targetxt string. Now that could have included substrings don't forget.
	If the output *
     * vector wasn't big enough, say so. Then output any substrings that were captured. */
    /* The output vector wasn't big enough */
    if (rc == 0) {
        rc = OVECCOUNT/3;
        printf("ovector only has room for %d captured substrings\n", rc - 1);
    }

    /* Show substrings stored in the output vector by number. Obviously, in a real
       application you might want to do things other than print them. */
    for (i = 1; i < rc; i++) { /* i=0 is the string for the full match ... but once we know there's been a full match, hardly worth the bother printing it. */
        char *substring_start = targetxt + ovector[2*i];
        int substring_length = ovector[2*i+1] - ovector[2*i];
        printf("%2d: %.*s\n", i, substring_length, substring_start);
    }

    /**************************************************************************
     * That concludes the basic part of this demonstration program. We have    *
     * compiled a pattern, and performed a single match. The code that follows *
     * shows first how to access named substrings, and then how to code for    *
     * repeated matches on the same targetxt.                                   *
     **************************************************************************/

    /* See if there are any named substrings, and if so, show them by name. First
       we have to extract the count of named parentheses from the pattern. */
    (void)pcre_fullinfo(
            re,                   /* the compiled pattern */
            NULL,                 /* no extra data - we didn't study the pattern */
            PCRE_INFO_NAMECOUNT,  /* number of named substrings */
            &namecount);          /* where to put the answer */

    if (namecount <= 0)
        printf("No named substrings\n");
    else {
        unsigned char *tabptr;
        printf("Named substrings\n");

        /* Before we can access the substrings, we must extract the table for
           translating names to numbers, and the size of each entry in the table. */
        (void)pcre_fullinfo(
                re,                       /* the compiled pattern */
                NULL,                     /* no extra data - we didn't study the pattern */
                PCRE_INFO_NAMETABLE,      /* address of the table */
                &name_table);             /* where to put the answer */

        (void)pcre_fullinfo(
                re,                       /* the compiled pattern */
                NULL,                     /* no extra data - we didn't study the pattern */
                PCRE_INFO_NAMEENTRYSIZE,  /* size of each entry in the table */
                &name_entry_size);        /* where to put the answer */

        /* Now we can scan the table and, for each entry, print the number, the name,
           and the substring itself. */
        tabptr = name_table;
        for (i = 0; i < namecount; i++) {
            int n = (tabptr[0] << 8) | tabptr[1];
            printf("(%d) %*s: %.*s\n", n, name_entry_size - 3, tabptr + 2,
                    ovector[2*n+1] - ovector[2*n], targetxt + ovector[2*n]);
            tabptr += name_entry_size;
        }
    }

    if (!find_all) {
        pcre_free(re);   /* Release the memory used for the compiled pattern */
        return 0;        /* Finish unless -g was given */
    }

    /*************************************************************************
     * If the "-g" option was given on the command line, we want to continue  *
     * to search for additional matches in the targetxt string, in a similar   *
     * way to the /g option in Perl. This turns out to be trickier than you   *
     * might think because of the possibility of matching an empty string.    *
     * What happens is as follows:                                            *
     *                                                                        *
     * If the previous match was NOT for an empty string, we can just start   *
     * the next match at the end of the previous one.                         *
     *                                                                        *
     * If the previous match WAS for an empty string, we can't do that, as it *
     * would lead to an infinite loop. Instead, a special call of pcre_exec() *
     * is made with the PCRE_NOTEMPTY_ATSTART and PCRE_ANCHORED flags set.    *
     * The first of these tells PCRE that an empty string at the start of the *
     * targetxt is not a valid match; other possibilities must be tried. The   *
     * second flag restricts PCRE to one match attempt at the initial string  *
     * position. If this match succeeds, an alternative to the empty string   *
     * match has been found, and we can proceed round the loop.               *
     *************************************************************************/
    for (;;) { /* Infinite loop for second and subsequent matches */
        int options = 0;                 /* Normally no options */
        int start_offset = ovector[1];   /* Start at end of previous match */

        /* If the previous match was for an empty string, we are finished if we are
           at the end of the targetxt. Otherwise, arrange to run another match at the
           same point to see if a non-empty match can be found. */

        if (ovector[0] == ovector[1]) {
            if (ovector[0] == ttxtlen) break;
//            options = PCRE_NOTEMPTY_ATSTART | PCRE_ANCHORED;
            options = PCRE_NOTEMPTY | PCRE_ANCHORED;
        }

        /* Run the next matching operation */
        rc = pcre_exec(
                re,                   /* the compiled pattern */
                NULL,                 /* no extra data - we didn't study the pattern */
                targetxt,              /* the targetxt string */
                ttxtlen,       /* the length of the targetxt */
                start_offset,         /* starting offset in the targetxt */
                options,              /* options */
                ovector,              /* output vector for substring information */
                OVECCOUNT);           /* number of elements in the output vector */

        /* This time, a result of NOMATCH isn't an error. If the value in "options"
           is zero, it just means we have found all possible matches, so the loop ends.
           Otherwise, it means we have failed to find a non-empty-string match at a
           point where there was a previous empty-string match. In this case, we do what
           Perl does: advance the matching position by one, and continue. We do this by
           setting the "end of previous match" offset, because that is picked up at the
           top of the loop as the point at which to start again. */
        if (rc == PCRE_ERROR_NOMATCH) {
            if (options == 0)
                break;
            ovector[1] = start_offset + 1;
            continue;    /* Go round the loop again */
        }

        /* Other matching errors are not recoverable. */
        if (rc < 0) {
            printf("Matching error %d\n", rc);
            pcre_free(re);    /* Release memory used for the compiled pattern */
            return 1;
        }

        /* Match succeded */
        printf("\nMatch succeeded again at offset %d\n", ovector[0]);

        /* The match succeeded, but the output vector wasn't big enough. */
        if (rc == 0) {
            rc = OVECCOUNT/3;
            printf("ovector only has room for %d captured substrings\n", rc - 1);
        }

        /* As before, show substrings stored in the output vector by number, and then
           also any named substrings. */
        for (i = 0; i < rc; i++) {
            char *substring_start = targetxt + ovector[2*i];
            int substring_length = ovector[2*i+1] - ovector[2*i];
            printf("%2d: %.*s\n", i, substring_length, substring_start);
        }

        if (namecount <= 0)
            printf("No named substrings\n");
        else {
            unsigned char *tabptr = name_table;
            printf("Named substrings\n");
            for (i = 0; i < namecount; i++) {
                int n = (tabptr[0] << 8) | tabptr[1];
                printf("(%d) %*s: %.*s\n", n, name_entry_size - 3, tabptr + 2, ovector[2*n+1] - ovector[2*n], targetxt + ovector[2*n]);
                tabptr += name_entry_size;
            }
        }
    } /* End of loop to find second and subsequent matches */

    printf("\n");
    pcre_free(re); /* Release memory used for the compiled pattern */
	for(i=0;i<numsq;++i) {
		free(sqisz[i].id);
		free(sqisz[i].sq);
	}
	free(sqisz);
    return 0;
}
