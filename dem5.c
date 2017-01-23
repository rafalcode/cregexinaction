/* example file, canonical
 Substrings: to you and me, this means groups, the stuff which you put in parenteses in your regex.
You can name a substring by following the first bracket with a question mark and embed the name in <>'s.

working with 
./dem5 "-*\d+\.\d+,\d+\.\d+,\d+\.\d+" alcabre-cuvi.kml "LineString><coordinates>"  "</coordinates></LineString"
on the kml file


By way of explanation: pcre.
ovector holds a single match, followed by the subgroups (denoted by parentheses) within the selfsame match
Both the full match and its subgroups required three numbers: start, end and length.
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <pcre.h>

#define OVECCOUNT 300    /* should be a multiple of 3 */
#define MXFSZ (1<<16)
#define BUF 4 /* tested for a low a value as 2, and is fine */

#define CONDREALLOC(x, b, c, a, t); \
    if((x)>=((b)-0)) { \
        (b) += (c); \
        (a)=realloc((a), (b)*sizeof(t)); \
    }

typedef struct /* gpx_t */
{
    float c[3]; // lat, lon, ele;
    int dt[6]; // yr, mon, day, hr, min, sec;
} gpx_t;

typedef struct /* sandestr: start and end string, this delimits and interior section of the file */
{
    char *ss; /* start string */
    int ssz; /* start string sz*/
    char *es; /* seq */
    int esz; /* seq sz*/
} sandestr; /* start and end strng struct */

typedef struct /* strblk_t: str blk type */
{
    char *sb; /* string block */
    int sbsz; /* start string sz*/
} strblk_t; /* start and end strng struct */

size_t fszfind(FILE *fp)
{
    fseek(fp, 0, SEEK_END);
    size_t fbytsz = (size_t)ftell(fp);
    rewind(fp);
    return fbytsz;
}

char *scas(char *s, int *sz) /* a function to scan the token strings and convert newlines double symbols to their proper 0x0A equivalents. Also the size is returned in the arguments */
{
    int j=0;
    size_t tsz=strlen(s);
    char *p=s, *ts= calloc(tsz+1, sizeof(char));
    char oldc='\0';
    int i;
    for(i=0; i<tsz; i++) {
        if(oldc != '\\') {
            if( *(p+i) != '\\') {
                *(ts+j) = *(p+i);
                j++;
            } else {
                oldc=*(p+i);
                continue;
            }
        } else {
            switch(*(p+i)) {
                case 'n':
                    *(ts+j) ='\n';
                    j++;
                    oldc='\n';
                    continue;
                default: /* placeholder for other symbols */
                    printf("Error. Only newline double symbols allowed.\n");
                    free(ts);
                    exit(EXIT_FAILURE);
            }
        }
    }
    *(ts+j)='\0';
    ts=realloc(ts, (j+1)*sizeof(char));
    *sz=(int)strlen(ts);
    return ts;
}

strblk_t *retmat_frfile(char *fname, sandestr *sestr, int *sbasz) /* return matrix from file */
{
    /* declarations */
    int i;
    FILE *fp=fopen(fname,"r");
    if(fp == NULL) {
        printf("The file you specified could not be opened. Wrong name, wrong directory?\n"); 
        exit(EXIT_FAILURE);
    }
    int fsz = fszfind(fp); /*  this size will not include the EOF, but yes,the final \n? */
    if(fsz> MXFSZ) {
        printf("Error. At %i input file size was too big by %3.4f%%. Sorry.\n", fsz, (float)(fsz - MXFSZ)/MXFSZ);
        exit(EXIT_FAILURE);
    }

    char *fslurp = malloc(sizeof(char)*(fsz+1));
    fread(fslurp, sizeof(char), (size_t)fsz, fp);
    fslurp[fsz]='\0'; /* makes sure fslurp is a string */
    fclose(fp); /* so now we can close out file */

    /*  now we get a pointer which will scan the big string */
    char seenss=0; /* marker for having seem the startstring */
    int nseens=0 /* Number of seen blocks */, oldnseens=0 /* previous value of it. */;
    char *p=fslurp; /* as it says, we shall slurp the input file into here */
    char *bp=NULL; // Ok, this is the back pointer that is 8 ssz characters back from p
    char *stcp=calloc(sestr->ssz+1, sizeof(char));
    char *etcp=calloc(sestr->esz+1, sizeof(char));
    size_t st_b=0, ed_b=0; // start byte end byte, important to initialize them */

    int sbbuf=BUF;
    strblk_t *sbt=malloc(sbbuf*sizeof(strblk_t));
    for(i=0;i<sbbuf;++i) {
        sbt[i].sb=NULL;
        sbt[i].sbsz=0;
    }

    for(i=0; i<fsz; i++) {
        /* assign the back pointer bp */
        switch(seenss) {
            case 0:
                if( ((p+i)-fslurp) < sestr->ssz)
                    continue;
                else if( ((p+i)-fslurp) >= sestr->ssz)
                    bp = (p+i)-sestr->ssz;
            case 1:
                if( ((p+i)-(fslurp+st_b)) < sestr->esz)
                    continue;
                else if( ((p+i)-(fslurp+st_b)) >= sestr->esz)
                    bp = (p+i)-sestr->esz;
        }

        /* the starting strng we're interested in */
        if( (seenss == 0) & (*(bp+sestr->ssz-3)==*(sestr->ss+sestr->ssz-3)) & (*(bp+sestr->ssz-2)==*(sestr->ss+sestr->ssz-2)) & (*(bp+sestr->ssz-1)==*(sestr->ss+sestr->ssz-1)) ) { // if these 3 chars are matched, thing look promising, but not guaranteed. 
            memcpy(stcp, bp, sestr->ssz);
            if( (strcmp(stcp, sestr->ss)) == 0 ) { /* got it, full start string match */
                // st_b=(size_t)((p+i)-fslurp);
                st_b=(size_t)((bp+sestr->ssz)-fslurp);
                seenss=1;
            } else
                continue; // false alarm
        }
        /* on the ending side, we're interested in */
        if( seenss & (*(bp+sestr->esz-3)==*(sestr->es+sestr->esz-3)) & (*(bp+sestr->esz-2)==*(sestr->es+sestr->esz-2)) & (*(bp+sestr->esz-1)==*(sestr->es+sestr->esz-1)) ) { // if these 3 chars are matched, thing look promising, but not guaranteed. 
            memcpy(etcp, bp, sestr->esz);
            if( (strcmp(etcp, sestr->es)) == 0 ) {
                ed_b=(size_t)(bp - fslurp); // double verified
                seenss=0;
                nseens++; /*Ok, we have just been though one block */
                // break; // I leave this here because it was there originally and caused quite a headache
            } else
                continue;
        }
        if(nseens != oldnseens) {
            CONDREALLOC(nseens, sbbuf, BUF, sbt, strblk_t);
            sbt[nseens-1].sbsz = ed_b - st_b;
            sbt[nseens-1].sb=realloc(sbt[nseens-1].sb, (sbt[nseens-1].sbsz+1)*sizeof(char));
            strncpy(sbt[nseens-1].sb, fslurp+st_b, sbt[nseens-1].sbsz);
            sbt[nseens-1].sb[sbt[nseens-1].sbsz]='\0';
#if defined DBG
            printf("Found blocksz=%d\n", sbt[nseens-1].sbsz);
            // printf("%s\n", sbt[nseens-1].sb);
#endif
            oldnseens=nseens;
        }
    }

    p=NULL;
    bp=NULL;
    *sbasz=nseens;
    sbt=realloc(sbt, (*sbasz)*sizeof(strblk_t));
    free(stcp);
    free(etcp);
    free(fslurp);

    return sbt;
}

int main(int argc, char **argv)
{
    unsigned char *name_table=NULL;
    int namecount;
    int name_entry_size;
    // int *ovector=calloc(OVECCOUNT, sizeof(int)); /* very important, this is for the substring ... the part of the targetxt which matches the pattern */
    int ovector[OVECCOUNT]={0};
    /* the ovec will appear to hold the start and end byte of a match, with reference to the total targetxt. The third of the triplets is the lenght
     * qhich appears to be put in there for convenience purposes */
    int i, j;

    /* argument accounting */
    if(argc!=5) {
        printf("Error. Pls supply 4 arguments (regex, name of text file, startstringtoken and endstringtoken.\n");
        printf("Typical usage on a KML file:\n");
        printf("./dem5 \"-*\\d+\\.\\d+,\\d+\\.\\d+,\\d+\\.\\d+\" alcabre-cuvi.kml \"LineString><coordinates>\"  \"</coordinates></LineString\"\n");
        //or
//        ./dem5_d 'lat="(-*\d+\.\d+)" lon="(-*\d+\.\d+)">\n\ +<ele>(\d+.\d+)</ele>\n\ +<time>(\d+)-(\d+)-(\d+)T(\d+):(\d+):(\d+)Z' ill.gpx "<trkseg>"  "</trkseg"
//        //for gpx
        exit(EXIT_FAILURE);
    }

    gpx_t *dholder=NULL;
    int dhbuf=BUF;
    unsigned char GPXYES=0;
    char *tstr=strrchr(argv[2], '.');
    if(!strcmp(tstr+1, "gpx")) {
            GPXYES=1;
            dholder=calloc(dhbuf, sizeof(gpx_t));
            // this is the famed derivative use of void pointer ... here is automaticaly gets assigned to gpx_t type!
    }

    char *pattern = argv[1]; /* OK, here is our pattern */
    /* allocate and populate the start-end-string struct properly */
    sandestr *sestr=malloc(sizeof(sandestr)); /* declare start-end-string struct */
    sestr->ss=scas(argv[3], &(sestr->ssz)); /* use the scas function to properly represent newliens and get the two string sizes. */
    sestr->es=scas(argv[4], &(sestr->esz));
    int sbasz=0; /// str block array size
    strblk_t *sbt=retmat_frfile(argv[2], sestr, &sbasz);
#if defined DBG
    printf("Interpreted sizes, startstr: %d, endstr: %d\n", sestr->ssz, sestr->esz); 
    printf("Numblocks=%d\n", sbasz);
#endif
    //char *targetxt=retmat_frfile(argv[2], sestr, &ttsz);

    /* OK, now we're goign for our first pcre function call, we'll need two error variables. Pattern gets compiled into "re" here and errors can be captured */
    pcre *re = NULL;
    int rc;
    const char *error;
    int erroffset;
    char *substring_start = NULL;
    int substring_length;
    unsigned char *tabptr=NULL;
    int n;
    int options;
    int start_offset;
    int nmatches=0;
    char tmpstr[128]={0};


    re = pcre_compile(
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

    for(j=0;j<sbasz;++j) {
        /* If the compilation succeeded, we can go on to execute. in order to do a     *
         * pattern match against the targetxt string. This does just ONE match. If *
         * further matching is needed, it will be done later. the returned integer
         reflects the amount of free space in the output vector because after the first three entries
         (which are the full pattern match itself), the subgroups are all stored. */
#if defined DBG
        // printf("%s\n", sbt[j].sb);
#endif
        rc = pcre_exec(
                re,                   /* the compiled pattern */
                NULL,                 /* no extra data - we didn't study the pattern */
                sbt[j].sb,              /* the targetxt string */
                sbt[j].sbsz,       /* the length of the targetxt */
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
        printf("Match at offset %i. ", ovector[0]);
        printf("Return val from pcre_exec was %i\n", rc);
        nmatches++;
        if(GPXYES)
            CONDREALLOC(nmatches, dhbuf, BUF, dholder, gpx_t);

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
            substring_start = sbt[j].sb + ovector[2*i];
            substring_length = ovector[2*i+1] - ovector[2*i];
            printf("%2d: %.*s\n", i, substring_length, substring_start);
            sprintf(tmpstr, "%.*s", substring_length, substring_start);
            if(GPXYES) {
                if(i<4)
                    dholder[nmatches-1].c[i-1]=atof(tmpstr);
                else if ((i>=4) &(i<10))
                    dholder[nmatches-1].dt[i-3-1]=atoi(tmpstr);
                else
                    printf("Shouldn't be getting more than 9 subgroup matches\n"); 
            }
        }

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

            /* Now we can scan the table and, for each entry, print the number, the name, and the substring itself. */
            tabptr = name_table;
            for (i = 0; i < namecount; i++) {
                n = (tabptr[0] << 8) | tabptr[1];
                printf("(%d) %*s: %.*s\n", n, name_entry_size - 3, tabptr + 2, ovector[2*n+1] - ovector[2*n], sbt[j].sb + ovector[2*n]);
                tabptr += name_entry_size;
            }
        }

        /* OK, first match exploration undertaken, next we look for subsequent matches: we use and infinite loop */
        for (;;) { /* Infinite loop for second and subsequent matches */
            options = 0;                 /* Normally no options */
            start_offset = ovector[1];   /* Start at end of previous match */

            /* If the previous match was for an empty string, we are finished if we are
               at the end of the targetxt. Otherwise, arrange to run another match at the
               same point to see if a non-empty match can be found. */

            if (ovector[0] == ovector[1]) {
                if (ovector[0] == sbt[j].sbsz)
                    break;
                //            options = PCRE_NOTEMPTY_ATSTART | PCRE_ANCHORED;
                options = PCRE_NOTEMPTY | PCRE_ANCHORED;
            }

            /* Run the next matching operation */
            rc = pcre_exec(
                    re,                   /* the compiled pattern */
                    NULL,                 /* no extra data - we didn't study the pattern */
                    sbt[j].sb,              /* the targetxt string */
                    sbt[j].sbsz,       /* the length of the targetxt */
                    start_offset,         /* starting offset in the targetxt */
                    options,              /* options */
                    ovector,              /* output vector for substring information */
                    OVECCOUNT);           /* number of elements in the output vector */

            if (rc == PCRE_ERROR_NOMATCH) { /* all possible matches found */
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
            printf("Match at offset %d; groups= ", ovector[0]);
            nmatches++;
            if(GPXYES)
                CONDREALLOC(nmatches, dhbuf, BUF, dholder, gpx_t);

            /* The match succeeded, but the output vector wasn't big enough. */
            if (rc == 0) {
                rc = OVECCOUNT/3;
                printf("ovector only has room for %d captured substrings\n", rc - 1);
            }

            /* As before, show substrings stored in the output vector by number, and then
               also any named substrings. */
            for (i = 1; i < rc; i++) { // avoid first three entries because these are the whole match
                substring_start = sbt[j].sb + ovector[2*i];
                substring_length = ovector[2*i+1] - ovector[2*i];
                printf("%d) %.*s ", i, substring_length, substring_start);
                sprintf(tmpstr, "%.*s", substring_length, substring_start);
                if(GPXYES) {
                    if(i<4)
                        dholder[nmatches-1].c[i-1]=atof(tmpstr);
                    else if ((i>=4) &(i<10))
                        dholder[nmatches-1].dt[i-3-1]=atoi(tmpstr);
                    else
                        printf("Shouldn't be getting more than 9 subgroup matches\n"); 
                }
            }
            printf("\n"); 

            if (namecount > 0) {
                tabptr = name_table;
                printf("Named substrings\n");
                for (i = 0; i < namecount; i++) {
                    n = (tabptr[0] << 8) | tabptr[1];
                    printf("(%d) %*s: %.*s\n", n, name_entry_size - 3, tabptr + 2, ovector[2*n+1] - ovector[2*n], sbt[j].sb + ovector[2*n]);
                    tabptr += name_entry_size;
                }
            }
        } /* End of loop to find second and subsequent matches */
        /* reset all the regex stuff */
        memset(ovector,0,OVECCOUNT*sizeof(int));
        tabptr=NULL;
        printf("Nummatches=%d\n", nmatches); 
        if(GPXYES) {
            for(i=0;i<nmatches;++i) {
                for(j=0;j<3;++j) 
                    printf("%4.6f ", dholder[i].c[j]); 
                for(j=3;j<9;++j) 
                    printf("%d ", dholder[i].dt[j-3]); 
                printf("\n"); 
            }
        }
    } /* end of loop for the different strblks */
    pcre_free(re);

    /* OK, closing down ... free up everything */
    for(i=0;i<sbasz;++i) 
        free(sbt[i].sb);
    if(GPXYES)
        free(dholder);
    free(sbt);
    free(sestr->ss);
    free(sestr->es);
    free(sestr);
    return 0;
}
