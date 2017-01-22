/* example file, canonical
 Substrings: to you and me, this means groups, the stuff which you put in parenteses in your regex.
You can name a substring by following the first bracket with a question mark and embed the name in <>'s.

working with 
./dem5 "-*\d+\.\d+,\d+\.\d+,\d+\.\d+" alcabre-cuvi.kml "LineString><coordinates>"  "</coordinates></LineString"
on the kml file

*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <pcre.h>

#define OVECCOUNT 300    /* should be a multiple of 3 */
#define MXFSZ (1<<16)
#define BUF 4 /* tested for a low a value as 2, and is fine */

size_t fszfind(FILE *fp)
{
    fseek(fp, 0, SEEK_END);
    size_t fbytsz = (size_t)ftell(fp);
    rewind(fp);
    return fbytsz;
}

typedef struct
{
    char *ss; /* start string */
    size_t ssz; /* start string sz*/
    char *es; /* seq */
    size_t esz; /* seq sz*/
} sandestr; /* start and end strng struct */

char *scas(char *s, size_t *sz) /* a function to scan the token strings and convert newlines double symbols to their proper 0x0A equivalents. Also the size is returned in the arguments */
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
    *sz=(size_t)strlen(ts);
    return ts;
}

char *retmat_frfile(char *fname, sandestr *sestr, size_t *stsz) /* return matrix from file */
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

    /* so now we can close out file */
    fclose(fp);

    /*  now we get a pointer which will scan the big string */
    char seenss=0;
    char *p=fslurp;
    char *bp=NULL; // Ok, this is the back pointer that is 8 ssz characters back from p
    char *stcp=calloc(sestr->ssz+1, sizeof(char));
    char *etcp=calloc(sestr->esz+1, sizeof(char));
    size_t st_b=0, ed_b=0; // start byte end byte, important to initialize them */

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
                continue;
        }
        /* the endingrng we're interested in */
        if( seenss & (*(bp+sestr->esz-3)==*(sestr->es+sestr->esz-3)) & (*(bp+sestr->esz-2)==*(sestr->es+sestr->esz-2)) & (*(bp+sestr->esz-1)==*(sestr->es+sestr->esz-1)) ) { // if these 3 chars are matched, thing look promising, but not guaranteed. 
            memcpy(etcp, bp, sestr->esz);
            if( (strcmp(etcp, sestr->es)) == 0 ) {
                ed_b=(size_t)(bp - fslurp); // double verified
                seenss=0;
                break;
            } else
                continue;
        }
    }
    size_t mstsz = ed_b - st_b;
    char *matstr=calloc(mstsz+1, sizeof(char));
    strncpy(matstr, fslurp+st_b, mstsz);

    p=NULL;
    bp=NULL;
    free(stcp);
    free(etcp);
    free(fslurp);

    *stsz=mstsz;
    return matstr;
}

int main(int argc, char **argv)
{
    unsigned char *name_table;
    int namecount;
    int name_entry_size;
    int ovector[OVECCOUNT]; /* very important, this is for the substring ... the part of the targetxt which matches the pattern */
    /* the ovec will appear to hold the start and end byte of a match, with reference to the total targetxt. The third of the triplets is the lenght
     * qhich appears to be put in there for convenience purposes */
    int ttxtlen;
    int i;


    /* argument accounting */
    if(argc!=5) {
        printf("Error. Pls supply 4 arguments (regex, name of text file, startstringtoken and endstringtoken.\n");
        exit(EXIT_FAILURE);
    }

    char *pattern = argv[1]; /* OK, here is our pattern */
    /* allocate and populate the start-end-string struct properly */
    sandestr *sestr=malloc(sizeof(sandestr)); /* declare start-end-string struct */
    sestr->ss=scas(argv[3], &(sestr->ssz)); /* use the scas function to properly represent newliens and get the two string sizes. */
    sestr->es=scas(argv[4], &(sestr->esz));
    size_t ttsz=0;
    char *targetxt=retmat_frfile(argv[2], sestr, &ttsz);

    ttxtlen = (int)ttsz;
    /* OK, now we're goign for our first pcre function call, we'll need two error variables. Pattern gets compiled into "re" here and errors can be captured */
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
    printf("Match at offset %i. ", ovector[0]);
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

    /* OK, first match exploration undertaken, next we look for subsequent matches: we use and infinite loop */
    for (;;) { /* Infinite loop for second and subsequent matches */
        int options = 0;                 /* Normally no options */
        int start_offset = ovector[1];   /* Start at end of previous match */

        /* If the previous match was for an empty string, we are finished if we are
           at the end of the targetxt. Otherwise, arrange to run another match at the
           same point to see if a non-empty match can be found. */

        if (ovector[0] == ovector[1]) {
            if (ovector[0] == ttxtlen)
                break;
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
        printf("Match at offset %d: ", ovector[0]);

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
            printf("%.*s\n", substring_length, substring_start);
        }

        if (namecount > 0) {
            unsigned char *tabptr = name_table;
            printf("Named substrings\n");
            for (i = 0; i < namecount; i++) {
                int n = (tabptr[0] << 8) | tabptr[1];
                printf("(%d) %*s: %.*s\n", n, name_entry_size - 3, tabptr + 2, ovector[2*n+1] - ovector[2*n], targetxt + ovector[2*n]);
                tabptr += name_entry_size;
            }
        }
    } /* End of loop to find second and subsequent matches */

    pcre_free(re); /* Release memory used for the compiled pattern */
    free(targetxt);
    free(sestr->ss);
    free(sestr->es);
    free(sestr);
    return 0;
}
