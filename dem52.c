/*************************************************
*           PCRE2 DEMONSTRATION PROGRAM          *
*************************************************/

/* This is a demonstration program to illustrate a straightforward way of
using the PCRE2 regular expression library from a C program. See the
pcre2sample documentation for a short discussion ("man pcre2sample" if you have
the PCRE2 man pages installed). PCRE2 is a revised API for the library, and is
incompatible with the original PCRE API.

There are actually three libraries, each supporting a different code unit
width. This demonstration program uses the 8-bit library. The default is to
process each code unit as a separate character, but if the pattern begins with
"(*UTF)", both it and the subject are treated as UTF-8 strings, where
characters may occupy multiple code units.

In Unix-like environments, if PCRE2 is installed in your standard system
libraries, you should be able to compile this program using this command:

cc -Wall pcre2demo.c -lpcre2-8 -o pcre2demo

If PCRE2 is not installed in a standard place, it is likely to be installed
with support for the pkg-config mechanism. If you have pkg-config, you can
compile this program using this command:

cc -Wall pcre2demo.c `pkg-config --cflags --libs libpcre2-8` -o pcre2demo

If you do not have pkg-config, you may have to use something like this:

cc -Wall pcre2demo.c -I/usr/local/include -L/usr/local/lib \
  -R/usr/local/lib -lpcre2-8 -o pcre2demo

Replace "/usr/local/include" and "/usr/local/lib" with wherever the include and
library files for PCRE2 are installed on your system. Only some operating
systems (Solaris is one) use the -R option.

Building under Windows:

If you want to statically link this program against a non-dll .a file, you must
define PCRE2_STATIC before including pcre2.h, so in this environment, uncomment
the following line. */

/* #define PCRE2_STATIC */

/* The PCRE2_CODE_UNIT_WIDTH macro must be defined before including pcre2.h.
For a program that uses only one code unit width, setting it to 8, 16, or 32
makes it possible to use generic function names such as pcre2_compile(). Note
that just changing 8 to 16 (for example) is not sufficient to convert this
program to process 16-bit characters. Even in a fully 16-bit environment, where
string-handling functions such as strcmp() and printf() work with 16-bit
characters, the code for handling the table of named substrings will still need
to be modified. */

#define PCRE2_CODE_UNIT_WIDTH 8

#include <stdio.h>
#include <string.h>
#include <pcre2.h>
#include <math.h>
#include <sys/stat.h>

#define MYRGX "(\\d+):(\\d+):(\\d+)"
#define GBUF 5
#define CONDREALLOC0(x, b, c, a, t); \
    if((x)>=((b)-1)) { \
        (b) += (c); \
        (a)=realloc((a), (b)*sizeof(t)); \
        memset(((a)+((b)-(c))), 0, (c)*sizeof(t)); \
    }

#define CONDREALLOC(x, b, c, a, t); \
    if((x)>=((b)-1)) { \
        (b) += (c); \
        (a)=realloc((a), (b)*sizeof(t)); \
        for(int i=((b)-(c)); i<(b); ++i) { \
            (a)[i].h = 0; \
            (a)[i].m = 0; \
            (a)[i].s = 0; \
            (a)[i].hh = 0; \
        } \
    }

// following is like atoi (natoi), but does not require \0 at end, only a simplesequence of chars of size n
#define NATOI(x, cs, n); \
    for(int i=0; i<(n); ++i) \
        (x) += ((int)(cs)[i]-48) * pow(10, ((n)-1-(i)));

/**************************************************************************
* Here is the program. The API includes the concept of "contexts" for     *
* setting up unusual interface requirements for compiling and matching,   *
* such as custom memory managers and non-standard newline definitions.    *
* This program does not do any of this, so it makes no use of contexts,   *
* always passing NULL where a context could be given.                     *
**************************************************************************/

typedef struct /* tt_c */
{
    int h, m, s, hh;
} tt_c;

typedef struct /* av_c */
{
    int vbf, vsz;
    tt_c *v;
} av_c;

av_c *crea_avc(int vbf)
{
    av_c *avc=malloc(sizeof(av_c));
    avc->vbf=vbf;
    avc->v=calloc(avc->vbf, sizeof(tt_c));
    avc->vsz=0;
    return avc;
}

void condrea_avc(av_c *avc)
{
    /* somewhat trivial, but idea is that, as avc is a container, it can be re-alloced inside a function */
    CONDREALLOC(avc->vsz, avc->vbf, GBUF, avc->v, tt_c);
    return;
}

void norm_avc(av_c *avc)
{
    avc->v=realloc(avc->v, avc->vsz*sizeof(tt_c));
    return;
}

void free_avc(av_c *avc)
{
    free(avc->v);
    free(avc);
    return;
}

void prtavec(av_c *avc)
{
    int i;
    for(i=0;i<avc->vsz;++i)
        printf("%i:%i:%i:%i\n", avc->v[i].h, avc->v[i].m, avc->v[i].s, avc->v[i].hh); 
    return;
}

int main(int argc, char **argv)
{
    if(argc !=2) {
        printf("The RGX is hardcoded. One argument only: the fil in which oto find the RGX.\n");
        exit(EXIT_FAILURE);
    }

    pcre2_code *re;
    PCRE2_SPTR pattern;     /* PCRE2_SPTR is a pointer to unsigned code units of */
    PCRE2_SPTR subject;     /* the appropriate width (in this case, 8 bits). */

    int crlf_is_newline;
    int errornumber;
    int find_all = 1; // this is for the -g option, which I generally want
    int i;
    int rc;
    int utf8;

    uint32_t option_bits;
    uint32_t namecount;
    uint32_t newline;

    PCRE2_SIZE erroroffset;
    PCRE2_SIZE *ovector;
    PCRE2_SIZE subject_length;

    pcre2_match_data *match_data;

    struct stat fsta;
    if(stat(argv[1], &fsta) == -1) {
        fprintf(stderr,"Can't open input file %s", argv[1]);
        exit(EXIT_FAILURE);
    }
#ifdef DBG
    printf("tot filesize=%zu\n", fsta.st_size); 
#endif
    FILE *fp=fopen(argv[1], "r");
    char *fslurp = malloc(sizeof(char)*(fsta.st_size+1));
    fread(fslurp, sizeof(char), fsta.st_size, fp);
    fslurp[fsta.st_size]='\0'; /* makes sure fslurp is a string */
    fclose(fp); /* so now we can close out file */

    /* Pattern and subject are char arguments, so they can be straightforwardly
       cast to PCRE2_SPTR because we are working in 8-bit code units. The subject
       length is cast to PCRE2_SIZE for completeness, though PCRE2_SIZE is in fact
       defined to be size_t. */
    pattern = (PCRE2_SPTR)MYRGX;
    // pattern = (PCRE2_SPTR)"[0-9]+:";
    // pattern = (PCRE2_SPTR)"\\d+:\\d+:\\d+";
    subject = (PCRE2_SPTR)fslurp;
    subject_length = (PCRE2_SIZE)fsta.st_size;

    /*************************************************************************
     * Now we are going to compile the regular expression pattern, and handle *
     * any errors that are detected.                                          *
     *************************************************************************/

    re = pcre2_compile(
            pattern,               /* the pattern */
            PCRE2_ZERO_TERMINATED, /* indicates pattern is zero-terminated */
            0,                     /* default options */
            &errornumber,          /* for error number */
            &erroroffset,          /* for error offset */
            NULL);                 /* use default compile context */

    /* Compilation failed: print the error message and exit. */
    if (re == NULL) {
        PCRE2_UCHAR buffer[256];
        pcre2_get_error_message(errornumber, buffer, sizeof(buffer));
        printf("PCRE2 compilation failed at offset %d: %s\n", (int)erroroffset, buffer);
        exit(EXIT_FAILURE);
    }

    /*************************************************************************
     * If the compilation succeeded, we call PCRE2 again, in order to do a    *
     * pattern match against the subject string. This does just ONE match. If *
     * further matching is needed, it will be done below. Before running the  *
     * match we must set up a match_data block for holding the result. Using  *
     * pcre2_match_data_create_from_pattern() ensures that the block is       *
     * exactly the right size for the number of capturing parentheses in the  *
     * pattern. If you need to know the actual size of a match_data block as  *
     * a number of bytes, you can find it like this:                          *
     *                                                                        *
     * PCRE2_SIZE match_data_size = pcre2_get_match_data_size(match_data);    *
     *************************************************************************/

    match_data = pcre2_match_data_create_from_pattern(re, NULL);

    /* Now run the match. */

    rc = pcre2_match(
            re,                   /* the compiled pattern */
            subject,              /* the subject string */
            subject_length,       /* the length of the subject */
            0,                    /* start at offset 0 in the subject */
            0,                    /* default options */
            match_data,           /* block for storing the result */
            NULL);                /* use default match context */

    /* Matching failed: handle error cases */

    if (rc < 0) {
        switch(rc) {
            case PCRE2_ERROR_NOMATCH: printf("No match\n"); break;
                                      /*
                                         Handle other special cases if you like
                                         */
            default: printf("Matching error %d\n", rc); break;
        }
        pcre2_match_data_free(match_data);   /* Release memory used for the match */
        pcre2_code_free(re);                 /*   data and the compiled pattern. */
        return 1;
    }

    /* Match succeded. Get a pointer to the output vector, where string offsets are
       stored. */

    ovector = pcre2_get_ovector_pointer(match_data);
    printf("Match succeeded at offset %d\n", (int)ovector[0]);


    /*************************************************************************
     * We have found the first match within the subject string. If the output *
     * vector wasn't big enough, say so. Then output any substrings that were *
     * captured.                                                              *
     *************************************************************************/

    /* The output vector wasn't big enough. This should not happen, because we used
       pcre2_match_data_create_from_pattern() above. */

    if (rc == 0)
        printf("ovector was not big enough for all the captured substrings\n");

    /* We must guard against patterns such as /(?=.\K)/ that use \K in an assertion
       to set the start of a match later than its end. In this demonstration program,
       we just detect this case and give up. */

    if (ovector[0] > ovector[1]) {
        printf("\\K was used in an assertion to set the match start after its end.\n"
                "From end to start the match was: %.*s\n", (int)(ovector[0] - ovector[1]),
                (char *)(subject + ovector[1]));
        printf("Run abandoned\n");
        pcre2_match_data_free(match_data);
        pcre2_code_free(re);
        return 1;
    }

    /* Show substrings stored in the output vector by number. Obviously, in a real
       application you might want to do things other than print them. */
    av_c *avc=crea_avc(GBUF);
    int cou=0;
    for (i = 0; i < rc; i++) {
        PCRE2_SPTR substring_start = subject + ovector[2*i];
        PCRE2_SIZE substring_length = ovector[2*i+1] - ovector[2*i];
        // printf("%2d: %.*s\n", i, (int)substring_length, (char *)substring_start);
        switch(i) {
            case 1:
                NATOI(avc->v[avc->vsz].h, (char*)substring_start, (int)substring_length);
                break;
            case 2:
                NATOI(avc->v[avc->vsz].m, (char*)substring_start, (int)substring_length);
                break;
            case 3:
                NATOI(avc->v[avc->vsz].s, (char*)substring_start, (int)substring_length);
                break;
        }
    }
    printf("first timerec: %i:%i:%i\n", avc->v[cou].h, avc->v[cou].m, avc->v[cou].s);
    avc->vsz++;

    /**************************************************************************
     * That concludes the basic part of this demonstration program. We have    *
     * compiled a pattern, and performed a single match. The code that follows *
     * shows first how to access named substrings, and then how to code for    *
     * repeated matches on the same subject.                                   *
     **************************************************************************/

    /* See if there are any named substrings, and if so, show them by name. First
       we have to extract the count of named parentheses from the pattern. */

    (void)pcre2_pattern_info(
            re,                   /* the compiled pattern */
            PCRE2_INFO_NAMECOUNT, /* get the number of named substrings */
            &namecount);          /* where to put the answer */

    if(namecount != 0) {
        printf("No named substrings allowed.\n");
        exit(EXIT_FAILURE);
    }

    /*************************************************************************
     * If the "-g" option was given on the command line, we want to continue  *
     * to search for additional matches in the subject string, in a similar   *
     * way to the /g option in Perl. This turns out to be trickier than you   *
     * might think because of the possibility of matching an empty string.    *
     * What happens is as follows:                                            *
     *                                                                        *
     * If the previous match was NOT for an empty string, we can just start   *
     * the next match at the end of the previous one.                         *
     *                                                                        *
     * If the previous match WAS for an empty string, we can't do that, as it *
     * would lead to an infinite loop. Instead, a call of pcre2_match() is    *
     * made with the PCRE2_NOTEMPTY_ATSTART and PCRE2_ANCHORED flags set. The *
     * first of these tells PCRE2 that an empty string at the start of the    *
     * subject is not a valid match; other possibilities must be tried. The   *
     * second flag restricts PCRE2 to one match attempt at the initial string *
     * position. If this match succeeds, an alternative to the empty string   *
     * match has been found, and we can print it and proceed round the loop,  *
     * advancing by the length of whatever was found. If this match does not  *
     * succeed, we still stay in the loop, advancing by just one character.   *
     * In UTF-8 mode, which can be set by (*UTF) in the pattern, this may be  *
     * more than one byte.                                                    *
     *                                                                        *
     * However, there is a complication concerned with newlines. When the     *
     * newline convention is such that CRLF is a valid newline, we must       *
     * advance by two characters rather than one. The newline convention can  *
     * be set in the regex by (*CR), etc.; if not, we must find the default.  *
     *************************************************************************/

    if (!find_all) {     /* Check for -g */
        pcre2_match_data_free(match_data);  /* Release the memory that was used */
        pcre2_code_free(re);                /* for the match data and the pattern. */
        return 0;                           /* Exit the program. */
    }

    /* Before running the loop, check for UTF-8 and whether CRLF is a valid newline
       sequence. First, find the options with which the regex was compiled and extract
       the UTF state. */
    (void)pcre2_pattern_info(re, PCRE2_INFO_ALLOPTIONS, &option_bits);
    utf8 = (option_bits & PCRE2_UTF) != 0;

    /* Now find the newline convention and see whether CRLF is a valid newline sequence. */
    (void)pcre2_pattern_info(re, PCRE2_INFO_NEWLINE, &newline);
    crlf_is_newline = newline == PCRE2_NEWLINE_ANY ||
        newline == PCRE2_NEWLINE_CRLF ||
        newline == PCRE2_NEWLINE_ANYCRLF;

    /* Loop for second and subsequent matches */
    for (;;) {
        uint32_t options = 0;                   /* Normally no options */
        PCRE2_SIZE start_offset = ovector[1];   /* Start at end of previous match */

        /* If the previous match was for an empty string, we are finished if we are
           at the end of the subject. Otherwise, arrange to run another match at the
           same point to see if a non-empty match can be found. */

        if (ovector[0] == ovector[1]) {
            if (ovector[0] == subject_length) break;
            options = PCRE2_NOTEMPTY_ATSTART | PCRE2_ANCHORED;
        /* If the previous match was not an empty string, there is one tricky case to
           consider. If a pattern contains \K within a lookbehind assertion at the
           start, the end of the matched string can be at the offset where the match
           started. Without special action, this leads to a loop that keeps on matching
           the same substring. We must detect this case and arrange to move the start on
           by one character. The pcre2_get_startchar() function returns the starting
           offset that was passed to pcre2_match(). */
        } else {
            PCRE2_SIZE startchar = pcre2_get_startchar(match_data);
            if (start_offset <= startchar) {
                if (startchar >= subject_length) break;   /* Reached end of subject.   */
                start_offset = startchar + 1;             /* Advance by one character. */
                if (utf8)                                 /* If UTF-8, it may be more  */
                {                                       /*   than one code unit.     */
                    for (; start_offset < subject_length; start_offset++)
                        if ((subject[start_offset] & 0xc0) != 0x80) break;
                }
            }
        }

        /* Run the next matching operation */

        rc = pcre2_match(
                re,                   /* the compiled pattern */
                subject,              /* the subject string */
                subject_length,       /* the length of the subject */
                start_offset,         /* starting offset in the subject */
                options,              /* options */
                match_data,           /* block for storing the result */
                NULL);                /* use default match context */

        /* This time, a result of NOMATCH isn't an error. If the value in "options"
           is zero, it just means we have found all possible matches, so the loop ends.
           Otherwise, it means we have failed to find a non-empty-string match at a
           point where there was a previous empty-string match. In this case, we do what
           Perl does: advance the matching position by one character, and continue. We
           do this by setting the "end of previous match" offset, because that is picked
           up at the top of the loop as the point at which to start again.

           There are two complications: (a) When CRLF is a valid newline sequence, and
           the current position is just before it, advance by an extra byte. (b)
           Otherwise we must ensure that we skip an entire UTF character if we are in
           UTF mode. */

        if (rc == PCRE2_ERROR_NOMATCH) {
            if (options == 0) break;                    /* All matches found */
            ovector[1] = start_offset + 1;              /* Advance one code unit */
            if (crlf_is_newline &&                      /* If CRLF is a newline & */
                    start_offset < subject_length - 1 &&    /* we are at CRLF, */
                    subject[start_offset] == '\r' &&
                    subject[start_offset + 1] == '\n')
                ovector[1] += 1;                          /* Advance by one more. */
            else if (utf8)                              /* Otherwise, ensure we */
            {                                         /* advance a whole UTF-8 */
                while (ovector[1] < subject_length)       /* character. */
                {
                    if ((subject[ovector[1]] & 0xc0) != 0x80) break;
                    ovector[1] += 1;
                }
            }
            continue;    /* Go round the loop again */
        }

        /* Other matching errors are not recoverable. */
        if (rc < 0) {
            printf("Matching error %d\n", rc);
            pcre2_match_data_free(match_data);
            pcre2_code_free(re);
            return 1;
        }

        /* Match succeded */
        printf("\nMatch succeeded again at offset %d\n", (int)ovector[0]);
        condrea_avc(avc); // get ready

        /* The match succeeded, but the output vector wasn't big enough. This
           should not happen. */

        if (rc == 0)
            printf("ovector was not big enough for all the captured substrings\n");

        /* We must guard against patterns such as /(?=.\K)/ that use \K in an
           assertion to set the start of a match later than its end. In this
           demonstration program, we just detect this case and give up. */

        if (ovector[0] > ovector[1]) {
            printf("\\K was used in an assertion to set the match start after its end.\n"
                    "From end to start the match was: %.*s\n", (int)(ovector[0] - ovector[1]),
                    (char *)(subject + ovector[1]));
            printf("Run abandoned\n");
            pcre2_match_data_free(match_data);
            pcre2_code_free(re);
            return 1;
        }

        /* As before, show substrings stored in the output vector by number, and then
           also any named substrings. */
        for (i = 0; i < rc; i++) {
            PCRE2_SPTR substring_start = subject + ovector[2*i];
            size_t substring_length = ovector[2*i+1] - ovector[2*i];
            // printf("%2d: %.*s\n", i, (int)substring_length, (char *)substring_start);
            switch(i) {
                case 1:
                    NATOI(avc->v[avc->vsz].h, (char*)substring_start, (int)substring_length);
                    break;
                case 2:
                    NATOI(avc->v[avc->vsz].m, (char*)substring_start, (int)substring_length);
                    break;
                case 3:
                    NATOI(avc->v[avc->vsz].s, (char*)substring_start, (int)substring_length);
                    break;
            }
        }
        avc->vsz++;

        if(namecount != 0) {
            printf("No named substrings allowed.\n");
            exit(EXIT_FAILURE);
        }
    }      /* End of loop to find second and subsequent matches */

    // avc->vsz=cou;
    norm_avc(avc);
    prtavec(avc);
    pcre2_match_data_free(match_data);
    pcre2_code_free(re);
    free_avc(avc);
    free(fslurp);
    return 0;
}
