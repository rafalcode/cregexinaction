/* example file, canonical */
#include <stdio.h>
#include <string.h>
#include <pcre.h>

#define OVECCOUNT 30    /* should be a multiple of 3 */

int main(int argc, char **argv)
{
    /* argument accounting */
    if(argc==1) {
        printf("Usage. SUpply any number of arguments. The final arguemnt is the subject string, all the rest are patterns.\n");
        exit(EXIT_FAILURE);
    }
    const char *error;
    unsigned char *name_table;
    int erroffset;
    int namecount;
    int name_entry_size;
    int ovector[OVECCOUNT]; /* very important, this is for the substring ... the part of the subject which matches the pattern */
    int subject_length;
    int i;

    char *subject = argv[argc-1];
    subject_length = (int)strlen(subject);
    /* Now we are going to compile the regular expression pattern, and handle errors that are detected. */
    int totpatt = argc -2;
    pcre **res=malloc(totpatt*sizeof(pcre*));

    for(i=0;i<totpatt;++i) {
        res[i] = pcre_compile(argv[i+1], 0, &error, &erroffset, NULL);
        if (res[i] == NULL) {
            printf("PCRE compilation failed at offset %d: %s\n", erroffset, error);
        }
    }

    /* If the compilation succeeded, now we execute, which is a pattern match. Only the first one will be done */
    int retx;
    for(i=0;i<totpatt;++i) {
        retx = pcre_exec(res[i], NULL, subject, subject_length, 0, 0, ovector, OVECCOUNT);
        if (retx == PCRE_ERROR_NOMATCH)
            printf("No match at patt %i\n", i+1);
        else
            printf("Yes match at patt %i\n", i+1);
    }

    for(i=0;i<totpatt;++i)
        pcre_free(res[i]);
    free(res);

    return 0;
}
