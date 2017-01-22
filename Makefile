CC=gcc
CFLAGS=-O3
DCFLAGS=-g -Wall -DDBG
SPECLIBS=-lpcre
EXECUTABLES=pcredemo dem2 dem3 dem4 dem5 dem5_d dem6 dem7 pcregrep nx2phy_cheap dem5a dem4a gen0

# mm1, an xercise in memory maps ... sys/ctypes-h need to be in place.
mm1: mm1.c
	${CC} ${CFLAGS} -o $@ $^

pcredemo: pcredemo.c
	${CC} ${CFLAGS} -o $@ $^ ${SPECLIBS}

dem2: dem2.c
	${CC} ${CFLAGS} -o $@ $^ ${SPECLIBS}

dem3: dem3.c
	${CC} ${CFLAGS} -o $@ $^ ${SPECLIBS}

dem4: dem4.c
	${CC} ${CFLAGS} -o $@ $^ ${SPECLIBS}


# dem5: this parses a KML file. Of course this can be done with libxml, this is just a raw way of doing it.
# best effect is seen on kml route file, you're able to search for the lat/long/alt triplet.
dem5: dem5.c
	${CC} ${CFLAGS} -o $@ $^ ${SPECLIBS}
dem5_d: dem5.c
	${CC} ${DCFLAGS} -o $@ $^ ${SPECLIBS}

# based on dem5, though this is more customised for the Nexusfile extraction, but either it doesn't workor is not very robust.
dem6: dem6.c
	${CC} ${CFLAGS} -o $@ $^ ${SPECLIBS}

# based on dem6
dem7: dem7.c
	${CC} ${CFLAGS} -o $@ $^ ${SPECLIBS}

# I call this cheap, because it uses pcre and also it should check ntxa and sqlen
nx2phy_cheap: nx2phy_cheap.c
	${CC} ${CFLAGS} -DDEBUG -o $@ $^ ${SPECLIBS}

gen0: gen0.c
	${CC} ${CFLAGS} -DDEBUG -o $@ $^ ${SPECLIBS}

pcregrep: pcregrep.c
	${CC} ${CFLAGS} -o $@ $^ ${SPECLIBS}


#dem4a, copied from dem4 ... _does not_ incorporate btwntxt routines.
dem4a: dem4a.c
	${CC} ${CFLAGS} -o $@ $^ ${SPECLIBS}

#dem5a, copied from dem5 ... incorporates btwntxt routines.
# Sample cmdline: ./dem5a "\[\d+\]" aio1.nex "\nMATRIX\n" "\n;\n"
dem5a: dem5a.c
	${CC} ${CFLAGS} -o $@ $^ ${SPECLIBS}

.PHONY: clean

clean:
	rm -f ${EXECUTABLES}
