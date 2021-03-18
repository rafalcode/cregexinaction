CC=gcc
CFLAGS=-O3
DCFLAGS=-g -Wall -DDBG
SPECLIBS=-lpcre
EXECUTABLES=pcredemo dem2 dem3 dem4 dem5 dem5_d dem6 dem7 pcregrep nx2phy_cheap dem5a dem4a gen0 simpd dem22 dem23

# mm1, an xercise in memory maps ... sys/ctypes-h need to be in place.
mm1: mm1.c
	${CC} ${CFLAGS} -o $@ $^

pcredemo: pcredemo.c
	${CC} ${CFLAGS} -o $@ $^ ${SPECLIBS}

dem2: dem2.c
	${CC} ${CFLAGS} -o $@ $^ ${SPECLIBS}


# so looking for a building brick approach
dem22: dem22.c
	${CC} ${CFLAGS} -o $@ $^ ${SPECLIBS}
dem23: dem23.c
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

# dem5 waht regex's work? Look at a GPX:
# ./dem5 'lon="-*\d+\.\d+">\n\ +<ele>\d+.\d+' ill.gpx "<trkseg>"  "</trkseg"
# note how this will capture a newline and starting spaces


# based on dem5, though this is more customised for the Nexusfile extraction, but either it doesn't workor is not very robust.
dem6: dem6.c
	${CC} ${CFLAGS} -o $@ $^ ${SPECLIBS}

# based on dem6
dem7: dem7.c
	${CC} ${CFLAGS} -o $@ $^ ${SPECLIBS}

# I call this cheap, because it uses pcre and also it should check ntxa and sqlen
nx2phy_cheap: nx2phy_cheap.c
	${CC} ${CFLAGS} -DDEBUG -o $@ $^ ${SPECLIBS}

pcregrep: pcregrep.c
	${CC} ${CFLAGS} -o $@ $^ ${SPECLIBS}


# simpd, copied from dem4a ... _does not_ incorporate btwntxt routines.
simpd: simpd.c
	${CC} ${DCFLAGS} -o $@ $^ ${SPECLIBS}

#dem5a, copied from dem5 ... incorporates btwntxt routines.
# Sample cmdline: ./dem5a "\[\d+\]" aio1.nex "\nMATRIX\n" "\n;\n"
dem5a: dem5a.c
	${CC} ${CFLAGS} -o $@ $^ ${SPECLIBS}

.PHONY: clean

clean:
	rm -f ${EXECUTABLES}
