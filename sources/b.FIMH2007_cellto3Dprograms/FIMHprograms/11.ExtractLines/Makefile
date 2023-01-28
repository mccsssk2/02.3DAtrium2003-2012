# Simple Makefile for the Sun 02 / avnode 02 Sun v 9a arch.
#
# C compiler and options
#
CC=    /opt/SUNWspro/bin/cc -xO3
LIB=	-lm

#
# Object files
#
OBJ=	extract_lines.o
#
# Compile
#
extract	:$(OBJ)
	$(CC) -o $@ $(OBJ) $(LIB)

.c.o:
	$(CC) -c $<

#
# Clean out object files and the executable.
#
clean:
	rm *.o extract
