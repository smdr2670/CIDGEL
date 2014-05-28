# ------------------------------------------------
# Makefile for the software which enumrates all 
# or degree compatible groebner bases of a code 
# ideal. (software-name not available yet)
#
# Author: daniel.rembold@tu-harburg.de
# Date  : 2014-05-27
#
# TODO: make the include path work
#
# Changelog :
#   2014-05-27 - first version
# ------------------------------------------------

# project name (generate executable with this name)
TARGET   = cidgel

CC       = gcc
# compiling flags here
CFLAGS   = -std=c99 -Wall -I.
MATH	 = -lm

LINKER   = gcc -o
# linking flags here
LFLAGS   = -Wall -lm -I.

# change these to set the proper directories where each files should be
SRCDIR   = src
INCDIR	 = inc
OBJDIR   = obj
BINDIR   = build

SOURCES  := $(wildcard $(SRCDIR)/*.c)
INCLUDES := $(wildcard $(INCDIR)/*.h)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.c=$(OBJDIR)/%.o)
rm       = rm -f


$(TARGET): $(OBJECTS)
	@$(LINKER) $@ $(LFLAGS) $(OBJECTS) $(MATH)
	@echo "Linking complete!"

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.c
	@$(CC) $(CFLAGS) -c $< -o $@ 
	@echo "Compiled "$<" successfully!"

.PHONEY: clean
clean:
	@$(rm) $(OBJECTS)
	@echo "Cleanup complete!"

.PHONEY: remove
remove: clean
	@$(rm) $(BINDIR)/$(TARGET)
	@echo "Executable removed!"



