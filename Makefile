CC = gcc
LIBS = -lm
LDFLAGS = -g
INC = -I./include
CFLAGS = -Wall -O2 -g $(INC)

SRC = $(shell ls src/*.c)
OBJ = $(SRC:.c=.o)
MKS = $(SRC:.c=.d)

TOOLS_SRC = $(shell ls tools/*.c)
TOOLS = $(TOOLS_SRC:.c=)
TOOLS_OBJ = $(TOOLS_SRC:.c=.o)
TOOLS_MKS = $(TOOLS_SRC:.c=.d)

LIBNAME = libad
OUT_STATIC = lib/$(LIBNAME).a
OUT_SHARED = lib/$(LIBNAME).so.1.0.1

.SUFFIXES: .c .o .d

all: $(OUT_STATIC) $(OUT_SHARED) $(TOOLS)

# Create a shared library
$(OUT_SHARED): $(OBJ) $(MKS)
	$(CC) -shared $(LIBS) -Wl,-soname,$(LIBNAME).so.1 -o \
		$(OUT_SHARED) $(OBJ)

# Create a static library
$(OUT_STATIC): $(OBJ) $(MKS)
	ar crs $(OUT_STATIC) $(OBJ)

# Create tools
$(TOOLS): $(TOOLS_OBJ) $(TOOLS_MKS)
	$(CC) $(CFLAGS) $(LIBS) $< -o $@

# Pattern rules
.c.o:
	$(CC) $(CFLAGS) -fPIC -c $< -o $@


# Creates a depend makefile for each source file
#%.d: %.c
.c.d: 
	$(SHELL) -ec '$(CC) -MM $(CFLAGS) $< | sed "s!$*\\.o[ :]*!& $@!g" > $@'

clean:
	rm -f $(OBJ) $(MKS) $(OUT_STATIC) $(OUT_SHARED) \
	$(TOOLS) $(TOOLS_OBJ) $(TOOLS_MKS)

# Include dependency makefiles
include $(MKS)
