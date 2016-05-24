PWD := $(shell pwd)
CFLAGS=-I"${PWD}/inc" -O3 -Wall -fPIC 
LDFLAGS= -lm

BUILDDIR=build

SRC_FLOAT := src/float_mat_ops.c \
	src/mat_types.c \

SRC_DOUBLE := src/double_mat_ops.c \
	src/mat_types.c \

OBJS_FLOAT := src/float_mat_ops.o \
	src/mat_types.o \

OBJS_DOUBLE := src/double_mat_ops.o \
	src/mat_types.o \

all : libmathf.so libmathd.so libmathf.a libmathd.a

libmathf.so : ${OBJS_FLOAT}
	gcc -shared -o $@ ${OBJS_FLOAT} ${LDFLAGS}
	rm ${OBJS_FLOAT}

libmathf.a : ${OBJS_FLOAT}
	gcc -static -o $@ ${OBJS_FLOAT} ${LDFLAGS}

libmathd.so : ${OBJS_DOUBLE}
	gcc -shared -o $@ ${OBJS_DOUBLE} ${LDFLAGS}
	rm ${OBJS_DOUBLE}

libmathd.a : ${OBJS_DOUBLE}
	gcc -static -o $@ ${OBJS_DOUBLE} ${LDFLAGS}

clean :
	rm ${OBJS_FLOAT} ${OBJS_DOUBLE}*.o

*.o : ${SRC_FLOAT} ${SRC_DOUBLE}
	gcc ${CFLAGS} -c -fmessage-length=0 -MMD -MP -o $@ $<
