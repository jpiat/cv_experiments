PWD := $(shell pwd)
CFLAGS=-I"${PWD}/inc" -O3 -Wall -fPIC 
LDFLAGS= -lm

BUILDDIR=build

SRC_FLOAT := src/float_ekf.c \
	src/float_mat_ops.c \
	src/mat_types.c \

SRC_DOUBLE := src/double_ekf.c \
	src/double_mat_ops.c \
	src/mat_types.c \

OBJS_FLOAT := src/float_ekf.o \
	src/float_mat_ops.o \
	src/mat_types.o \

OBJS_DOUBLE := src/double_ekf.o \
	src/double_mat_ops.o \
	src/mat_types.o \

test_gps_float : test_gps_float.o libekff.so
	gcc -o $@ -L${PWD} $< -lekff -lm
	rm test_gps_float.o

test_gps_double : test_gps_double.o libekf.so
	gcc -o $@ -L${PWD} $< -lekf -lm
	rm test_gps_double.o

test_gps_float.o : test/test_gps.c
	gcc -I"inc/" -DFLOAT_EKF -c -o $@ $<

test_gps_double.o : test/test_gps.c
	gcc -I"inc/" -DDOUBLE_EKF -c -o $@ $<

libekff.so : ${OBJS_FLOAT}
	gcc -shared -o $@ ${OBJS_FLOAT} ${LDFLAGS}
	rm ${OBJS_FLOAT}

libekff.a : ${OBJS_FLOAT}	
	gcc -static -o $@ ${OBJS_FLOAT} ${LDFLAGS}

libekf.so : ${OBJS_DOUBLE}
	gcc -shared -o $@ ${OBJS_DOUBLE} ${LDFLAGS}
	rm ${OBJS_DOUBLE}

libekf.a : ${OBJS_DOUBLE}
	gcc -static -o $@ ${OBJS_DOUBLE} ${LDFLAGS}

clean :
	rm ${OBJS_FLOAT} ${OBJS_DOUBLE}*.o

*.o : ${SRC_FLOAT} ${SRC_DOUBLE}
	gcc ${CFLAGS} -c -fmessage-length=0 -MMD -MP -o $@ $<
