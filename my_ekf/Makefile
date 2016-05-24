PWD := $(shell pwd)
CFLAGS=-I"${PWD}/inc" -I"${PWD}/../lib/math/inc" -O3 -Wall -fPIC 
LDFLAGS= -lm -lmathd -lmathf -L"${PWD}" -L"${PWD}/../lib/math"

BUILDDIR=build

SRC_FLOAT := src/float_ekf.c

SRC_DOUBLE := src/double_ekf.c

OBJS_FLOAT := src/float_ekf.o

OBJS_DOUBLE := src/double_ekf.o

test_gps_float : test_gps_float.o libekff.so
	gcc -o $@ $< ${LDFLAGS} -lekff
	rm test_gps_float.o

test_gps_double : test_gps_double.o libekf.so
	gcc -o $@  $< ${LDFLAGS} -lekf
	rm test_gps_double.o

test_gps_float.o : test/test_gps.c
	gcc ${CFLAGS} -DFLOAT_EKF -c -o $@ $<

test_gps_double.o : test/test_gps.c
	gcc ${CFLAGS} -DDOUBLE_EKF -c -o $@ $<

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
