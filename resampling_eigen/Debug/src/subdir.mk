################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/dlt_plane.cpp \
../src/resampling.cpp 

OBJS += \
./src/dlt_plane.o \
./src/resampling.o 

CPP_DEPS += \
./src/dlt_plane.d \
./src/resampling.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I"/home/jpiat/development/SOFTWARE/cv_experiments/resampling_eigen/inc" -O0 -msse4 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


