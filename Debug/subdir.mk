################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../allocation_mw.cpp \
../convex_hull.cpp \
../envelope.cpp \
../global_problem.cpp \
../instance.cpp \
../main.cpp \
../subproblem.cpp 

OBJS += \
./allocation_mw.o \
./convex_hull.o \
./envelope.o \
./global_problem.o \
./instance.o \
./main.o \
./subproblem.o 

CPP_DEPS += \
./allocation_mw.d \
./convex_hull.d \
./envelope.d \
./global_problem.d \
./instance.d \
./main.d \
./subproblem.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


