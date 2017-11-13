CC = gcc
CFLAGS = -g -Wall -O2 -Wno-unused-variable -Wno-unused-result -Wno-unused-function
#CFLAGS = -g -Wall -O0 -Wno-unused-variable -Wno-unused-result -Wno-unused-function
LIB = -lz -lpthread

BIN_DIR = .
SRC_DIR = ./src

SOURCE = $(wildcard ${SRC_DIR}/*.c) 

OBJS = $(SOURCE:.c=.o)

deMFC : $(OBJS) 
	$(CC) $(CFLAGS) -o deMFC $(OBJS) $(LIB)

clean : 
		rm deMFC $(OBJS) 
