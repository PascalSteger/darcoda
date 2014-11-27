default: all

all: gravimage.py
	cython -3 --embed gravimage.py
	gcc -Wall -Werror -I/usr/include/python3.2mu/ gravimage.c -lpython3.2mu -o gravimage
