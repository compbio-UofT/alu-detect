BOOST=/home/matei/code/get-regions/boost
all:
	BOOST=$(BOOST) make -C src
clean:
	make -C src clean
