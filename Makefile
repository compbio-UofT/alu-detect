BOOST=/filer/harun/boost_1_49_0/boost
all:
	BOOST=$(BOOST) make -C src
clean:
	make -C src clean
