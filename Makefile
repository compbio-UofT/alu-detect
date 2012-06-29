all:
	BOOST=$(BOOST) make -C src
distclean: clean
	rm -rf settings data
clean:
	make -C src clean
