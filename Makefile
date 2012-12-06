all:
	BOOST=$(BOOST) make -C src
distclean: clean
	rm -rf settings \
	data/*fa data/*fa.gz data/*fai data/*bt2 \
	data/alus.* data/deletions.* data/targets.*
clean:
	make -C src clean
