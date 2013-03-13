all:
	BOOST=$(BOOST) make -C src

distclean: clean
	rm -rf settings \
	data/*.fa data/*.fa.fai data/*.bt2 \
	data/*.fa.cat data/*.fa.ref data/*.fa.out data/*.fa.tbl data/*.fa.masked \
	data/alus.* data/deletions.* data/targets.*

clean:
	make -C src clean

VERSION := $(shell cat VERSION)
bin-package: all
	rm -f bin/*.pyc bin/*~
	./get_git_version >GIT_VERSION
	tar cvzf alu-detect-${VERSION}.lx26.x86_64.tar.gz --transform "s,^,alu-detect-${VERSION}/," \
		Makefile README HISTORY VERSION GIT_VERSION get_git_version alu-detect \
		bin/[a-zA-Z0-9]* data/known-novel-alus.*.bed

GIT_VERSION := $(shell ./get_git_version)
test-package: bin-package
	tar xf ~/git/alu-detect/alu-detect-${VERSION}.lx26.x86_64.tar.gz -C ~/opt \
	  --transform "s,alu-detect-${VERSION}/,alu-detect-${VERSION}-${GIT_VERSION}/,"
	rm -rf ~/opt/alu-detect-${VERSION}-${GIT_VERSION}/data
	ln -s ~/git/alu-detect/data ~/opt/alu-detect-${VERSION}-${GIT_VERSION}/data
