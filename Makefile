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

test-package: bin-package
	cd /tmp
	tar xf ~/git/alu-detect/alu-detect-${VERSION}.lx26.x86_64.tar.gz
	pwd
	ls -l
	GIT_VERSION := $(shell cat alu-detect-${VERSION}/GIT_VERSION)
	mv alu-detect-${VERSION} alu-detect-${VERSION}-${GIT_VERSION}
	mv alu-detect-${VERSION}-${GIT_VERSION} ~/opt
	cd ~/opt/alu-detect-${VERSION}-${GIT_VERSION}
	mv data data-orig
	ln -s ~/git/alu-detect/data
