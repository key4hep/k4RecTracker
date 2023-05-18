make:
	@ mkdir -p build install ; \
	cd build ; \
	source /cvmfs/sw.hsf.org/key4hep/setup.sh ; \
	cmake .. -DCMAKE_INSTALL_PREFIX=../install ; \
	make install -j4 ; \
	cd .. ; \
	export LD_LIBRARY_PATH=${CURDIR}/install/lib:${CURDIR}/install/lib64:$$LD_LIBRARY_PATH ; \
	export PYTHONPATH=${CURDIR}/install/python:$$PYTHONPATH ; \
	printf "#!/bin/bash\nsource /cvmfs/sw.hsf.org/key4hep/setup.sh\nexport LD_LIBRARY_PATH=${CURDIR}/install/lib:${CURDIR}/install/lib64:\$$LD_LIBRARY_PATH\nexport PYTHONPATH=${CURDIR}/install/python:\$$PYTHONPATH\n" > ${CURDIR}/setup.sh ; \
	chmod +x ${CURDIR}/setup.sh
.PHONY: clean
clean:
	@ (rm -r build install && rm setup.sh) || true
