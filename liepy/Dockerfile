FROM konstin2/maturin

WORKDIR /
RUN yum makecache
RUN yum -y install openblas-devel # for system installs
RUN yum -y install gcc-gfortran # for static installs

WORKDIR /io
ENTRYPOINT ["/usr/bin/maturin"]
