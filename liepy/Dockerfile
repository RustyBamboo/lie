FROM quay.io/pypa/manylinux2014_x86_64

RUN curl https://sh.rustup.rs -sSf | sh -s -- --default-toolchain stable -y
RUN yum install -y openblas-devel
RUN yum install -y python3
RUN python3 -mpip install maturin==0.12.7_beta.8
RUN echo 'source $HOME/.cargo/env' >> $HOME/.bashrc
ENV PATH="/root/.cargo/bin:${PATH}"

WORKDIR /io
ENTRYPOINT ["/usr/local/bin/maturin"]
