FROM dongli/metmpi-tcp:0.0.1

LABEL description="This is an implementation of barotropic model using quadratic-conservation finite difference scheme."
LABEL version="0.0.1"
LABEL maintainer="Li Dong <dongli@lasg.iap.ac.cn>"

RUN cd && git clone https://github.com/dongli/fortran-datetime && \
    mkdir fortran-datetime/build && cd fortran-datetime/build && \
    cmake .. -DCMAKE_BUILD_TYPE=Release && make && sudo make install && \
    cd && rm -rf fortran-datetime
RUN cd && git clone https://github.com/dongli/qconswm && \
    mkdir qconswm/build && cd qconswm/build && \
    cmake .. -DCMAKE_BUILD_TYPE=Release && make

WORKDIR qconswm/build
