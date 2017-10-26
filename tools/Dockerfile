FROM dongli/metenv
MAINTAINER Li Dong <dongli@lasg.iap.ac.cn>

RUN cd && git clone https://github.com/dongli/fortran-datetime && \
    mkdir fortran-datetime/build && cd fortran-datetime/build && \
    cmake .. -DCMAKE_BUILD_TYPE=Release && make && make install && \
    cd && rm -rf fortran-datetime
RUN cd && git clone https://github.com/dongli/qconswm && \
    mkdir qconswm/build && cd qconswm/build && \
    cmake .. -DCMAKE_BUILD_TYPE=Release && make

WORKDIR /root/qconswm/build
