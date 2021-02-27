FROM nanjiang/common-ubuntu
LABEL maintainer "Nanjiang Shu (nanjiang.shu@scilifelab.se)"
LABEL version "1.4"

#================================
# Install basics
#===============================
RUN apt-get update -y &&\
    apt-get install -y apt-utils  \
                       curl wget bc \
                       build-essential  \
                       make  \
                       locales-all \
                       blast2

#================================
#  Add predzinc source code
#===============================
WORKDIR /app
# add the source code to WORKDIR /app
ADD predzinc ./predzinc

RUN /scratch/ /static/
# building predzinc
RUN cd /app/predzinc/&& \
    make && make install

#================================
# Setting ENVs
#===============================
ENV USER_DIRS "/app"
ENV BLASTDB /data/blastdb
ENV BLASTBIN /usr/bin

CMD ["/bin/bash" ]
