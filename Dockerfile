FROM dannywilson/gcat-project
LABEL app="GCAT with omegaMap library"
LABEL description="General computational analysis tool with omegaMap library"
LABEL maintainer="Daniel Wilson"
LABEL build-type="From source"
ENV MKDIR /tmp/libgcat_omegaMap
RUN mkdir $MKDIR
COPY . $MKDIR
WORKDIR $MKDIR
RUN make
RUN mv lib* /usr/lib/
RUN mv src/* /usr/include/gcat/
RUN rm *.o
WORKDIR /home/ubuntu
ENTRYPOINT ["/usr/bin/gcat"]