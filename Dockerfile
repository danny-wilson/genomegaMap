FROM dannywilson/gcat-project
LABEL app="GenomegaMap"
LABEL description="Within-species genome-wide dN/dS estimation from very many genomes"
LABEL maintainer="Daniel Wilson"
LABEL build-type="From source"
ENV MKDIR /tmp/libgenomegaMap
RUN mkdir $MKDIR
COPY . $MKDIR
WORKDIR $MKDIR
RUN make
RUN mv lib* /usr/lib/
RUN mv src/* /usr/include/gcat/
RUN rm *.o
WORKDIR /home/ubuntu
ENTRYPOINT ["/usr/bin/gcat"]
