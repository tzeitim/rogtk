FROM ubuntu:20.04
RUN apt-get update
RUN apt-get install -y vim
RUN apt-get install -y wget
RUN apt-get install -y rustc
RUN apt-get install -y cargo
RUN apt-get install -y git
RUN mkdir src
RUN cd src && git clone https://github.com/lskatz/fasten
RUN cd src/fasten && cargo build --release 
RUN echo hello
RUN rustc --version
CMD ["/bin/bash"]
ENTRYPOINT bash
