FROM ubuntu:20.04
ENV TZ=Europe/Berlin
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
RUN apt-get update ; \
	apt-get install -y vim wget rustc cargo git cmake sudo

RUN useradd --create-home --shell /bin/bash --groups sudo appuser; \
    echo "appuser:1" | chpasswd

USER appuser

RUN mkdir /home/appuser/src ; \
    cd /home/appuser/src ; \ 
    git clone https://github.com/lskatz/fasten; \
    cd fasten && cargo build --release 

RUN cd /home/appuser/src; \
	git clone https://github.com/rust-bio/rust-htslib.git; \
	cd rust-htslib 
RUN cd /home/appuser/src; \
	git clone https://github.com/tzeitim/rogtk.git ;\
	cd rogtk

RUN echo 'export PATH="$PATH:/home/appuser/src/fasten/target/release/"' >> /home/appuser/.bashrc
	
WORKDIR /home/appuser
CMD ["/bin/bash"]
ENTRYPOINT bash
