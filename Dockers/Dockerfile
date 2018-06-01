# Pull base image.
FROM ubuntu:16.04
RUN apt-get update
RUN apt-get dist-upgrade -y

RUN apt-get -y dist-upgrade
RUN apt-get -y install python-software-properties
RUN apt-get -y install software-properties-common
RUN apt-get update

RUN apt-get -y install python3-dev



RUN add-apt-repository ppa:symengine/ppa
RUN apt-get update
RUN apt-get -y install libsymengine-dev


RUN apt-get -y install python3-pip
COPY ./requirements_for_docker.txt / 
COPY . /

RUN pip3 install -r ./requirements_for_docker.txt 

