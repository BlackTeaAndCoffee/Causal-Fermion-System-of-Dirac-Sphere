# Pull base image.
FROM ubuntu:16.04
RUN apt-get update
RUN apt-get dist-upgrade -y
 
RUN apt-get -y dist-upgrade
RUN apt-get -y install python-software-properties
RUN apt-get -y install software-properties-common
RUN apt-get update
#ADD /home/mustafa/Regensburg/Causal-Fermion-System-of-Dirac-Sphere/dist /dist
#RUN ls /dist
#RUN pip3 install /dist/Numerical_CFS-0.1.9.5.tar.gz 
#RUN pip3 install --upgrade Numerical-CFS
  
RUN apt-get -y install python3-dev



RUN add-apt-repository ppa:symengine/ppa
RUN apt-get update
RUN apt-get -y install libsymengine-dev


RUN apt-get -y install python3-pip
#der nächste Teil sollte bald raus, ich habs nur jetzt hier weil cache veraltet ist.
RUN pip3 install --upgrade setuptools 
#
COPY ./requirements_for_docker.txt / 

RUN pip3 install -r ./requirements_for_docker.txt 
