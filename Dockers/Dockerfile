FROM cfsprebase

ADD ./dist /dist
ADD ./build /build
RUN ls /dist
RUN pip3 install /dist/Numerical_CFS_Docker-0.1.9.8.1.tar.gz 
#RUN pip3 install --upgrade Numerical-CFS

VOLUME /output
VOLUME /usr/local/lib/python3.5/dist-packages/Numerical_CFS_Docker/config

#CMD python3 foo.py
#CMD cat  /usr/local/lib/python3.5/dist-packages/Numerical_CFS/config/Settings.cfg
CMD python3 -m Numerical_CFS_Docker.Lib_Action_Minimum
