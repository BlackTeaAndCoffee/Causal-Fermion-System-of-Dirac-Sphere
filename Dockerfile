FROM cfsprebase

ADD ./dist /dist
RUN ls /dist
RUN pip3 install /dist/Numerical_CFS-0.1.9.5.tar.gz 
#RUN pip3 install --upgrade Numerical-CFS

VOLUME /output
VOLUME /usr/local/lib/python3.5/dist-packages/Numerical_CFS/config

#CMD python3 foo.py
#CMD cat  /usr/local/lib/python3.5/dist-packages/Numerical_CFS/config/Settings.cfg
CMD python3 -m Numerical_CFS.Lib_Action_Minimum
