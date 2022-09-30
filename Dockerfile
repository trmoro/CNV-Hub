FROM ubuntu:20.04
ADD . /
RUN apt-get update -y
RUN apt-get -y install python python3 python3-pip
RUN pip3 install -r requirements.txt
RUN DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get -y install tzdata
RUN apt-get -y install dirmngr gnupg apt-transport-https ca-certificates software-properties-common
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'
RUN apt-get -y install r-base
RUN R -e "install.packages('data.table')"
RUN R -e "install.packages('xgboost')"
ENTRYPOINT ["sh","/start.sh"]
