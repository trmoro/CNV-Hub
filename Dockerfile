FROM public.ecr.aws/lambda/python:3.9

COPY . .
RUN yum -y install wget git tar make

RUN pip3 install -r requirements.txt --target "${LAMBDA_TASK_ROOT}"

#Install R and dependencies
ENV R_VERSION=4.2.0

RUN yum -y install https://dl.fedoraproject.org/pub/epel/epel-release-latest-7.noarch.rpm \
  && wget https://cdn.rstudio.com/r/centos-7/pkgs/R-${R_VERSION}-1-1.x86_64.rpm \
  && yum -y install R-${R_VERSION}-1-1.x86_64.rpm \
  && rm R-${R_VERSION}-1-1.x86_64.rpm

ENV PATH="${PATH}:/opt/R/${R_VERSION}/bin/"

# System requirements for R packages
RUN yum -y install openssl-devel
RUN Rscript -e "install.packages(c('data.table','xgboost'),repos = 'https://packagemanager.rstudio.com/all/__linux__/centos7/latest' )"

RUN git clone https://github.com/kbvstmd/XCNV.git
WORKDIR XCNV
RUN chmod -R 755 . && ./Install.sh
WORKDIR ..

CMD ["cnv_app.handler"]
