FROM python:3.6.5
MAINTAINER Charles Tapley Hoyt "cthoyt@gmail.com"

RUN pip3 install --upgrade pip
RUN pip3 install gunicorn mysqlclient

COPY . /app
WORKDIR /app

RUN pip3 install .[web]
RUN bio2bel_hgnc populate --skip-hcop
