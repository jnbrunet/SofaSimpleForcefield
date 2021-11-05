FROM dolfinx/dolfinx:latest
ENV DEBIAN_FRONTEND noninteractive
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN apt-get -qq update && apt-get -qq upgrade

RUN apt-get -qq --no-install-recommends install libpython3.8 python3.8 python3-pip libpython3.8-dev \
&&  apt-get -qq --no-install-recommends install qtbase5-dev libqt5charts5-dev libqt5opengl5-dev qtwebengine5-dev libopengl0 libeigen3-dev libglew-dev zlib1g-dev \
&&  apt-get -qq --no-install-recommends install libboost1.71-dev libboost-system1.71-dev libboost-filesystem1.71-dev libboost-program-options1.71-dev libboost-thread1.71-dev \
&&  apt-get -qq --no-install-recommends install g++ ccache ninja-build curl unzip git \
&&  python3.8 -m pip install numpy \
&&  apt-get clean

# Install CMake
ADD https://github.com/Kitware/CMake/releases/download/v3.20.1/cmake-3.20.1-linux-x86_64.sh /tmp
RUN chmod a+x /tmp/cmake-3.20.1-linux-x86_64.sh \
&&  /tmp/cmake-3.20.1-linux-x86_64.sh --skip-license --prefix=/usr/local \
&&  rm /tmp/cmake-3.20.1-linux-x86_64.sh

WORKDIR /opt
ENV SOFA_ROOT=/opt/sofa