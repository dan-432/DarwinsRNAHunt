# Force x86_64 platform
FROM --platform=linux/amd64 ubuntu:22.04

# Install system dependencies
RUN apt-get update && apt-get install -y \
    wget \
    make \
    gcc \
    g++ \
    perl \
    && rm -rf /var/lib/apt/lists/*

# Download and install CMfinder
WORKDIR /opt

RUN wget -O cmfinder-0.4.1.18.tgz https://sourceforge.net/projects/weinberg-cmfinder/files/cmfinder-0.4.1.18.tgz \
    && tar -xzvf cmfinder-0.4.1.18.tgz \
    && cd cmfinder-0.4.1.18 \
    && ./configure --build=x86_64-unknown-linux-gnu CFLAGS="-O3 -fcommon" CXXFLAGS="-O3 -fcommon" \
    && make

# Set CMfinder environment variable and add to PATH
ENV CMfinder=/opt/cmfinder-0.4.1.18
ENV PATH="/opt/cmfinder-0.4.1.18/bin:/opt/cmfinder-0.4.1.18/cmfinder03:/opt/cmfinder-0.4.1.18/cmfinder04:${PATH}"

WORKDIR /data
# No ENTRYPOINT - you specify the command each time
CMD ["/bin/bash"]