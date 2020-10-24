# #######################################################################
# ## BUILDTIME_IMAGE                                                   ##
# #######################################################################

FROM ubuntu:20.04 AS BUILDTIME_IMAGE

RUN apt update \
    && apt -y install \
    g++ \
    cmake \
    libboost-dev \
    libboost-program-options-dev \
    libarmadillo-dev \
    libgtest-dev \
    && apt clean

# #######################################################################
# ## RUNTIME_IMAGE                                                     ##
# #######################################################################

FROM ubuntu:20.04 AS RUNTIME_IMAGE

RUN apt update \
    && apt -y install \
    libgomp1 \
    libboost-program-options1.71.0 \
    libarmadillo9 \
    && apt clean

# #######################################################################
# ## BUILD_IMAGE                                                       ##
# #######################################################################

FROM BUILDTIME_IMAGE AS BUILD_IMAGE

COPY . /bfpt
RUN which gcc
RUN gcc --version
RUN mkdir /bfpt/build_release_from_dockerfile
RUN cd /bfpt/build_release_from_dockerfile && cmake -DCMAKE_BUILD_TYPE=Release ..
RUN cd /bfpt/build_release_from_dockerfile && make "-j$(nproc)"

# #######################################################################
# ## PRODUCT_IMAGE                                                     ##
# #######################################################################

FROM RUNTIME_IMAGE AS PRODUCT_IMAGE

COPY --from=BUILD_IMAGE /bfpt/build_release_from_dockerfile/bin/monostar_app /bfpt/bin/
ENTRYPOINT ["/bfpt/bin/monostar_app"]