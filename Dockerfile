FROM rust:1-buster AS build
RUN apt-get update && apt-get install -y git libclang-dev llvm-dev clang cmake zlib1g-dev
RUN rustup update stable
WORKDIR /project
COPY . /project
RUN cargo test --release
RUN cargo build --release

FROM debian:buster-slim AS download-bcftools
RUN apt-get update && apt-get install -y curl lbzip2 bzip2
ARG BCFTOOLS_VERSION=1.15.1
RUN curl -OL https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2
RUN tar xf bcftools-${BCFTOOLS_VERSION}.tar.bz2

FROM debian:buster-slim AS download-samtools
RUN apt-get update && apt-get install -y curl lbzip2 bzip2
ARG BCFTOOLS_VERSION=1.15.1
RUN curl -OL https://github.com/samtools/samtools/releases/download/${BCFTOOLS_VERSION}/samtools-${BCFTOOLS_VERSION}.tar.bz2
RUN tar xf samtools-${BCFTOOLS_VERSION}.tar.bz2

FROM debian:buster-slim AS download-htslib
RUN apt-get update && apt-get install -y curl lbzip2 bzip2
ARG BCFTOOLS_VERSION=1.15.1
RUN curl -OL https://github.com/samtools/htslib/releases/download/${BCFTOOLS_VERSION}/htslib-${BCFTOOLS_VERSION}.tar.bz2
RUN tar xf htslib-${BCFTOOLS_VERSION}.tar.bz2

FROM debian:buster-slim AS buildenv-bcftools
RUN apt-get update && apt-get install -y build-essential ncurses-dev libbz2-dev zlib1g-dev libcurl4-openssl-dev curl liblzma-dev
ARG BCFTOOLS_VERSION=1.15.1
COPY --from=download-bcftools /bcftools-${BCFTOOLS_VERSION} /bcftools-${BCFTOOLS_VERSION}
WORKDIR /bcftools-${BCFTOOLS_VERSION}
RUN ./configure --prefix=/usr
RUN make -j4
RUN make install DESTDIR=/dest

FROM debian:buster-slim AS buildenv-samtools
RUN apt-get update && apt-get install -y build-essential ncurses-dev libbz2-dev zlib1g-dev libcurl4-openssl-dev curl liblzma-dev
ARG BCFTOOLS_VERSION=1.15.1
COPY --from=download-samtools /samtools-${BCFTOOLS_VERSION} /bcftools-${BCFTOOLS_VERSION}
WORKDIR /bcftools-${BCFTOOLS_VERSION}
RUN ./configure --prefix=/usr
RUN make -j4
RUN make install DESTDIR=/dest

FROM debian:buster-slim AS buildenv-htslib
RUN apt-get update && apt-get install -y build-essential ncurses-dev libbz2-dev zlib1g-dev libcurl4-openssl-dev curl liblzma-dev
ARG BCFTOOLS_VERSION=1.15.1
COPY --from=download-htslib /htslib-${BCFTOOLS_VERSION} /htslib-${BCFTOOLS_VERSION}
WORKDIR /htslib-${BCFTOOLS_VERSION}
RUN ./configure --prefix=/usr
RUN make -j4
RUN make install DESTDIR=/dest

FROM debian:buster-slim
RUN apt-get update && apt-get install -y bash libbz2-1.0 libcurl4 liblzma5 && apt-get clean -y && rm -rf /var/lib/apt/lists/*
COPY --from=buildenv-bcftools /dest /
COPY --from=buildenv-samtools /dest /
COPY --from=buildenv-htslib /dest /
COPY --from=build /project/target/release/sequencetoolkit /usr/local/bin/
