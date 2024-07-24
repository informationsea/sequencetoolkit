FROM rust:1-bookworm AS build
RUN apt-get update && apt-get install -y git libclang-dev llvm-dev clang cmake zlib1g-dev
RUN rustup update stable
WORKDIR /project
COPY . /project
RUN cargo test --release
RUN cargo build --release

FROM debian:bookworm-slim AS download-bcftools
RUN apt-get update && apt-get install -y curl lbzip2 bzip2
ARG BCFTOOLS_VERSION=1.20
RUN curl --fail -OL https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2
RUN tar xf bcftools-${BCFTOOLS_VERSION}.tar.bz2

FROM debian:bookworm-slim AS download-samtools
RUN apt-get update && apt-get install -y curl lbzip2 bzip2
ARG SAMTOOLS_VERSION=1.20
RUN curl --fail -OL https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2
RUN tar xf samtools-${SAMTOOLS_VERSION}.tar.bz2

FROM debian:bookworm-slim AS download-htslib
RUN apt-get update && apt-get install -y curl lbzip2 bzip2
ARG HTSLIB_VERSION=1.20
RUN curl --fail -OL https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2
RUN tar xf htslib-${HTSLIB_VERSION}.tar.bz2

FROM debian:bookworm-slim AS buildenv-bcftools
RUN apt-get update && apt-get install -y build-essential ncurses-dev libbz2-dev zlib1g-dev libcurl4-openssl-dev curl liblzma-dev
ARG BCFTOOLS_VERSION=1.20
COPY --from=download-bcftools /bcftools-${BCFTOOLS_VERSION} /bcftools-${BCFTOOLS_VERSION}
WORKDIR /bcftools-${BCFTOOLS_VERSION}
RUN ./configure --prefix=/usr
RUN make -j4
RUN make install DESTDIR=/dest

FROM debian:bookworm-slim AS buildenv-samtools
RUN apt-get update && apt-get install -y build-essential ncurses-dev libbz2-dev zlib1g-dev libcurl4-openssl-dev curl liblzma-dev
ARG SAMTOOLS_VERSION=1.20
COPY --from=download-samtools /samtools-${SAMTOOLS_VERSION} /bcftools-${SAMTOOLS_VERSION}
WORKDIR /bcftools-${SAMTOOLS_VERSION}
RUN ./configure --prefix=/usr
RUN make -j4
RUN make install DESTDIR=/dest

FROM debian:bookworm-slim AS buildenv-htslib
RUN apt-get update && apt-get install -y build-essential ncurses-dev libbz2-dev zlib1g-dev libcurl4-openssl-dev curl liblzma-dev
ARG HTSLIB_VERSION=1.20
COPY --from=download-htslib /htslib-${HTSLIB_VERSION} /htslib-${HTSLIB_VERSION}
WORKDIR /htslib-${HTSLIB_VERSION}
RUN ./configure --prefix=/usr
RUN make -j4
RUN make install DESTDIR=/dest

FROM debian:bookworm-slim
LABEL maintainer="okamura@informationsea.info"
LABEL org.opencontainers.image.source=https://github.com/informationsea/sequencetoolkit
LABEL org.opencontainers.image.description="Toolkit for Genome Sequence Analysis"
LABEL org.opencontainers.image.licenses=GPL-3.0
RUN apt-get update && apt-get install -y bash libbz2-1.0 libcurl4 liblzma5 && apt-get clean -y && rm -rf /var/lib/apt/lists/*
COPY --from=buildenv-bcftools /dest /
COPY --from=buildenv-samtools /dest /
COPY --from=buildenv-htslib /dest /
COPY --from=build /project/target/release/sequencetoolkit /usr/local/bin/
