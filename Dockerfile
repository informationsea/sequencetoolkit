FROM rust:1 AS build
RUN apt-get update && apt-get install -y git musl-dev musl-tools libclang-dev llvm-dev clang
RUN rustup update stable && rustup target add x86_64-unknown-linux-musl
WORKDIR /project
COPY . /project
RUN cargo test --release --target=x86_64-unknown-linux-musl
RUN cargo build --release --target=x86_64-unknown-linux-musl

FROM alpine:3.10 AS download-bcftools
RUN apk add curl libarchive-tools
ENV BCFTOOLS_VERSION=1.10.2
RUN curl -OL https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2
RUN bsdtar xf bcftools-${BCFTOOLS_VERSION}.tar.bz2

FROM alpine:3.10 AS download-htslib
RUN apk add curl libarchive-tools
ENV BCFTOOLS_VERSION=1.10.2
RUN curl -OL https://github.com/samtools/htslib/releases/download/${BCFTOOLS_VERSION}/htslib-${BCFTOOLS_VERSION}.tar.bz2
RUN bsdtar xf htslib-${BCFTOOLS_VERSION}.tar.bz2

FROM alpine:3.10 AS buildenv-bcftools
RUN apk add gcc make libc-dev ncurses-dev bzip2-dev zlib-dev curl-dev curl xz-dev
ENV BCFTOOLS_VERSION=1.10.2
COPY --from=download-bcftools /bcftools-${BCFTOOLS_VERSION} /bcftools-${BCFTOOLS_VERSION}
WORKDIR /bcftools-${BCFTOOLS_VERSION}
RUN ./configure --prefix=/usr
RUN make -j4
RUN make install DESTDIR=/dest

FROM alpine:3.10 AS buildenv-htslib
RUN apk add gcc make libc-dev ncurses-dev bzip2-dev zlib-dev curl-dev curl xz-dev
ENV BCFTOOLS_VERSION=1.10.2
COPY --from=download-htslib /htslib-${BCFTOOLS_VERSION} /htslib-${BCFTOOLS_VERSION}
WORKDIR /htslib-${BCFTOOLS_VERSION}
RUN ./configure --prefix=/usr
RUN make -j4
RUN make install DESTDIR=/dest

FROM alpine:3.10
RUN apk add --no-cache bash libbz2 libcurl xz-libs 
COPY --from=buildenv-bcftools /dest /
COPY --from=buildenv-htslib /dest /
COPY --from=build /project/target/x86_64-unknown-linux-musl/release/sequencetoolkit /usr/local/bin/
