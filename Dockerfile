FROM rust:1 AS build
RUN apt-get update && apt-get install -y git musl-dev musl-tools libclang-dev llvm-dev clang
RUN rustup update stable && rustup target add x86_64-unknown-linux-musl
WORKDIR /project
COPY . /project
RUN cargo test --release --target=x86_64-unknown-linux-musl
RUN cargo build --release --target=x86_64-unknown-linux-musl

FROM alpine:3.10
RUN apk add --no-cache bash
COPY --from=build /project/target/x86_64-unknown-linux-musl/release/sequencetoolkit /usr/local/bin/