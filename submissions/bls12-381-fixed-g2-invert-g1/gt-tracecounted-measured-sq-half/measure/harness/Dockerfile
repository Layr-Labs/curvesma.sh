FROM python:3.12-slim@sha256:a6e34c598f2467ed0e9a8d349809fcd8b5c603269512df273a0bb1784edc11b1
RUN apt-get update && apt-get install -y --no-install-recommends gcc libc6-dev valgrind perl make && rm -rf /var/lib/apt/lists/*
WORKDIR /work
