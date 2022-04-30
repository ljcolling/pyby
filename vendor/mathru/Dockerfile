FROM archlinux:latest

RUN pwd
RUN pacman -Sy

RUN pacman -S rustup --noconfirm
RUN pacman -S binutils --noconfirm
RUN pacman -S gcc --noconfirm
RUN pacman -S cargo --noconfirm
RUN pacman -S cmake make pkgconf freetype2 python python-pip --noconfirm

RUN rustup default nightly

RUN rustup component add llvm-tools-preview

RUN cargo install grcov

RUN export RUSTFLAGS="-Zinstrument-coverage"


