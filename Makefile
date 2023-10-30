##
# Simple IMPPROVE application container
#
# @file
# @version 0.1

OCI_AGENT = podman

.PHONY: build tarball sif clean all tarinstall uninstall

build:
	$(OCI_AGENT) build -t apply-impprove:latest .

apply-impprove.tar:
	rm -f apply-impprove.tar
	$(OCI_AGENT) save -o apply-impprove.tar localhost/apply-impprove

apply-impprove.tar.gz: apply-impprove.tar
	gzip --keep apply-impprove.tar

apply-impprove.sif: apply-impprove.tar
	apptainer build -F apply-impprove.sif docker-archive:apply-impprove.tar

tarball: apply-impprove.tar apply-impprove.tar.gz

sif: apply-impprove.sif

clean:
	rm -f apply-impprove.tar apply-impprove.tar.gz apply-impprove.sif

all: build tarball sif

tarinstall:
	$(OCI_AGENT) load -i apply-impprove.tar

uninstall:
	$(OCI_AGENT) rmi localhost/apply-impprove

# end
