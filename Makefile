ARCH=$(shell uname)
ifeq ($(ARCH),Darwin)
	include src/Makefile.osx
else ifeq ($(ARCH),Linux)
	include src/Makefile.linux
else
	$(warning Unsupported architecture. Trying Makefile.linux)
	include src/Makefile.linux
endif
