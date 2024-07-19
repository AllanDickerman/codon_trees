TOP_DIR = ../..
DEPLOY_RUNTIME ?= /disks/patric-common/runtime
TARGET ?= /tmp/deployment
include $(TOP_DIR)/tools/Makefile.common

APP_SERVICE = app_service

SERVICE_SPEC = 
SERVICE_NAME = codon_trees
SERVICE_DIR  = $(SERVICE_NAME)
SERVICE_APP_DIR      = $(TARGET)/lib/$(SERVICE_NAME)

WRAP_PYTHON_TOOL = wrap_python3
WRAP_PYTHON_SCRIPT = bash $(TOOLS_DIR)/$(WRAP_PYTHON3_TOOL).sh

all: bin

bin: $(BIN_PYTHON) $(BIN_PERL) $(BIN_SERVICE_PERL)

dist: 

test: 

deploy: deploy-all

deploy-all: deploy-client

deploy-client: deploy-scripts deploy-libs

deploy-service: deploy-libs deploy-scripts deploy-service-scripts deploy-specs

deploy-docs:


build-libs:

include $(TOP_DIR)/tools/Makefile.common.rules
