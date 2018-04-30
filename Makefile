TOP_DIR = ../..
DEPLOY_RUNTIME ?= /disks/patric-common/runtime
TARGET ?= /tmp/deployment
include $(TOP_DIR)/tools/Makefile.common

SERVICE_SPEC = 
SERVICE_NAME = codon_trees
SERVICE_DIR  = $(SERVICE_NAME)
SERVICE_APP_DIR      = $(TARGET)/lib/$(SERVICE_NAME)

APP_DIR = .
APP_COMPONENTS = 

DATA_API_URL = https://p3.theseed.org/services/data_api
TEMPLATE_DIR = $(shell pwd)/templates
DEPLOY_TEMPLATE_DIR = $(SERVICE_APP_DIR)/templates

TPAGE_BUILD_ARGS =  \
	--define kb_top=$(TARGET) \
	--define kb_runtime=$(DEPLOY_RUNTIME) \
	--define template_dir=$(TEMPLATE_DIR)

TPAGE_DEPLOY_ARGS =  \
	--define kb_top=$(DEPLOY_TARGET) \
	--define kb_runtime=$(DEPLOY_RUNTIME) \
	--define template_dir=$(DEPLOY_TEMPLATE_DIR)

TPAGE_ARGS = \
	--define kb_service_name=$(SERVICE_NAME) \
	--define kb_service_dir=$(SERVICE_DIR) \
	--define kb_service_port=$(SERVICE_PORT) \
	--define kb_psgi=$(SERVICE_PSGI) \
	--define kb_app_dir=$(SERVICE_APP_DIR) \
	--define kb_app_script=$(APP_SCRIPT) \
	--define data_api_url=$(DATA_API_URL)

default: $(BIN_PYTHON)

dist: 

test: 

deploy: deploy-client

deploy-all: deploy-client

deploy-client: deploy-scripts deploy-libs

deploy-service: 

deploy-docs:


build-libs:

include $(TOP_DIR)/tools/Makefile.common.rules
