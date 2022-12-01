TOP_DIR = ../..
DEPLOY_RUNTIME ?= /disks/patric-common/runtime
TARGET ?= /tmp/deployment
include $(TOP_DIR)/tools/Makefile.common

APP_SERVICE = app_service

SERVICE_SPEC = 
SERVICE_NAME = codon_trees
SERVICE_DIR  = $(SERVICE_NAME)
SERVICE_APP_DIR      = $(TARGET)/lib/$(SERVICE_NAME)

SRC_SERVICE_PERL = $(wildcard service-scripts/*.pl)
BIN_SERVICE_PERL = $(addprefix $(BIN_DIR)/,$(basename $(notdir $(SRC_SERVICE_PERL))))
DEPLOY_SERVICE_PERL = $(addprefix $(SERVICE_DIR)/bin/,$(basename $(notdir $(SRC_SERVICE_PERL))))

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

all: bin

bin: $(BIN_PYTHON) $(BIN_PERL) $(BIN_SERVICE_PERL)

dist: 

test: 

deploy: deploy-all

deploy-all: deploy-client

deploy-client: deploy-scripts deploy-libs

deploy-service: deploy-libs deploy-scripts deploy-service-scripts deploy-specs

deploy-specs:
	mkdir -p $(TARGET)/services/$(APP_SERVICE)
	rsync -arv app_specs $(TARGET)/services/$(APP_SERVICE)/.

deploy-service-scripts:
	export KB_TOP=$(TARGET); \
	export KB_RUNTIME=$(DEPLOY_RUNTIME); \
	export KB_PERL_PATH=$(TARGET)/lib ; \
	export KB_PYTHON_PATH=$(TARGET)/lib ; \
	export PATH_ADDITIONS=$(TARGET_VENV)/app-bin; \
	for src in $(SRC_SERVICE_PYTHON) ; do \
	        basefile=`basename $$src`; \
	        base=`basename $$src .py`; \
	        echo install $$src $$base ; \
	        cp $$src $(TARGET)/pybin ; \
	        $(WRAP_PYTHON3_SCRIPT) "$(TARGET)/pybin/$$basefile" $(TARGET)/bin/$$base ; \
	done; \
	for src in $(SRC_SERVICE_PERL) ; do \
	        basefile=`basename $$src`; \
	        base=`basename $$src .pl`; \
	        echo install $$src $$base ; \
	        cp $$src $(TARGET)/plbin ; \
	        $(WRAP_PERL_SCRIPT) "$(TARGET)/plbin/$$basefile" $(TARGET)/bin/$$base ; \
	done



deploy-docs:


build-libs:

$(BIN_DIR)/%: service-scripts/%.pl $(TOP_DIR)/user-env.sh
	$(WRAP_PERL_SCRIPT) '$$KB_TOP/modules/$(CURRENT_DIR)/$<' $@

$(BIN_DIR)/%: service-scripts/%.py $(TOP_DIR)/user-env.sh
	$(WRAP_PYTHON3_SCRIPT) '$$KB_TOP/modules/$(CURRENT_DIR)/$<' $@

include $(TOP_DIR)/tools/Makefile.common.rules
