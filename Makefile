TARGETS := all clean

.PHONY: $(TARGETS)

$(TARGETS):
	$(MAKE) -C sumatra $@
	$(MAKE) -C sumaclust $@
