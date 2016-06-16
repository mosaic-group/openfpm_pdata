SUBDIRS = src images openfpm_data openfpm_io openfpm_devices openfpm_vcluster openfpm_numerics

bin_PROGRAMS =

test_pdata:
	cd src && make test

test_data:
	cd openfpm_data/src && make test

test_devices:
	cd openfpm_devices/src && make test

test_vcluster:
	cd openfpm_vcluster/src && make test

test_io:
	cd openfpm_io/src && make test

test_numerics:
	cd openfpm_numerics/src && make test

test: test_devices test_data test_vcluster test_pdata test_io test_numerics

.PHONY: test_pdata test_data test_devices test_vcluster test_io test_numerics

