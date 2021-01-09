#!/bin/bash 

sed "s/se path = /se path =${PWD}/" gc_sync_test.conf > gc_synch_test_v1.conf
