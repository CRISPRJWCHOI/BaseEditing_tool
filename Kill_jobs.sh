#!/bin/bash

# Confirm the jobs.
# ps aux | grep $USER | grep BaseEdit_freq_ver1.0.py | less

kill -9 $(ps aux | grep $USER | grep BaseEdit_freq_ver1.0.py | awk '{print$2}')
