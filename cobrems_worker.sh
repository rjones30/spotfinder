#!/bin/bash
cd ~/spotfinder_dev
source setup.sh
pkill -9 celery
celery -A cobrems_worker worker
