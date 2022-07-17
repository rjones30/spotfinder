#!/bin/bash
cd $(dirname $0)
source setup.sh
if [ "$1" = "-h" -o "$1" = "--help" -o "$1" = "-?" ]; then
    echo "Usage: cobrems_worker.sh [-f]"
    exit 1
elif [ "$1" = "-f" ]; then
    pkill -9 -u $USER celery
    sleep 1
fi
if ps -u $USER | grep -v grep | grep -q celery; then
    echo "already running"
else
    celery -A cobrems_worker worker -l INFO >$(pwd)/cobrems_worker_$(hostname).log 2>&1 & 
fi
